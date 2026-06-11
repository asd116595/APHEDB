#define PROFILE
#include "aggregation.h"

using namespace lbcrypto;

// ===== SUM =====
// Implemented manually with shift-and-add to handle both even and odd counts.
// Algorithm:
//   For step = 1, 2, 4, 8, ... while step < numSlots:
//     shifted = rotate(ctxt, -step)
//     ctxt = ctxt + shifted
//   After log2(numSlots) steps, every slot contains the total sum.
//
// For odd numbers (e.g. 5), after step=1: (a0+a1, a1+a2, a2+a3, a3+a4, a4+a0)
//   After step=2: (a0+a1+a2+a3, a1+a2+a3+a4, ...)
//   After step=4: all slots = a0+a1+a2+a3+a4 ✓
//
// The algorithm naturally handles odd counts because the wrap-around rotation
// accumulates all elements correctly (n rotations of step size cover all elements
// when step is a power of 2 and step < n).
Ciphertext<DCRTPoly> HomSum(
    const std::vector<Poly>& ctxtData,
    CryptoContextCKKSRNS::ContextType& cc,
    BigInteger QBFVInit,
    BigInteger Q,
    BigInteger Bigq,
    KeyPair<DCRTPoly>& keyPair,
    uint32_t numSlotsCKKS,
    uint32_t numSlots,
    uint32_t depth,
    uint32_t levelsAvailableBeforeBootstrap) {

    // Convert RLWE to CKKS
    auto ctxtCKKS = SchemeletRLWEMP::ConvertRLWEToCKKS(
        *cc, ctxtData, keyPair.publicKey, QBFVInit, numSlotsCKKS,
        depth - (levelsAvailableBeforeBootstrap > 0));

    // Manual shift-and-add sum
    // For n slots, need ceil(log2(n)) steps
    // After step k, each slot contains sum of 2^k consecutive elements
    // After all steps, each slot = total sum of all n elements
    auto result = ctxtCKKS;

    for (uint32_t step = 1; step < numSlots; step <<= 1) {
        auto rotated = cc->EvalRotate(result, -static_cast<int32_t>(step));
        result = cc->EvalAdd(result, rotated);
    }

    return result;
}


// ===== COUNT =====
Ciphertext<DCRTPoly> HomCount(
    const std::vector<Poly>& ctxtData,
    CryptoContextCKKSRNS::ContextType& cc,
    BigInteger QBFVInit,
    BigInteger PInput,
    BigInteger POutput,
    BigInteger Q,
    BigInteger Bigq,
    uint64_t scaleTHI,
    size_t order,
    uint32_t numSlots,
    uint32_t ringDim,
    uint32_t dcrtBits,
    std::shared_ptr<M4DCRTParams>& ep,
    KeyPair<DCRTPoly>& keyPair,
    uint32_t numSlotsCKKS,
    uint32_t depth,
    uint32_t levelsAvailableBeforeBootstrap,
    uint32_t levelsComputation) {

    // Make a copy of ctxtData for EqualBootstrapping (takes non-const ref)
    std::vector<Poly> ctxtCopy = ctxtData;

    // EqualBootstrapping checks: is value == 0? Returns 1 if equal to 0.
    // For COUNT, we want 1 if value != 0.
    auto eqZero = EqualBootstrapping(ctxtCopy, cc, QBFVInit, PInput, POutput, Q, Bigq,
                                      scaleTHI, order, numSlots, ringDim, dcrtBits,
                                      ep, keyPair, numSlotsCKKS, depth,
                                      levelsAvailableBeforeBootstrap, levelsComputation);

    // Compute 1 - eqZero to get "is non-zero"
    std::vector<double> x_one(numSlots, 1);
    uint32_t resultLevel = eqZero->GetLevel();
    Plaintext ptxtone = cc->MakeCKKSPackedPlaintext(x_one, 1, resultLevel, nullptr, numSlotsCKKS);
    auto cone = cc->Encrypt(keyPair.publicKey, ptxtone);
    auto nonZeroMask = cc->EvalSub(cone, eqZero);

    // Manual shift-and-add sum to count non-zero entries
    // Same approach as HomSum: log2(numSlots) steps of rotate + add
    auto countResult = nonZeroMask;
    for (uint32_t step = 1; step < numSlots; step <<= 1) {
        auto rotated = cc->EvalRotate(countResult, -static_cast<int32_t>(step));
        countResult = cc->EvalAdd(countResult, rotated);
    }

    return countResult;
}


// ===== AVG =====
std::pair<Ciphertext<DCRTPoly>, Ciphertext<DCRTPoly>> HomAvg(
    const std::vector<Poly>& ctxtData,
    CryptoContextCKKSRNS::ContextType& cc,
    BigInteger QBFVInit,
    BigInteger PInput,
    BigInteger POutput,
    BigInteger Q,
    BigInteger Bigq,
    uint64_t scaleTHI,
    size_t order,
    uint32_t numSlots,
    uint32_t ringDim,
    uint32_t dcrtBits,
    std::shared_ptr<M4DCRTParams>& ep,
    KeyPair<DCRTPoly>& keyPair,
    uint32_t numSlotsCKKS,
    uint32_t depth,
    uint32_t levelsAvailableBeforeBootstrap,
    uint32_t levelsComputation) {

    auto sumResult = HomSum(ctxtData, cc, QBFVInit, Q, Bigq, keyPair,
                            numSlotsCKKS, numSlots, depth, levelsAvailableBeforeBootstrap);

    auto countResult = HomCount(ctxtData, cc, QBFVInit, PInput, POutput, Q, Bigq,
                                scaleTHI, order, numSlots, ringDim, dcrtBits,
                                ep, keyPair, numSlotsCKKS, depth,
                                levelsAvailableBeforeBootstrap, levelsComputation);

    return {sumResult, countResult};
}


// ===== MIN/MAX: tournament reduction =====
// Strategy:
//   1. Convert initial RLWE data to CKKS once
//   2. For each tournament step, extract pairs via EvalRotate in CKKS
//   3. For each pair, convert both to RLWE, call GreaterThan/LessThan (returns CKKS),
//      apply selector in CKKS: result = sel*(a-b) + b
//   4. Repeat until all slots contain the min/max
//
// This works for any number of values <= numSlots.
// For ringDim=4096, you can have up to 4096 values per ciphertext.

namespace {

/**
 * @brief Compare two CKKS ciphertexts and select max element-wise
 *
 * Converts CKKS -> RLWE for comparison (using GreaterThan),
 * gets CKKS selector, applies in CKKS.
 */
Ciphertext<DCRTPoly> MaxOfTwoCKKS(
    Ciphertext<DCRTPoly> a,        // by value - we may ModReduce
    Ciphertext<DCRTPoly> b,        // by value - we may ModReduce
    CryptoContextCKKSRNS::ContextType& cc,
    BigInteger QBFVInit,
    BigInteger PInput,
    BigInteger PDigit,
    BigInteger Q,
    BigInteger Bigq,
    uint64_t scaleTHI,
    uint64_t scaleStepTHI,
    size_t order,
    uint32_t numSlots,
    uint32_t ringDim,
    uint32_t dcrtBits,
    std::shared_ptr<M4DCRTParams>& ep,
    KeyPair<DCRTPoly>& keyPair,
    uint32_t numSlotsCKKS,
    uint32_t depth,
    uint32_t levelsAvailableBeforeBootstrap,
    uint32_t levelsComputation) {

    // Convert CKKS a, b to RLWE for comparison
    auto aRLWE = SchemeletRLWEMP::ConvertCKKSToRLWE(a, Q);
    auto bRLWE = SchemeletRLWEMP::ConvertCKKSToRLWE(b, Q);

    // sel = GreaterThan(a, b) -> CKKS ciphertext, 1 where a > b
    auto sel = GreaterThan(aRLWE, bRLWE, cc, QBFVInit, PInput, PDigit, Q, Bigq,
                           scaleTHI, scaleStepTHI, order, numSlots, ringDim, dcrtBits,
                           ep, keyPair, numSlotsCKKS, depth,
                           levelsAvailableBeforeBootstrap, levelsComputation);

    // Ensure a and b are at compatible levels with sel (non-const copies)
    while (a->GetLevel() > sel->GetLevel()) {
        cc->ModReduceInPlace(a);
    }
    while (b->GetLevel() > sel->GetLevel()) {
        cc->ModReduceInPlace(b);
    }

    // max = sel * a + (1 - sel) * b
    // = sel * (a - b) + b
    auto diff = cc->EvalSub(a, b);
    auto selTimesDiff = cc->EvalMult(sel, diff);
    cc->ModReduceInPlace(selTimesDiff);
    auto maxResult = cc->EvalAdd(selTimesDiff, b);

    return maxResult;
}

/**
 * @brief Compare two CKKS ciphertexts and select min element-wise
 */
Ciphertext<DCRTPoly> MinOfTwoCKKS(
    Ciphertext<DCRTPoly> a,        // by value
    Ciphertext<DCRTPoly> b,        // by value
    CryptoContextCKKSRNS::ContextType& cc,
    BigInteger QBFVInit,
    BigInteger PInput,
    BigInteger PDigit,
    BigInteger Q,
    BigInteger Bigq,
    uint64_t scaleTHI,
    uint64_t scaleStepTHI,
    size_t order,
    uint32_t numSlots,
    uint32_t ringDim,
    uint32_t dcrtBits,
    std::shared_ptr<M4DCRTParams>& ep,
    KeyPair<DCRTPoly>& keyPair,
    uint32_t numSlotsCKKS,
    uint32_t depth,
    uint32_t levelsAvailableBeforeBootstrap,
    uint32_t levelsComputation) {

    // Convert CKKS a, b to RLWE for comparison
    auto aRLWE = SchemeletRLWEMP::ConvertCKKSToRLWE(a, Q);
    auto bRLWE = SchemeletRLWEMP::ConvertCKKSToRLWE(b, Q);

    // sel = LessThan(a, b) -> CKKS ciphertext, 1 where a < b
    auto sel = LessThan(aRLWE, bRLWE, cc, QBFVInit, PInput, PDigit, Q, Bigq,
                        scaleTHI, scaleStepTHI, order, numSlots, ringDim, dcrtBits,
                        ep, keyPair, numSlotsCKKS, depth,
                        levelsAvailableBeforeBootstrap, levelsComputation);

    // Ensure a and b are at compatible levels with sel (non-const copies)
    while (a->GetLevel() > sel->GetLevel()) {
        cc->ModReduceInPlace(a);
    }
    while (b->GetLevel() > sel->GetLevel()) {
        cc->ModReduceInPlace(b);
    }

    // min = sel * a + (1 - sel) * b = sel * (a - b) + b
    auto diff = cc->EvalSub(a, b);
    auto selTimesDiff = cc->EvalMult(sel, diff);
    cc->ModReduceInPlace(selTimesDiff);
    auto minResult = cc->EvalAdd(selTimesDiff, b);

    return minResult;
}

} // anonymous namespace


// ===== MAX =====
Ciphertext<DCRTPoly> HomMax(
    const std::vector<Poly>& ctxtData,
    CryptoContextCKKSRNS::ContextType& cc,
    BigInteger QBFVInit,
    BigInteger PInput,
    BigInteger PDigit,
    BigInteger Q,
    BigInteger Bigq,
    uint64_t scaleTHI,
    uint64_t scaleStepTHI,
    size_t order,
    uint32_t numSlots,
    uint32_t ringDim,
    uint32_t dcrtBits,
    std::shared_ptr<M4DCRTParams>& ep,
    KeyPair<DCRTPoly>& keyPair,
    uint32_t numSlotsCKKS,
    uint32_t depth,
    uint32_t levelsAvailableBeforeBootstrap,
    uint32_t levelsComputation) {

    std::cerr << "\n=== HomMax: tournament on " << numSlots << " slots ===" << std::endl;

    // Step 1: Convert RLWE data to CKKS once
    auto ctxtCKKS = SchemeletRLWEMP::ConvertRLWEToCKKS(
        *cc, ctxtData, keyPair.publicKey, QBFVInit, numSlotsCKKS,
        depth - (levelsAvailableBeforeBootstrap > 0));

    // Step 2: Tournament reduction
    // log2(numSlots) rounds, each round reduces by half
    uint32_t n = numSlots;
    auto current = ctxtCKKS;

    for (uint32_t step = 1; step < n; step *= 2) {
        std::cerr << "  Tournament step=" << step << std::endl;

        // Rotate current by -step to align adjacent pairs
        auto rotated = cc->EvalRotate(current, -static_cast<int32_t>(step));

        // Compute max of current and rotated element-wise
        current = MaxOfTwoCKKS(current, rotated, cc, QBFVInit, PInput, PDigit, Q, Bigq,
                               scaleTHI, scaleStepTHI, order, numSlots, ringDim, dcrtBits,
                               ep, keyPair, numSlotsCKKS, depth,
                               levelsAvailableBeforeBootstrap, levelsComputation);
    }

    // After log2(n) steps, all slots contain the maximum value
    return current;
}


// ===== MIN =====
Ciphertext<DCRTPoly> HomMin(
    const std::vector<Poly>& ctxtData,
    CryptoContextCKKSRNS::ContextType& cc,
    BigInteger QBFVInit,
    BigInteger PInput,
    BigInteger PDigit,
    BigInteger Q,
    BigInteger Bigq,
    uint64_t scaleTHI,
    uint64_t scaleStepTHI,
    size_t order,
    uint32_t numSlots,
    uint32_t ringDim,
    uint32_t dcrtBits,
    std::shared_ptr<M4DCRTParams>& ep,
    KeyPair<DCRTPoly>& keyPair,
    uint32_t numSlotsCKKS,
    uint32_t depth,
    uint32_t levelsAvailableBeforeBootstrap,
    uint32_t levelsComputation) {

    std::cerr << "\n=== HomMin: tournament on " << numSlots << " slots ===" << std::endl;

    // Convert RLWE data to CKKS once
    auto ctxtCKKS = SchemeletRLWEMP::ConvertRLWEToCKKS(
        *cc, ctxtData, keyPair.publicKey, QBFVInit, numSlotsCKKS,
        depth - (levelsAvailableBeforeBootstrap > 0));

    uint32_t n = numSlots;
    auto current = ctxtCKKS;

    for (uint32_t step = 1; step < n; step *= 2) {
        std::cerr << "  Tournament step=" << step << std::endl;

        auto rotated = cc->EvalRotate(current, -static_cast<int32_t>(step));

        current = MinOfTwoCKKS(current, rotated, cc, QBFVInit, PInput, PDigit, Q, Bigq,
                               scaleTHI, scaleStepTHI, order, numSlots, ringDim, dcrtBits,
                               ep, keyPair, numSlotsCKKS, depth,
                               levelsAvailableBeforeBootstrap, levelsComputation);
    }

    return current;
}
