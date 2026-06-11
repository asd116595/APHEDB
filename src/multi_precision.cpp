#define PROFILE
#include "multi_precision.h"

using namespace lbcrypto;

// ===== EqualBootstrapping: checks if value == 0 =====
// Uses the equal-LUT: func(x) = 1 if (x % PInput) == 0, else 0
Ciphertext<DCRTPoly> EqualBootstrapping(
    std::vector<Poly> ctxtBFV,
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

    auto a = PInput.ConvertToInt<int64_t>();
    auto func = [a](int64_t x) -> int64_t {
        return (x % a) == 0 ? 1 : 0;
    };

    std::vector<int64_t> coeffint;
    std::vector<std::complex<double>> coeffcomp;
    bool binaryLUT = (PInput.ConvertToInt() == 2) && (order == 1);

    if (binaryLUT) {
        coeffint = {func(1), func(0) - func(1), 0};
    } else {
        coeffcomp = GetHermiteTrigCoefficients(func, PInput.ConvertToInt(), order, scaleTHI);
    }

    SchemeletRLWEMP::ModSwitch(ctxtBFV, Q, QBFVInit);

    auto encryptedDigit = ctxtBFV;
    encryptedDigit[0].SwitchModulus(Bigq, 1, 0, 0);
    encryptedDigit[1].SwitchModulus(Bigq, 1, 0, 0);

    auto ctxt = SchemeletRLWEMP::ConvertRLWEToCKKS(
        *cc, encryptedDigit, keyPair.publicKey, Bigq, numSlotsCKKS,
        depth - (levelsAvailableBeforeBootstrap > 0));

    Ciphertext<DCRTPoly> ctxtAfterFBT;
    if (binaryLUT)
        ctxtAfterFBT = cc->EvalFBTNoDecoding(ctxt, coeffint, PInput.GetMSB() - 1, ep->GetModulus(), order);
    else
        ctxtAfterFBT = cc->EvalFBTNoDecoding(ctxt, coeffcomp, PInput.GetMSB() - 1, ep->GetModulus(), order);

    return cc->EvalHomDecoding(ctxtAfterFBT, scaleTHI, levelsComputation);
}

// Local digit-sign LUT for multi-precision comparison only.
// diff = X - Y mod PInput, so X > Y iff diff is in [1, PInput / 2).
static Ciphertext<DCRTPoly> DigitGreaterBootstrapping(
    std::vector<Poly> ctxtBFV,
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

    (void)numSlots;
    (void)ringDim;
    (void)dcrtBits;

    const auto positiveMod = [](int64_t x, int64_t m) -> int64_t {
        int64_t r = x % m;
        return (r < 0) ? (r + m) : r;
    };

    auto p = PInput.ConvertToInt<int64_t>();
    auto func = [positiveMod, p](int64_t x) -> int64_t {
        int64_t r = positiveMod(x, p);
        return (r > 0 && r < (p / 2)) ? 1 : 0;
    };

    std::vector<int64_t> coeffint;
    std::vector<std::complex<double>> coeffcomp;
    bool binaryLUT = (PInput.ConvertToInt() == 2) && (order == 1);

    if (binaryLUT) {
        coeffint = {func(1), func(0) - func(1), 0};
    } else {
        coeffcomp = GetHermiteTrigCoefficients(func, PInput.ConvertToInt(), order, scaleTHI);
    }

    std::vector<uint32_t> lvlb = {1, 1};
    uint32_t lvlsAfterBoot = 5;
    if (binaryLUT)
        cc->EvalFBTSetup(coeffint, numSlotsCKKS, PInput, POutput, Bigq, keyPair.publicKey, {0, 0}, lvlb,
                         lvlsAfterBoot, levelsComputation, order);
    else
        cc->EvalFBTSetup(coeffcomp, numSlotsCKKS, PInput, POutput, Bigq, keyPair.publicKey, {0, 0}, lvlb,
                         lvlsAfterBoot, levelsComputation, order);

    SchemeletRLWEMP::ModSwitch(ctxtBFV, Q, QBFVInit);

    auto encryptedDigit = ctxtBFV;
    encryptedDigit[0].SwitchModulus(Bigq, 1, 0, 0);
    encryptedDigit[1].SwitchModulus(Bigq, 1, 0, 0);

    auto ctxt = SchemeletRLWEMP::ConvertRLWEToCKKS(
        *cc, encryptedDigit, keyPair.publicKey, Bigq, numSlotsCKKS,
        depth - (levelsAvailableBeforeBootstrap > 0));

    Ciphertext<DCRTPoly> ctxtAfterFBT;
    if (binaryLUT)
        ctxtAfterFBT = cc->EvalFBTNoDecoding(ctxt, coeffint, PInput.GetMSB() - 1, ep->GetModulus(), order);
    else
        ctxtAfterFBT = cc->EvalFBTNoDecoding(ctxt, coeffcomp, PInput.GetMSB() - 1, ep->GetModulus(), order);

    return cc->EvalHomDecoding(ctxtAfterFBT, scaleTHI, levelsComputation);
}

static void AlignCiphertextLevels(Ciphertext<DCRTPoly>& a, Ciphertext<DCRTPoly>& b,
                                  CryptoContextCKKSRNS::ContextType& cc) {
    while (a->GetLevel() < b->GetLevel()) {
        cc->ModReduceInPlace(a);
    }
    while (b->GetLevel() < a->GetLevel()) {
        cc->ModReduceInPlace(b);
    }
}

// ===== MultiPrecisionSign: core recursive comparison =====
// Uses the Multi_Precision_Sign framework:
//   result = equal * prev_result + (1 - equal) * sign
// Process from LSB to MSB
// No inverse sign (Negatediff SignCipher) is computed.
// Equal is computed using EqualBootstrapping.
Ciphertext<DCRTPoly> MultiPrecisionSign(
    const std::vector<std::vector<Poly>>& ctxtX,
    const std::vector<std::vector<Poly>>& ctxtY,
    uint32_t numDigits,
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

    auto start = std::chrono::high_resolution_clock::now();

    std::cerr << "\n=== MultiPrecisionSign: " << numDigits << " digits ===" << std::endl;
    (void)PDigit;
    (void)scaleStepTHI;

    // --- Step 1: For each digit, compute sign(diff) and equal(diff) ---
    std::vector<Ciphertext<DCRTPoly>> sign_block(numDigits);
    std::vector<Ciphertext<DCRTPoly>> equal_block(numDigits);

    for (uint32_t i = 0; i < numDigits; ++i) {
        std::cerr << "  Processing digit " << i << " ("
                  << (i == numDigits - 1 ? "LSB" : (i == 0 ? "MSB" : "mid"))
                  << ")" << std::endl;

        // diff = X[i] - Y[i]
        std::vector<Poly> diff = {
            ctxtX[i][0] - ctxtY[i][0],
            ctxtX[i][1] - ctxtY[i][1]
        };

        // sign(diff): 1 when X[i] > Y[i]
        BigInteger Q_copy = Q;
        BigInteger PInput_copy = PInput;
        sign_block[i] = DigitGreaterBootstrapping(diff, cc, QBFVInit, PInput_copy, PInput_copy, Q_copy, Bigq,
                                                  scaleTHI, order, numSlots, ringDim, dcrtBits,
                                                  ep, keyPair, numSlotsCKKS, depth,
                                                  levelsAvailableBeforeBootstrap, levelsComputation);

        // equal(diff): 1 when X[i] == Y[i], using EqualBootstrapping
        // Recompute diff since SignCipher modifies it
        std::vector<Poly> diffForEqual = {
            ctxtX[i][0] - ctxtY[i][0],
            ctxtX[i][1] - ctxtY[i][1]
        };
        BigInteger Q_copy2 = Q;
        BigInteger PInput_copy2 = PInput;
        equal_block[i] = EqualBootstrapping(diffForEqual, cc, QBFVInit, PInput_copy2, PInput_copy2, Q_copy2, Bigq,
                                            scaleTHI, order, numSlots, ringDim, dcrtBits,
                                            ep, keyPair, numSlotsCKKS, depth,
                                            levelsAvailableBeforeBootstrap, levelsComputation);

        std::cerr << "  digit[" << i << "] sign_block level=" << sign_block[i]->GetLevel()
                  << ", equal_block level=" << equal_block[i]->GetLevel() << std::endl;
    }

    // --- Step 2: Combine results from LSB to MSB ---
    // Recursion: result[i] = sign_i + equal_i * result[i+1]
    // - If digits are NOT equal (equal=0): result = sign (current digit decides)
    // - If digits ARE equal (equal=1): result = result[i+1] (look at lower digits)
    // LSB (i = numDigits-1): result = sign (no lower digits to look at)

    Ciphertext<DCRTPoly> result;
    for (int i = (int)numDigits - 1; i >= 0; --i) {
        if (i == (int)numDigits - 1) {
            // LSB: result = sign (equal * 0 = 0, so sign + 0 = sign)
            result = sign_block[i];
        } else {
            // result = sign + equal * prev_result
            AlignCiphertextLevels(equal_block[i], result, cc);
            auto term = cc->EvalMult(equal_block[i], result);
            AlignCiphertextLevels(sign_block[i], term, cc);
            result = cc->EvalAdd(sign_block[i], term);
        }

        // Debug: print levels
        std::cerr << "  digit[" << i << "] sign_block level=" << sign_block[i]->GetLevel()
                  << ", equal_block level=" << equal_block[i]->GetLevel()
                  << ", result level=" << result->GetLevel() << std::endl;
    }

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> elapsed = end - start;
    std::cout << "MultiPrecisionSign time: " << elapsed.count() << " ms" << std::endl;

    // ---- Debug: decrypt result via RLWE to verify plaintext ----
    auto resultDecoded = cc->EvalHomDecoding(result, scaleTHI, levelsComputation);
    std::vector<Poly> ctxtResultRLWE = SchemeletRLWEMP::ConvertCKKSToRLWE(resultDecoded, Q);
    auto plainResult = SchemeletRLWEMP::DecryptCoeff(ctxtResultRLWE, Q, PInput, keyPair.secretKey, ep, numSlotsCKKS, numSlots);
    std::cerr << "  [DEBUG] MultiPrecisionSign plaintext (first " << std::min((size_t)8, plainResult.size()) << " elements): [";
    for (size_t s = 0; s < std::min((size_t)8, plainResult.size()); s++) {
        std::cerr << plainResult[s] << (s < std::min((size_t)8, plainResult.size()) - 1 ? ", " : "");
    }
    std::cerr << "]" << std::endl;
    return result;
}

// ===== MultiPrecisionEqual: all digits equal =====
Ciphertext<DCRTPoly> MultiPrecisionEqual(
    const std::vector<std::vector<Poly>>& ctxtX,
    const std::vector<std::vector<Poly>>& ctxtY,
    uint32_t numDigits,
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

    auto start = std::chrono::high_resolution_clock::now();
    std::cerr << "\n=== MultiPrecisionEqual: " << numDigits << " digits ===" << std::endl;

    // equal = equal_0 * equal_1 * ... * equal_{k-1}
    Ciphertext<DCRTPoly> result;

    for (uint32_t i = 0; i < numDigits; i++) {
        std::vector<Poly> diff = {
            ctxtX[i][0] - ctxtY[i][0],
            ctxtX[i][1] - ctxtY[i][1]
        };

        BigInteger Q_copy = Q;
        BigInteger PInput_copy = PInput;

        auto equal_i = EqualBootstrapping(diff, cc, QBFVInit, PInput_copy, PInput_copy, Q_copy, Bigq,
                                          scaleTHI, order, numSlots, ringDim, dcrtBits,
                                          ep, keyPair, numSlotsCKKS, depth,
                                          levelsAvailableBeforeBootstrap, levelsComputation);

        if (i == 0) {
            result = equal_i;
        } else {
            AlignCiphertextLevels(result, equal_i, cc);
            result = cc->EvalMult(result, equal_i);
        }
    }

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> elapsed = end - start;
    std::cout << "MultiPrecisionEqual time: " << elapsed.count() << " ms" << std::endl;

    

    return result;
}

// ===== MultiPrecisionGreaterThan: X > Y =====
Ciphertext<DCRTPoly> MultiPrecisionGreaterThan(
    const std::vector<std::vector<Poly>>& ctxtX,
    const std::vector<std::vector<Poly>>& ctxtY,
    uint32_t numDigits,
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

    std::cerr << "\n=== MultiPrecisionGreaterThan (X > Y) ===" << std::endl;
    return MultiPrecisionSign(ctxtX, ctxtY, numDigits, cc, QBFVInit, PInput, PDigit, Q, Bigq,
                              scaleTHI, scaleStepTHI, order, numSlots, ringDim, dcrtBits,
                              ep, keyPair, numSlotsCKKS, depth,
                              levelsAvailableBeforeBootstrap, levelsComputation);
}

// ===== MultiPrecisionGreaterEqual: X >= Y = 1 - (Y > X) =====
Ciphertext<DCRTPoly> MultiPrecisionGreaterEqual(
    const std::vector<std::vector<Poly>>& ctxtX,
    const std::vector<std::vector<Poly>>& ctxtY,
    uint32_t numDigits,
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

    std::cerr << "\n=== MultiPrecisionGreaterEqual (X >= Y) ===" << std::endl;

    // X >= Y  <=>  NOT (Y > X)
    auto y_gt_x = MultiPrecisionSign(ctxtY, ctxtX, numDigits, cc, QBFVInit, PInput, PDigit, Q, Bigq,
                                     scaleTHI, scaleStepTHI, order, numSlots, ringDim, dcrtBits,
                                     ep, keyPair, numSlotsCKKS, depth,
                                     levelsAvailableBeforeBootstrap, levelsComputation);

    std::vector<double> x_one = {1};
    if (x_one.size() < numSlots)
        x_one = Fill<double>(x_one, numSlots);

    Plaintext ptxt_one = cc->MakeCKKSPackedPlaintext(x_one, 1, y_gt_x->GetLevel(), nullptr, numSlotsCKKS);
    return cc->EvalSub(ptxt_one, y_gt_x);
}

// ===== MultiPrecisionLessThan: X < Y = Y > X =====
Ciphertext<DCRTPoly> MultiPrecisionLessThan(
    const std::vector<std::vector<Poly>>& ctxtX,
    const std::vector<std::vector<Poly>>& ctxtY,
    uint32_t numDigits,
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

    std::cerr << "\n=== MultiPrecisionLessThan (X < Y) ===" << std::endl;
    return MultiPrecisionSign(ctxtY, ctxtX, numDigits, cc, QBFVInit, PInput, PDigit, Q, Bigq,
                              scaleTHI, scaleStepTHI, order, numSlots, ringDim, dcrtBits,
                              ep, keyPair, numSlotsCKKS, depth,
                              levelsAvailableBeforeBootstrap, levelsComputation);
}

// ===== MultiPrecisionLessEqual: X <= Y = 1 - (X > Y) =====
Ciphertext<DCRTPoly> MultiPrecisionLessEqual(
    const std::vector<std::vector<Poly>>& ctxtX,
    const std::vector<std::vector<Poly>>& ctxtY,
    uint32_t numDigits,
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

    std::cerr << "\n=== MultiPrecisionLessEqual (X <= Y) ===" << std::endl;

    // X <= Y  <=>  NOT (X > Y)
    auto x_gt_y = MultiPrecisionSign(ctxtX, ctxtY, numDigits, cc, QBFVInit, PInput, PDigit, Q, Bigq,
                                     scaleTHI, scaleStepTHI, order, numSlots, ringDim, dcrtBits,
                                     ep, keyPair, numSlotsCKKS, depth,
                                     levelsAvailableBeforeBootstrap, levelsComputation);

    std::vector<double> x_one = {1};
    if (x_one.size() < numSlots)
        x_one = Fill<double>(x_one, numSlots);

    Plaintext ptxt_one = cc->MakeCKKSPackedPlaintext(x_one, 1, x_gt_y->GetLevel(), nullptr, numSlotsCKKS);
    return cc->EvalSub(ptxt_one, x_gt_y);
}
