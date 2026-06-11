#define PROFILE
#include "group_equal.h"

#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <map>
#include <set>

using namespace lbcrypto;

/**
 * @brief Helper: decrypt CKKS ciphertext and print
 */
void DecryptAndPrintCKKS(
    const lbcrypto::Ciphertext<lbcrypto::DCRTPoly>& ctxt,
    lbcrypto::CryptoContextCKKSRNS::ContextType& cc,
    const BigInteger& Q,
    uint64_t scaleTHI,
    uint32_t levelsComputation,
    const BigInteger& POutput,
    lbcrypto::KeyPair<lbcrypto::DCRTPoly>& keyPair,
    std::shared_ptr<lbcrypto::M4DCRTParams>& ep,
    uint32_t numSlotsCKKS,
    uint32_t numSlots,
    const std::string& label) {

    auto ctxtDecoded = cc->EvalHomDecoding(ctxt, scaleTHI, levelsComputation);
    auto polys = SchemeletRLWEMP::ConvertCKKSToRLWE(ctxtDecoded, Q);
    auto computed = SchemeletRLWEMP::DecryptCoeff(polys, Q, POutput, keyPair.secretKey, ep, numSlotsCKKS, numSlots);

    std::cerr << "  " << label << ": [";
    for (size_t i = 0; i < numSlots; ++i) {
        std::cerr << computed[i] << " ";
    }
    std::cerr << "]" << std::endl;
}

/**
 * @brief Helper: decrypt RLWE and print
 */
void DecryptAndPrintRLWE(
    const lbcrypto::Ciphertext<lbcrypto::DCRTPoly>& ctxt,
    lbcrypto::CryptoContextCKKSRNS::ContextType& cc,
    const BigInteger& Q,
    uint64_t scaleTHI,
    uint32_t levelsComputation,
    const BigInteger& POutput,
    lbcrypto::KeyPair<lbcrypto::DCRTPoly>& keyPair,
    std::shared_ptr<lbcrypto::M4DCRTParams>& ep,
    uint32_t numSlotsCKKS,
    uint32_t numSlots,
    const std::string& label) {

    auto ctxtDecoded = cc->EvalHomDecoding(ctxt, scaleTHI, levelsComputation);
    auto polys = SchemeletRLWEMP::ConvertCKKSToRLWE(ctxtDecoded, Q);
    auto computed = SchemeletRLWEMP::DecryptCoeff(polys, Q, POutput, keyPair.secretKey, ep, numSlotsCKKS, numSlots);

    std::cerr << "  " << label << ": [";
    for (size_t i = 0; i < numSlots; ++i) {
        std::cerr << computed[i] << " ";
    }
    std::cerr << "]" << std::endl;
}


/**
 * @brief Print expected equal results for given data and groups
 */
void PrintExpectedGroupEqual(
    const std::vector<int64_t>& data,
    const std::vector<int64_t>& groups,
    uint32_t numSlots,
    const std::string& label) {

    std::map<int64_t, int64_t> groupRef;  // group -> reference value
    for (size_t i = 0; i < groups.size() && i < numSlots; ++i) {
        if (groupRef.find(groups[i]) == groupRef.end()) {
            groupRef[groups[i]] = data[i];
        }
    }

    std::cerr << "  " << label << ": [";
    for (size_t i = 0; i < numSlots; ++i) {
        if (i < groups.size()) {
            int64_t g = groups[i];
            int64_t ref = groupRef[g];
            std::cerr << (data[i] == ref ? 1 : 0) << " ";
        } else {
            std::cerr << "0 ";
        }
    }
    std::cerr << "]" << std::endl;
}


/**
 * @brief Run a single GroupEqual test
 */
void RunGroupEqualTest(
    lbcrypto::CryptoContextCKKSRNS::ContextType& cc,
    lbcrypto::KeyPair<lbcrypto::DCRTPoly>& keyPair,
    std::shared_ptr<lbcrypto::M4DCRTParams>& ep,
    const BigInteger& PInput,
    const BigInteger& POutput,
    const BigInteger& Q,
    const BigInteger& Bigq,
    uint64_t scaleTHI,
    size_t order,
    uint32_t numSlots,
    uint32_t ringDim,
    uint32_t dcrtBits,
    uint32_t numSlotsCKKS,
    uint32_t depth,
    uint32_t levelsAvailableBeforeBootstrap,
    uint32_t levelsComputation,
    const std::vector<int64_t>& data,
    const std::vector<int64_t>& groups,
    const std::string& testName) {

    std::cerr << "\n===== " << testName << " =====" << std::endl;
    std::cerr << "  Data:   " << data << std::endl;
    std::cerr << "  Groups: " << groups << std::endl;

    // Prepare data
    auto x = data;
    if (x.size() < numSlots) x = Fill<int64_t>(x, numSlots);

    // Encrypt
    auto ctxtData = SchemeletRLWEMP::EncryptCoeff(x, QBFVINIT, PInput, keyPair.secretKey, ep);

    // Run GroupEqual
    auto start = std::chrono::high_resolution_clock::now();
    auto result = GroupEqual(ctxtData, groups, cc, QBFVINIT, PInput, POutput, Q, Bigq,
                              scaleTHI, order, numSlots, ringDim, dcrtBits,
                              ep, keyPair, numSlotsCKKS, depth,
                              levelsAvailableBeforeBootstrap, levelsComputation);
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> elapsed = end - start;

    // Print result
    DecryptAndPrintRLWE(result, cc, Q, scaleTHI, levelsComputation, BigInteger(2), keyPair, ep,
                        numSlotsCKKS, numSlots, "Computed (0/1)");

    PrintExpectedGroupEqual(data, groups, numSlots, "Expected (0/1)");

    std::cout << "  Time: " << elapsed.count() << " ms" << std::endl;
}


/**
 * @brief Run a GroupEqualMultiGroup test
 */
void RunMultiGroupTest(
    lbcrypto::CryptoContextCKKSRNS::ContextType& cc,
    lbcrypto::KeyPair<lbcrypto::DCRTPoly>& keyPair,
    std::shared_ptr<lbcrypto::M4DCRTParams>& ep,
    const BigInteger& PInput,
    const BigInteger& POutput,
    const BigInteger& Q,
    const BigInteger& Bigq,
    uint64_t scaleTHI,
    size_t order,
    uint32_t numSlots,
    uint32_t ringDim,
    uint32_t dcrtBits,
    uint32_t numSlotsCKKS,
    uint32_t depth,
    uint32_t levelsAvailableBeforeBootstrap,
    uint32_t levelsComputation,
    const std::vector<int64_t>& data,
    const std::vector<int64_t>& groups,
    int64_t numGroups,
    const std::string& testName) {

    std::cerr << "\n===== " << testName << " =====" << std::endl;
    std::cerr << "  Data:   " << data << std::endl;
    std::cerr << "  Groups: " << groups << std::endl;
    std::cerr << "  NumGroups: " << numGroups << std::endl;

    // Prepare data
    auto x = data;
    if (x.size() < numSlots) x = Fill<int64_t>(x, numSlots);

    // Encrypt
    auto ctxtData = SchemeletRLWEMP::EncryptCoeff(x, QBFVINIT, PInput, keyPair.secretKey, ep);

    // Get reference values per group
    std::map<int64_t, int64_t> groupRef;
    for (size_t i = 0; i < groups.size(); ++i) {
        if (groupRef.find(groups[i]) == groupRef.end()) {
            groupRef[groups[i]] = data[i];
        }
    }

    // Run GroupEqualMultiGroup
    auto start = std::chrono::high_resolution_clock::now();
    auto results = GroupEqualMultiGroup(ctxtData, groups, numGroups, cc, QBFVINIT, PInput, POutput, Q, Bigq,
                                         scaleTHI, order, numSlots, ringDim, dcrtBits,
                                         ep, keyPair, numSlotsCKKS, depth,
                                         levelsAvailableBeforeBootstrap, levelsComputation);
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> elapsed = end - start;

    // Print results per group
    for (size_t g = 0; g < results.size(); ++g) {
        DecryptAndPrintRLWE(results[g], cc, Q, scaleTHI, levelsComputation, BigInteger(2), keyPair, ep,
                            numSlotsCKKS, numSlots,
                            "Computed group " + std::to_string(g) + " (0/1)");

        // Expected: for each slot i, 1 if data[i] == groupRef[g], else 0
        std::cerr << "  Expected group " << g << " (0/1): [";
        for (size_t i = 0; i < numSlots; ++i) {
            if (i < data.size()) {
                std::cerr << (data[i] == groupRef[g] ? 1 : 0) << " ";
            } else {
                std::cerr << "0 ";
            }
        }
        std::cerr << "]" << std::endl;
    }

    std::cout << "  Time: " << elapsed.count() << " ms" << std::endl;
}


int main() {
    std::cerr << "\n===== GroupEqual Tests =====" << std::endl;
    std::cerr << "Testing group-based equality comparison" << std::endl;
    std::cerr << "================================\n" << std::endl;

    // ===== Setup =====
    BigInteger PInput(BigInteger(4096));
    BigInteger POutput(BigInteger(4096));
    BigInteger Q(BigInteger(1) << 47);
    BigInteger Bigq(BigInteger(1) << 47);
    uint64_t scaleTHI     = 32;
    size_t order          = 1;
    uint32_t ringDim      = 4096;
    uint32_t numSlots     = 8;

    bool flagSP       = (numSlots <= ringDim / 2);
    auto numSlotsCKKS = flagSP ? numSlots : numSlots / 2;

    uint32_t dcrtBits                       = Bigq.GetMSB() - 1;
    uint32_t firstMod                       = Bigq.GetMSB() - 1;
    uint32_t levelsAvailableAfterBootstrap  = 5;
    uint32_t levelsAvailableBeforeBootstrap = 0;
    uint32_t levelsComputation              = 5;
    uint32_t dnum                           = 3;
    SecretKeyDist secretKeyDist             = SPARSE_ENCAPSULATED;
    std::vector<uint32_t> lvlb              = {1, 1};

    auto a = PInput.ConvertToInt<int64_t>();
    auto func = [a](int64_t x) -> int64_t {
        return (x % a) == 0 ? 1 : 0;
    };

    std::vector<int64_t> coeffint;
    std::vector<std::complex<double>> coeffcomp;
    bool binaryLUT = (PInput.ConvertToInt() == 2) && (order == 1);

    if (binaryLUT) {
        coeffint = {func(1), func(0) - func(1), 0};
    }
    else {
        coeffcomp = GetHermiteTrigCoefficients(func, PInput.ConvertToInt(), order, scaleTHI);
    }

    CCParams<CryptoContextCKKSRNS> parameters;
    parameters.SetSecretKeyDist(secretKeyDist);
    parameters.SetSecurityLevel(HEStd_NotSet);
    parameters.SetScalingModSize(dcrtBits);
    parameters.SetScalingTechnique(FIXEDMANUAL);
    parameters.SetFirstModSize(firstMod);
    parameters.SetNumLargeDigits(dnum);
    parameters.SetBatchSize(numSlotsCKKS);
    parameters.SetRingDim(ringDim);

    uint32_t depth = levelsAvailableAfterBootstrap + levelsComputation;
    if (binaryLUT)
        depth += FHECKKSRNS::GetFBTDepth(lvlb, coeffint, PInput, order, secretKeyDist);
    else
        depth += FHECKKSRNS::GetFBTDepth(lvlb, coeffcomp, PInput, order, secretKeyDist);
    parameters.SetMultiplicativeDepth(depth);

    auto cc = GenCryptoContext(parameters);
    cc->Enable(PKE);
    cc->Enable(KEYSWITCH);
    cc->Enable(LEVELEDSHE);
    cc->Enable(ADVANCEDSHE);
    cc->Enable(FHE);

    auto keyPair = cc->KeyGen();

    std::cout << "CKKS scheme: ringDim=" << cc->GetRingDimension()
              << ", depth=" << depth << std::endl;

    cc->EvalMultKeyGen(keyPair.secretKey);

    if (binaryLUT) {
        cc->EvalFBTSetup(coeffint, numSlotsCKKS, PInput, POutput, Bigq, keyPair.publicKey, {0, 0}, lvlb,
                         levelsAvailableAfterBootstrap, levelsComputation, order);
    }
    else {
        cc->EvalFBTSetup(coeffcomp, numSlotsCKKS, PInput, POutput, Bigq, keyPair.publicKey, {0, 0}, lvlb,
                         levelsAvailableAfterBootstrap, levelsComputation, order);
    }

    cc->EvalBootstrapKeyGen(keyPair.secretKey, numSlotsCKKS);

    // Generate EvalSum keys (needed for broadcasting reference value)
    cc->EvalSumKeyGen(keyPair.secretKey);

    // Generate rotation keys needed for group operations
    std::vector<int32_t> rotIndices;
    for (int32_t r = -static_cast<int32_t>(numSlots); r <= static_cast<int32_t>(numSlots); ++r) {
        if (r != 0) rotIndices.push_back(r);
    }
    cc->EvalAtIndexKeyGen(keyPair.secretKey, rotIndices);

    auto ep = SchemeletRLWEMP::GetElementParams(keyPair.secretKey, depth - (levelsAvailableBeforeBootstrap > 0));

    // ===== GroupEqual Tests =====

    // Test 1: 2 groups, same values within each group
    RunGroupEqualTest(cc, keyPair, ep, PInput, POutput, Q, Bigq,
                      scaleTHI, order, numSlots, ringDim, dcrtBits,
                      numSlotsCKKS, depth, levelsAvailableBeforeBootstrap, levelsComputation,
                      {100, 100, 200, 200, 100, 100, 200, 200},
                      {0, 0, 1, 1, 0, 0, 1, 1},
                      "Test 1: 2 groups, same values in each group");

    // Test 2: 2 groups, mixed values
    RunGroupEqualTest(cc, keyPair, ep, PInput, POutput, Q, Bigq,
                      scaleTHI, order, numSlots, ringDim, dcrtBits,
                      numSlotsCKKS, depth, levelsAvailableBeforeBootstrap, levelsComputation,
                      {50, 100, 50, 100, 50, 100, 50, 100},
                      {0, 1, 0, 1, 0, 1, 0, 1},
                      "Test 2: 2 groups alternating, different values");

    // Test 3: 3 groups (mapped to 1/2/3)
    RunGroupEqualTest(cc, keyPair, ep, PInput, POutput, Q, Bigq,
                      scaleTHI, order, numSlots, ringDim, dcrtBits,
                      numSlotsCKKS, depth, levelsAvailableBeforeBootstrap, levelsComputation,
                      {10, 10, 20, 20, 30, 30, 20, 10},
                      {0, 0, 1, 1, 2, 2, 1, 0},
                      "Test 3: 3 groups (0,1,2)");

    // Test 4: All same group, all same values
    RunGroupEqualTest(cc, keyPair, ep, PInput, POutput, Q, Bigq,
                      scaleTHI, order, numSlots, ringDim, dcrtBits,
                      numSlotsCKKS, depth, levelsAvailableBeforeBootstrap, levelsComputation,
                      {42, 42, 42, 42, 42, 42, 42, 42},
                      {0, 0, 0, 0, 0, 0, 0, 0},
                      "Test 4: Single group, all equal");

    // Test 5: All same group, different values
    RunGroupEqualTest(cc, keyPair, ep, PInput, POutput, Q, Bigq,
                      scaleTHI, order, numSlots, ringDim, dcrtBits,
                      numSlotsCKKS, depth, levelsAvailableBeforeBootstrap, levelsComputation,
                      {10, 20, 30, 40, 50, 60, 70, 80},
                      {0, 0, 0, 0, 0, 0, 0, 0},
                      "Test 5: Single group, all different");

    // ===== GroupEqualMultiGroup Tests =====

    // Test 6: MultiGroup with 2 groups
    RunMultiGroupTest(cc, keyPair, ep, PInput, POutput, Q, Bigq,
                      scaleTHI, order, numSlots, ringDim, dcrtBits,
                      numSlotsCKKS, depth, levelsAvailableBeforeBootstrap, levelsComputation,
                      {100, 100, 200, 200, 100, 200, 100, 200},
                      {0, 0, 1, 1, 0, 1, 0, 1},
                      2,
                      "Test 6: MultiGroup 2 groups");

    // Test 7: MultiGroup with 3 groups
    RunMultiGroupTest(cc, keyPair, ep, PInput, POutput, Q, Bigq,
                      scaleTHI, order, numSlots, ringDim, dcrtBits,
                      numSlotsCKKS, depth, levelsAvailableBeforeBootstrap, levelsComputation,
                      {10, 10, 20, 20, 30, 30, 20, 10},
                      {0, 0, 1, 1, 2, 2, 1, 0},
                      3,
                      "Test 7: MultiGroup 3 groups");

    std::cerr << "\n===== All GroupEqual Tests Completed! =====" << std::endl;

    return 0;
}
