#define PROFILE
#include "comparsion.h"

#include <iostream>
#include <iomanip>
#include <vector>
#include <string>

using namespace lbcrypto;

const BigInteger QBFVINIT(BigInteger(1) << 60);
const BigInteger QBFVINITLARGE(BigInteger(1) << 80);





/**
 * @brief Test all 4 comparison operations on a single dataset
 */
void RunAllComparisonTests(
    lbcrypto::CryptoContextCKKSRNS::ContextType& cc,
    lbcrypto::KeyPair<lbcrypto::DCRTPoly>& keyPair,
    std::shared_ptr<lbcrypto::M4DCRTParams>& ep,
    const BigInteger& PInput,
    const BigInteger& PDigit,
    const BigInteger& Q,
    const BigInteger& Bigq,
    uint64_t scaleTHI,
    uint64_t scaleStepTHI,
    size_t order,
    uint32_t numSlots,
    uint32_t ringDim,
    uint32_t dcrtBits,
    uint32_t numSlotsCKKS,
    uint32_t depth,
    uint32_t levelsAvailableBeforeBootstrap,
    uint32_t levelsComputation,
    const std::vector<int64_t>& x_input,
    const std::vector<int64_t>& y_input,
    const std::string& testName) {

    std::cerr << "\n===== " << testName << " =====" << std::endl;
    std::cerr << "x = " << x_input << std::endl;
    std::cerr << "y = " << y_input << std::endl;

    // Prepare inputs
    auto x = x_input;
    if (x.size() < numSlots) x = Fill<int64_t>(x, numSlots);
    auto y = y_input;
    if (y.size() < numSlots) y = Fill<int64_t>(y, numSlots);

    // Encrypt
    auto ctxtX = SchemeletRLWEMP::EncryptCoeff(x, QBFVINIT, PInput, keyPair.secretKey, ep);
    auto ctxtY = SchemeletRLWEMP::EncryptCoeff(y, QBFVINIT, PInput, keyPair.secretKey, ep);

    // ===== 1. GreaterThan: x > y =====
    auto start = std::chrono::high_resolution_clock::now();
    auto ctxtGT = GreaterThan(ctxtX, ctxtY, cc, QBFVINIT, PInput, PDigit, Q, Bigq,
                               scaleTHI, scaleStepTHI, order, numSlots, ringDim, dcrtBits,
                               ep, keyPair, numSlotsCKKS, depth,
                               levelsAvailableBeforeBootstrap, levelsComputation);
    auto endGT = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> elapsedGT = endGT - start;


    std::cout << "  GT time: " << elapsedGT.count() << " ms" << std::endl;

    // ===== 2. GreaterEqual: x >= y =====
    start = std::chrono::high_resolution_clock::now();
    auto ctxtGE = GreaterEqual(ctxtX, ctxtY, cc, QBFVINIT, PInput, PDigit, Q, Bigq,
                                scaleTHI, scaleStepTHI, order, numSlots, ringDim, dcrtBits,
                                ep, keyPair, numSlotsCKKS, depth,
                                levelsAvailableBeforeBootstrap, levelsComputation);
    auto endGE = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> elapsedGE = endGE - start;


    std::cout << "  GE time: " << elapsedGE.count() << " ms" << std::endl;

    // ===== 3. LessThan: x < y =====
    start = std::chrono::high_resolution_clock::now();
    auto ctxtLT = LessThan(ctxtX, ctxtY, cc, QBFVINIT, PInput, PDigit, Q, Bigq,
                            scaleTHI, scaleStepTHI, order, numSlots, ringDim, dcrtBits,
                            ep, keyPair, numSlotsCKKS, depth,
                            levelsAvailableBeforeBootstrap, levelsComputation);
    auto endLT = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> elapsedLT = endLT - start;


    std::cout << "  LT time: " << elapsedLT.count() << " ms" << std::endl;

    // ===== 4. LessEqual: x <= y =====
    start = std::chrono::high_resolution_clock::now();
    auto ctxtLE = LessEqual(ctxtX, ctxtY, cc, QBFVINIT, PInput, PDigit, Q, Bigq,
                             scaleTHI, scaleStepTHI, order, numSlots, ringDim, dcrtBits,
                             ep, keyPair, numSlotsCKKS, depth,
                             levelsAvailableBeforeBootstrap, levelsComputation);
    auto endLE = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> elapsedLE = endLE - start;



    std::cout << "  LE time: " << elapsedLE.count() << " ms" << std::endl;

    std::cerr << std::endl;
}


int main() {
    std::cerr << "\n===== Homomorphic Comparison Tests =====" << std::endl;
    std::cerr << "Testing: GreaterThan, GreaterEqual, LessThan, LessEqual" << std::endl;
    std::cerr << "========================================\n" << std::endl;

    // ===== Setup parameters (shared across all tests) =====
    BigInteger PInput(BigInteger(4096));
    BigInteger PDigit(BigInteger(2));
    BigInteger Q(BigInteger(1) << 46);
    BigInteger Bigq(BigInteger(1) << 35);
    uint64_t scaleTHI     = 1;
    uint64_t scaleStepTHI = 1;
    size_t order          = 1;
    uint32_t ringDim      = 2048;
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
    auto b = PDigit.ConvertToInt<int64_t>();

    auto funcMod = [b](int64_t x) -> int64_t { return (x % b); };
    auto funcStep = [a, b](int64_t x) -> int64_t { return (x % a) >= (b / 2); };

    std::vector<int64_t> coeffintMod;
    std::vector<std::complex<double>> coeffcompMod;
    std::vector<std::complex<double>> coeffcompStep;

    bool binaryLUT = (PDigit.ConvertToInt() == 2) && (order == 1);

    if (binaryLUT) {
        coeffintMod = {funcMod(1), funcMod(0) - funcMod(1)};
    }
    else {
        coeffcompMod  = GetHermiteTrigCoefficients(funcMod, PDigit.ConvertToInt(), order, scaleTHI);
        coeffcompStep = GetHermiteTrigCoefficients(funcStep, PDigit.ConvertToInt(), order, scaleStepTHI);
    }

    // ===== CryptoContext =====
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
        depth += FHECKKSRNS::GetFBTDepth(lvlb, coeffintMod, PDigit, order, secretKeyDist);
    else
        depth += FHECKKSRNS::GetFBTDepth(lvlb, coeffcompMod, PDigit, order, secretKeyDist);

    parameters.SetMultiplicativeDepth(depth);

    auto cc = GenCryptoContext(parameters);
    cc->Enable(PKE);
    cc->Enable(KEYSWITCH);
    cc->Enable(LEVELEDSHE);
    cc->Enable(ADVANCEDSHE);
    cc->Enable(FHE);

    // ===== KeyGen =====
    auto keyPair = cc->KeyGen();

    std::cout << "CKKS scheme: ringDim=" << cc->GetRingDimension()
              << ", depth=" << depth << std::endl;

    cc->EvalMultKeyGen(keyPair.secretKey);

    if (binaryLUT)
        cc->EvalFBTSetup(coeffintMod, numSlotsCKKS, PDigit, PInput, Bigq, keyPair.publicKey, {0, 0}, lvlb,
                         levelsAvailableAfterBootstrap, levelsComputation, order);
    else
        cc->EvalFBTSetup(coeffcompMod, numSlotsCKKS, PDigit, PInput, Bigq, keyPair.publicKey, {0, 0}, lvlb,
                         levelsAvailableAfterBootstrap, levelsComputation, order);

    cc->EvalBootstrapKeyGen(keyPair.secretKey, numSlotsCKKS);
    cc->EvalAtIndexKeyGen(keyPair.secretKey, std::vector<int32_t>({-2}));

    auto ep = SchemeletRLWEMP::GetElementParams(keyPair.secretKey, depth - (levelsAvailableBeforeBootstrap > 0));

    // ===== Test Data Sets =====
    // Each test compares vector x against vector y

    // Test 1: Original data
    RunAllComparisonTests(cc, keyPair, ep, PInput, PDigit, Q, Bigq,
                          scaleTHI, scaleStepTHI, order, numSlots, ringDim, dcrtBits,
                          numSlotsCKKS, depth, levelsAvailableBeforeBootstrap, levelsComputation,
                          {100, 98, 97, 5, 16, 32, 128, 96},
                          {97},
                          "Test 1: Original data");

    // Test 2: All x > y
    RunAllComparisonTests(cc, keyPair, ep, PInput, PDigit, Q, Bigq,
                          scaleTHI, scaleStepTHI, order, numSlots, ringDim, dcrtBits,
                          numSlotsCKKS, depth, levelsAvailableBeforeBootstrap, levelsComputation,
                          {200, 150, 180, 120, 90, 77, 66, 55},
                          {50},
                          "Test 2: All x > y");

    // Test 3: All x < y
    RunAllComparisonTests(cc, keyPair, ep, PInput, PDigit, Q, Bigq,
                          scaleTHI, scaleStepTHI, order, numSlots, ringDim, dcrtBits,
                          numSlotsCKKS, depth, levelsAvailableBeforeBootstrap, levelsComputation,
                          {1, 2, 3, 4, 5, 6, 7, 8},
                          {100},
                          "Test 3: All x < y");

    // Test 4: All x == y
    RunAllComparisonTests(cc, keyPair, ep, PInput, PDigit, Q, Bigq,
                          scaleTHI, scaleStepTHI, order, numSlots, ringDim, dcrtBits,
                          numSlotsCKKS, depth, levelsAvailableBeforeBootstrap, levelsComputation,
                          {50, 50, 50, 50, 50, 50, 50, 50},
                          {50},
                          "Test 4: All x == y");

    // Test 5: Mixed (some greater, some less, some equal)
    RunAllComparisonTests(cc, keyPair, ep, PInput, PDigit, Q, Bigq,
                          scaleTHI, scaleStepTHI, order, numSlots, ringDim, dcrtBits,
                          numSlotsCKKS, depth, levelsAvailableBeforeBootstrap, levelsComputation,
                          {30, 80, 50, 60, 45, 50, 10, 70},
                          {50},
                          "Test 5: Mixed (vs 50)");

    // Test 6: Edge cases - zero and small numbers
    RunAllComparisonTests(cc, keyPair, ep, PInput, PDigit, Q, Bigq,
                          scaleTHI, scaleStepTHI, order, numSlots, ringDim, dcrtBits,
                          numSlotsCKKS, depth, levelsAvailableBeforeBootstrap, levelsComputation,
                          {0, 1, 2, 3, 0, 1, 2, 3},
                          {1},
                          "Test 6: Edge cases (0,1,2,3 vs 1)");

    // Test 7: x == y (all equal to 0)
    RunAllComparisonTests(cc, keyPair, ep, PInput, PDigit, Q, Bigq,
                          scaleTHI, scaleStepTHI, order, numSlots, ringDim, dcrtBits,
                          numSlotsCKKS, depth, levelsAvailableBeforeBootstrap, levelsComputation,
                          {0, 0, 0, 0, 0, 0, 0, 0},
                          {0},
                          "Test 7: All zero");

    // Test 8: Larger range values
    RunAllComparisonTests(cc, keyPair, ep, PInput, PDigit, Q, Bigq,
                          scaleTHI, scaleStepTHI, order, numSlots, ringDim, dcrtBits,
                          numSlotsCKKS, depth, levelsAvailableBeforeBootstrap, levelsComputation,
                          {500, 1000, 250, 750, 125, 875, 50, 950},
                          {500},
                          "Test 8: Larger range (vs 500)");

    // Test 9: x == y with varied values
    RunAllComparisonTests(cc, keyPair, ep, PInput, PDigit, Q, Bigq,
                          scaleTHI, scaleStepTHI, order, numSlots, ringDim, dcrtBits,
                          numSlotsCKKS, depth, levelsAvailableBeforeBootstrap, levelsComputation,
                          {42, 99, 7, 13, 56, 88, 21, 77},
                          {42, 99, 7, 13, 56, 88, 21, 77},
                          "Test 9: Exactly equal vectors");

    // Test 10: x vs y with alternating comparison
    RunAllComparisonTests(cc, keyPair, ep, PInput, PDigit, Q, Bigq,
                          scaleTHI, scaleStepTHI, order, numSlots, ringDim, dcrtBits,
                          numSlotsCKKS, depth, levelsAvailableBeforeBootstrap, levelsComputation,
                          {10, 20, 30, 40, 50, 60, 70, 80},
                          {80, 70, 60, 50, 40, 30, 20, 10},
                          "Test 10: Alternating (10..80 vs 80..10)");

    std::cerr << "\n===== All Comparison Tests Completed! =====" << std::endl;

    return 0;
}
