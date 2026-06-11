#define PROFILE
#include "multi_precision.h"

#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <chrono>
#include <cmath>

using namespace lbcrypto;

const BigInteger QBFVINIT(BigInteger(1) << 60);

/**
 * @brief RunMultiPrecisionTest - runs a multi-precision comparison test
 *
 * @param x_input  Plaintext value of X (full integer, will be split into digits)
 * @param y_input  Plaintext value of Y (full integer, will be split into digits)
 * @param digitBits Number of bits per digit (e.g., 8)
 * @param numDigits Number of digits
 * @param testName  Name of this test
 */
void RunMultiPrecisionTest(
    CryptoContextCKKSRNS::ContextType& cc,
    KeyPair<DCRTPoly>& keyPair,
    std::shared_ptr<M4DCRTParams>& ep,
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
    const std::vector<int64_t>& x_values,
    const std::vector<int64_t>& y_values,
    uint32_t digitBits,
    uint32_t numDigits,
    const std::string& testName) {

    std::cerr << "\n===== " << testName << " =====" << std::endl;

    int64_t digitMask = (1LL << digitBits) - 1;
    size_t N = x_values.size();

    // Split each value into digitBits-bit digits, MSB first
    // Each digit[i][slot] = (value[slot] >> ((numDigits-1-i) * digitBits)) & digitMask
    std::vector<std::vector<int64_t>> x_digits(numDigits, std::vector<int64_t>(N, 0));
    std::vector<std::vector<int64_t>> y_digits(numDigits, std::vector<int64_t>(N, 0));

    for (size_t s = 0; s < N; s++) {
        for (uint32_t d = 0; d < numDigits; d++) {
            int shift = (numDigits - 1 - d) * digitBits;
            x_digits[d][s] = (x_values[s] >> shift) & digitMask;
            y_digits[d][s] = (y_values[s] >> shift) & digitMask;
        }
    }

    // Print expected comparisons
    std::cerr << "X values: [";
    for (size_t s = 0; s < N; s++) std::cerr << x_values[s] << (s < N-1 ? ", " : "");
    std::cerr << "]" << std::endl;

    std::cerr << "Y values: [";
    for (size_t s = 0; s < N; s++) std::cerr << y_values[s] << (s < N-1 ? ", " : "");
    std::cerr << "]" << std::endl;

    std::cerr << "Expected X > Y: [";
    for (size_t s = 0; s < N; s++) std::cerr << (x_values[s] > y_values[s] ? 1 : 0) << (s < N-1 ? ", " : "");
    std::cerr << "]" << std::endl;

    std::cerr << "Expected X == Y: [";
    for (size_t s = 0; s < N; s++) std::cerr << (x_values[s] == y_values[s] ? 1 : 0) << (s < N-1 ? ", " : "");
    std::cerr << "]" << std::endl;

    // Print digits for debugging
    for (uint32_t d = 0; d < numDigits; d++) {
        std::cerr << "  Digit[" << d << "] (shift=" << ((numDigits-1-d)*digitBits) << "): X=[";
        for (size_t s = 0; s < N; s++) std::cerr << x_digits[d][s] << (s < N-1 ? ", " : "");
        std::cerr << "], Y=[";
        for (size_t s = 0; s < N; s++) std::cerr << y_digits[d][s] << (s < N-1 ? ", " : "");
        std::cerr << "]" << std::endl;
    }

    // Fill to numSlots
    for (uint32_t d = 0; d < numDigits; d++) {
        if (x_digits[d].size() < numSlots)
            x_digits[d] = Fill<int64_t>(x_digits[d], numSlots);
        if (y_digits[d].size() < numSlots)
            y_digits[d] = Fill<int64_t>(y_digits[d], numSlots);
    }

    // Encrypt each digit
    std::vector<std::vector<Poly>> ctxtX(numDigits);
    std::vector<std::vector<Poly>> ctxtY(numDigits);

    for (uint32_t d = 0; d < numDigits; d++) {
        ctxtX[d] = SchemeletRLWEMP::EncryptCoeff(x_digits[d], QBFVINIT, PInput, keyPair.secretKey, ep);
        ctxtY[d] = SchemeletRLWEMP::EncryptCoeff(y_digits[d], QBFVINIT, PInput, keyPair.secretKey, ep);
    }

    // Start case timer
    auto case_start = std::chrono::high_resolution_clock::now();

    // 1. MultiPrecisionGreaterThan: X > Y
    auto t0 = std::chrono::high_resolution_clock::now();
    MultiPrecisionGreaterThan(ctxtX, ctxtY, numDigits, cc, QBFVINIT, PInput, PDigit, Q, Bigq,
                              scaleTHI, scaleStepTHI, order, numSlots, ringDim, dcrtBits,
                              ep, keyPair, numSlotsCKKS, depth,
                              levelsAvailableBeforeBootstrap, levelsComputation);
    auto t1 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> dt_gt = t1 - t0;
    std::cerr << "  Timing: GreaterThan = " << dt_gt.count() << " ms" << std::endl;

    // 2. MultiPrecisionEqual: X == Y
    t0 = std::chrono::high_resolution_clock::now();
    MultiPrecisionEqual(ctxtX, ctxtY, numDigits, cc, QBFVINIT, PInput, PDigit, Q, Bigq,
                        scaleTHI, scaleStepTHI, order, numSlots, ringDim, dcrtBits,
                        ep, keyPair, numSlotsCKKS, depth,
                        levelsAvailableBeforeBootstrap, levelsComputation);
    t1 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> dt_eq = t1 - t0;
    std::cerr << "  Timing: Equal = " << dt_eq.count() << " ms" << std::endl;

    // 3. MultiPrecisionLessThan: X < Y
    t0 = std::chrono::high_resolution_clock::now();
    MultiPrecisionLessThan(ctxtX, ctxtY, numDigits, cc, QBFVINIT, PInput, PDigit, Q, Bigq,
                           scaleTHI, scaleStepTHI, order, numSlots, ringDim, dcrtBits,
                           ep, keyPair, numSlotsCKKS, depth,
                           levelsAvailableBeforeBootstrap, levelsComputation);
    t1 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> dt_lt = t1 - t0;
    std::cerr << "  Timing: LessThan = " << dt_lt.count() << " ms" << std::endl;

    // 4. MultiPrecisionGreaterEqual: X >= Y
    t0 = std::chrono::high_resolution_clock::now();
    MultiPrecisionGreaterEqual(ctxtX, ctxtY, numDigits, cc, QBFVINIT, PInput, PDigit, Q, Bigq,
                               scaleTHI, scaleStepTHI, order, numSlots, ringDim, dcrtBits,
                               ep, keyPair, numSlotsCKKS, depth,
                               levelsAvailableBeforeBootstrap, levelsComputation);
    t1 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> dt_ge = t1 - t0;
    std::cerr << "  Timing: GreaterEqual = " << dt_ge.count() << " ms" << std::endl;

    // 5. MultiPrecisionLessEqual: X <= Y
    t0 = std::chrono::high_resolution_clock::now();
    MultiPrecisionLessEqual(ctxtX, ctxtY, numDigits, cc, QBFVINIT, PInput, PDigit, Q, Bigq,
                            scaleTHI, scaleStepTHI, order, numSlots, ringDim, dcrtBits,
                            ep, keyPair, numSlotsCKKS, depth,
                            levelsAvailableBeforeBootstrap, levelsComputation);
    t1 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> dt_le = t1 - t0;
    std::cerr << "  Timing: LessEqual = " << dt_le.count() << " ms" << std::endl;

    auto case_end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> dt_case = case_end - case_start;
    std::cerr << "  Timing: Total case = " << dt_case.count() << " ms" << std::endl;

    std::cerr << std::endl;
}


int main() {
    std::cerr << "\n===== Multi-Precision Homomorphic Comparison Tests =====" << std::endl;
    std::cerr << "Testing: multi-digit sign comparison using recursive digit-wise logic" << std::endl;
    std::cerr << "Algorithm: result[i] = greater_i + equal_i * result[i+1]" << std::endl;
    std::cerr << "============================================================\n" << std::endl;

    // ===== Setup parameters =====
    uint32_t digitBits = 8;              // bits per digit
    uint32_t numDigits = 4;              // total digits -> 32-bit numbers
    BigInteger PInput(BigInteger(1) << digitBits);  // 256
    BigInteger PDigit(BigInteger(2));    // binary digit for sign
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

    // For sign (binary LUT)
    auto a = PInput.ConvertToInt<int64_t>();
    auto b = PDigit.ConvertToInt<int64_t>();

    auto funcMod = [b](int64_t x) -> int64_t { return (x % b); };
    auto funcStep = [a, b](int64_t x) -> int64_t { return (x % a) >= (b / 2); };

    // Equal LUT
    auto funcEqual = [a](int64_t x) -> int64_t { return (x % a) == 0 ? 1 : 0; };

    std::vector<int64_t> coeffintMod;
    std::vector<std::complex<double>> coeffcompMod;
    std::vector<std::complex<double>> coeffcompStep;

    // Equal coefficients
    std::vector<int64_t> coeffintEqual;
    std::vector<std::complex<double>> coeffcompEqual;

    bool binaryLUT = (PDigit.ConvertToInt() == 2) && (order == 1);
    bool binaryLUTEqual = (PInput.ConvertToInt() == 2) && (order == 1);  // equal uses PInput as digit

    if (binaryLUT) {
        coeffintMod = {funcMod(1), funcMod(0) - funcMod(1)};
    } else {
        coeffcompMod  = GetHermiteTrigCoefficients(funcMod, PDigit.ConvertToInt(), order, scaleTHI);
        coeffcompStep = GetHermiteTrigCoefficients(funcStep, PDigit.ConvertToInt(), order, scaleStepTHI);
    }

    // Equal LUT: uses PInput=256 as digit, never binary LUT (PInput != 2)
    if (binaryLUTEqual) {
        coeffintEqual = {funcEqual(1), funcEqual(0) - funcEqual(1), 0};
    } else {
        coeffcompEqual = GetHermiteTrigCoefficients(funcEqual, PInput.ConvertToInt(), order, scaleTHI);
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
    std::cerr << "DEBUG: initial depth = " << depth << std::endl;

    // Sign FBT depth
    if (binaryLUT) {
        uint32_t signFBTDepth = FHECKKSRNS::GetFBTDepth(lvlb, coeffintMod, PDigit, order, secretKeyDist);
        std::cerr << "DEBUG: sign FBT depth = " << signFBTDepth << std::endl;
        depth += signFBTDepth;
    } else {
        uint32_t signFBTDepth = FHECKKSRNS::GetFBTDepth(lvlb, coeffcompMod, PDigit, order, secretKeyDist);
        std::cerr << "DEBUG: sign FBT depth = " << signFBTDepth << std::endl;
        depth += signFBTDepth;
    }

    // Equal FBT depth
    if (binaryLUTEqual) {
        uint32_t equalFBTDepth = FHECKKSRNS::GetFBTDepth(lvlb, coeffintEqual, PInput, order, secretKeyDist);
        std::cerr << "DEBUG: equal FBT depth = " << equalFBTDepth << std::endl;
        depth += equalFBTDepth;
    } else {
        uint32_t equalFBTDepth = FHECKKSRNS::GetFBTDepth(lvlb, coeffcompEqual, PInput, order, secretKeyDist);
        std::cerr << "DEBUG: equal FBT depth = " << equalFBTDepth << std::endl;
        depth += equalFBTDepth;
    }

    // Extra depth for EvalMult in recursion (1 mult per digit)
    depth += numDigits;

    std::cerr << "DEBUG: final depth = " << depth << std::endl;

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
              << ", depth=" << depth
              << ", digitBits=" << digitBits
              << ", numDigits=" << numDigits << std::endl;

    cc->EvalMultKeyGen(keyPair.secretKey);

    // Setup for sign FBT
    if (binaryLUT)
        cc->EvalFBTSetup(coeffintMod, numSlotsCKKS, PDigit, PInput, Bigq, keyPair.publicKey, {0, 0}, lvlb,
                         levelsAvailableAfterBootstrap, levelsComputation, order);
    else
        cc->EvalFBTSetup(coeffcompMod, numSlotsCKKS, PDigit, PInput, Bigq, keyPair.publicKey, {0, 0}, lvlb,
                         levelsAvailableAfterBootstrap, levelsComputation, order);

    // Setup for equal FBT (using PInput as digit for the full-width equal check)
    if (binaryLUTEqual)
        cc->EvalFBTSetup(coeffintEqual, numSlotsCKKS, PInput, PInput, Bigq, keyPair.publicKey, {0, 0}, lvlb,
                         levelsAvailableAfterBootstrap, levelsComputation, order);
    else
        cc->EvalFBTSetup(coeffcompEqual, numSlotsCKKS, PInput, PInput, Bigq, keyPair.publicKey, {0, 0}, lvlb,
                         levelsAvailableAfterBootstrap, levelsComputation, order);

    cc->EvalBootstrapKeyGen(keyPair.secretKey, numSlotsCKKS);
    cc->EvalAtIndexKeyGen(keyPair.secretKey, std::vector<int32_t>({-2}));

    auto ep = SchemeletRLWEMP::GetElementParams(keyPair.secretKey, depth - (levelsAvailableBeforeBootstrap > 0));

    // ===== Test Data Sets =====
    // 32-bit values (4 digits of 8 bits each)

    // Test 1: Simple comparison, X > Y
    RunMultiPrecisionTest(cc, keyPair, ep, PInput, PDigit, Q, Bigq,
                          scaleTHI, scaleStepTHI, order, numSlots, ringDim, dcrtBits,
                          numSlotsCKKS, depth, levelsAvailableBeforeBootstrap, levelsComputation,
                          {1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000},
                          {500,  1500, 2500, 3500, 4500, 5500, 6500, 7500},
                          digitBits, numDigits,
                          "Test 1: X > Y for all slots");

    // Test 2: X == Y for some slots
    RunMultiPrecisionTest(cc, keyPair, ep, PInput, PDigit, Q, Bigq,
                          scaleTHI, scaleStepTHI, order, numSlots, ringDim, dcrtBits,
                          numSlotsCKKS, depth, levelsAvailableBeforeBootstrap, levelsComputation,
                          {100, 200, 300, 400, 500, 600, 700, 800},
                          {100, 150, 300, 350, 500, 550, 700, 750},
                          digitBits, numDigits,
                          "Test 2: X == Y for some slots (slots 0,2,4,6 equal)");

    // Test 3: X < Y
    RunMultiPrecisionTest(cc, keyPair, ep, PInput, PDigit, Q, Bigq,
                          scaleTHI, scaleStepTHI, order, numSlots, ringDim, dcrtBits,
                          numSlotsCKKS, depth, levelsAvailableBeforeBootstrap, levelsComputation,
                          {100, 200, 300, 400, 500, 600, 700, 800},
                          {900, 900, 900, 900, 900, 900, 900, 900},
                          digitBits, numDigits,
                          "Test 3: All X < Y");

    // Test 4: Large values near 32-bit max
    RunMultiPrecisionTest(cc, keyPair, ep, PInput, PDigit, Q, Bigq,
                          scaleTHI, scaleStepTHI, order, numSlots, ringDim, dcrtBits,
                          numSlotsCKKS, depth, levelsAvailableBeforeBootstrap, levelsComputation,
                          {1000000, 2000000, 500000, 100000, 9999999, 1, 65536, 256},
                          {500000,  1500000, 800000, 50000,  8888888, 2, 32768, 512},
                          digitBits, numDigits,
                          "Test 4: Large values (up to 10M)");

    // Test 5: Equal values (all slots)
    RunMultiPrecisionTest(cc, keyPair, ep, PInput, PDigit, Q, Bigq,
                          scaleTHI, scaleStepTHI, order, numSlots, ringDim, dcrtBits,
                          numSlotsCKKS, depth, levelsAvailableBeforeBootstrap, levelsComputation,
                          {12345, 54321, 11111, 22222, 33333, 44444, 55555, 66666},
                          {12345, 54321, 11111, 22222, 33333, 44444, 55555, 66666},
                          digitBits, numDigits,
                          "Test 5: All X == Y");

    std::cerr << "\n===== All Multi-Precision Tests Completed! =====" << std::endl;

    return 0;
}
