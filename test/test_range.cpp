#define PROFILE
#include "range.h"

#include <iostream>
#include <iomanip>
#include <vector>
#include <string>

using namespace lbcrypto;

const BigInteger QBFVINIT(BigInteger(1) << 60);
const BigInteger QBFVINITLARGE(BigInteger(1) << 80);


void RunRangeTest(
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
    const std::vector<int64_t>& x_input,
    const std::vector<int64_t>& a_input,
    const std::vector<int64_t>& b_input,
    const std::string& testName) {

    std::cerr << "\n===== " << testName << " =====" << std::endl;
    std::cerr << "x = " << x_input << std::endl;
    std::cerr << "a = " << a_input << ", b = " << b_input << std::endl;

    auto x = x_input;
    if (x.size() < numSlots) x = Fill<int64_t>(x, numSlots);
    auto a = a_input;
    if (a.size() < numSlots) a = Fill<int64_t>(a, numSlots);
    auto b = b_input;
    if (b.size() < numSlots) b = Fill<int64_t>(b, numSlots);

    auto ctxtX = SchemeletRLWEMP::EncryptCoeff(x, QBFVINIT, PInput, keyPair.secretKey, ep);
    auto ctxtA = SchemeletRLWEMP::EncryptCoeff(a, QBFVINIT, PInput, keyPair.secretKey, ep);
    auto ctxtB = SchemeletRLWEMP::EncryptCoeff(b, QBFVINIT, PInput, keyPair.secretKey, ep);

    // 1. HomRangeStrict: a < x < b
    HomRangeStrict(ctxtX, ctxtA, ctxtB, cc, QBFVINIT, PInput, PDigit, Q, Bigq,
                   scaleTHI, scaleStepTHI, order, numSlots, ringDim, dcrtBits,
                   ep, keyPair, numSlotsCKKS, depth,
                   levelsAvailableBeforeBootstrap, levelsComputation);

    // 2. HomRangeLowerInc: a <= x < b
    HomRangeLowerInc(ctxtX, ctxtA, ctxtB, cc, QBFVINIT, PInput, PDigit, Q, Bigq,
                     scaleTHI, scaleStepTHI, order, numSlots, ringDim, dcrtBits,
                     ep, keyPair, numSlotsCKKS, depth,
                     levelsAvailableBeforeBootstrap, levelsComputation);

    // 3. HomRangeUpperInc: a < x <= b
    HomRangeUpperInc(ctxtX, ctxtA, ctxtB, cc, QBFVINIT, PInput, PDigit, Q, Bigq,
                     scaleTHI, scaleStepTHI, order, numSlots, ringDim, dcrtBits,
                     ep, keyPair, numSlotsCKKS, depth,
                     levelsAvailableBeforeBootstrap, levelsComputation);

    // 4. HomRangeBothInc: a <= x <= b
    HomRangeBothInc(ctxtX, ctxtA, ctxtB, cc, QBFVINIT, PInput, PDigit, Q, Bigq,
                    scaleTHI, scaleStepTHI, order, numSlots, ringDim, dcrtBits,
                    ep, keyPair, numSlotsCKKS, depth,
                    levelsAvailableBeforeBootstrap, levelsComputation);

    std::cerr << std::endl;
}


int main() {
    std::cerr << "\n===== Homomorphic Range Tests =====" << std::endl;
    std::cerr << "Testing: HomRangeStrict, HomRangeLowerInc, HomRangeUpperInc, HomRangeBothInc" << std::endl;
    std::cerr << "========================================\n" << std::endl;

    // ===== Setup parameters =====
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
    // Test 1: Original data (a=18, b=97)
    RunRangeTest(cc, keyPair, ep, PInput, PDigit, Q, Bigq,
                 scaleTHI, scaleStepTHI, order, numSlots, ringDim, dcrtBits,
                 numSlotsCKKS, depth, levelsAvailableBeforeBootstrap, levelsComputation,
                 {100, 98, 97, 5, 16, 32, 128, 96},
                 {18},
                 {97},
                 "Test 1: Original data (a=18, b=97)");

    // Test 2: All x in range
    RunRangeTest(cc, keyPair, ep, PInput, PDigit, Q, Bigq,
                 scaleTHI, scaleStepTHI, order, numSlots, ringDim, dcrtBits,
                 numSlotsCKKS, depth, levelsAvailableBeforeBootstrap, levelsComputation,
                 {30, 40, 50, 60, 70, 80, 90, 95},
                 {20},
                 {100},
                 "Test 2: All x in range (20 < x < 100)");

    // Test 3: All x out of range (below)
    RunRangeTest(cc, keyPair, ep, PInput, PDigit, Q, Bigq,
                 scaleTHI, scaleStepTHI, order, numSlots, ringDim, dcrtBits,
                 numSlotsCKKS, depth, levelsAvailableBeforeBootstrap, levelsComputation,
                 {1, 2, 3, 4, 5, 6, 7, 8},
                 {10},
                 {100},
                 "Test 3: All x below range (10 < x < 100)");

    // Test 4: All x out of range (above)
    RunRangeTest(cc, keyPair, ep, PInput, PDigit, Q, Bigq,
                 scaleTHI, scaleStepTHI, order, numSlots, ringDim, dcrtBits,
                 numSlotsCKKS, depth, levelsAvailableBeforeBootstrap, levelsComputation,
                 {101, 150, 200, 250, 300, 350, 400, 500},
                 {10},
                 {100},
                 "Test 4: All x above range (10 < x < 100)");

    // Test 5: x exactly at boundaries
    RunRangeTest(cc, keyPair, ep, PInput, PDigit, Q, Bigq,
                 scaleTHI, scaleStepTHI, order, numSlots, ringDim, dcrtBits,
                 numSlotsCKKS, depth, levelsAvailableBeforeBootstrap, levelsComputation,
                 {10, 100, 50, 50, 50, 50, 50, 50},
                 {10},
                 {100},
                 "Test 5: Boundary values (x=10, x=100)");

    // Test 6: Mixed in/out
    RunRangeTest(cc, keyPair, ep, PInput, PDigit, Q, Bigq,
                 scaleTHI, scaleStepTHI, order, numSlots, ringDim, dcrtBits,
                 numSlotsCKKS, depth, levelsAvailableBeforeBootstrap, levelsComputation,
                 {5, 25, 50, 75, 95, 105, 30, 80},
                 {20},
                 {100},
                 "Test 6: Mixed values (20 < x < 100)");

    std::cerr << "\n===== All Range Tests Completed! =====" << std::endl;

    return 0;
}
