#define PROFILE
#include "aggregation.h"

#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <numeric>

using namespace lbcrypto;

// ===== Helper: decrypt CKKS and print =====
void DecryptCKKS(
    const Ciphertext<DCRTPoly>& ctxt,
    CryptoContextCKKSRNS::ContextType& cc,
    KeyPair<DCRTPoly>& keyPair,
    uint32_t numSlots,
    const std::string& label) {

    Plaintext ptxt;
    cc->Decrypt(keyPair.secretKey, ctxt, &ptxt);
    ptxt->SetLength(numSlots);
    auto values = ptxt->GetRealPackedValue();

    std::cerr << "  " << label << ": [";
    for (size_t i = 0; i < values.size() && i < numSlots; ++i) {
        std::cerr << std::fixed << std::setprecision(1) << values[i] << " ";
    }
    std::cerr << "]" << std::endl;
}

// ===== Helper: decrypt RLWE and print =====
void DecryptRLWE(
    const Ciphertext<DCRTPoly>& ctxt,
    CryptoContextCKKSRNS::ContextType& cc,
    const BigInteger& Q,
    uint64_t scaleTHI,
    uint32_t levelsComputation,
    const BigInteger& POutput,
    KeyPair<DCRTPoly>& keyPair,
    std::shared_ptr<M4DCRTParams>& ep,
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

// ===== SUM Test =====
void TestSum(
    CryptoContextCKKSRNS::ContextType& cc,
    KeyPair<DCRTPoly>& keyPair,
    const BigInteger& QBFVInit,
    const BigInteger& PInput,
    const BigInteger& Q,
    const BigInteger& Bigq,
    uint32_t numSlots,
    uint32_t numSlotsCKKS,
    uint32_t depth,
    uint32_t levelsAvailableBeforeBootstrap,
    const std::vector<int64_t>& data,
    const std::string& testName) {

    std::cerr << "\n===== " << testName << " =====" << std::endl;
    std::cerr << "  Data: " << data << std::endl;

    auto x = data;
    if (x.size() < numSlots) x = Fill<int64_t>(x, numSlots);

    auto ep = SchemeletRLWEMP::GetElementParams(keyPair.secretKey, depth - (levelsAvailableBeforeBootstrap > 0));
    auto ctxtData = SchemeletRLWEMP::EncryptCoeff(x, QBFVINIT, PInput, keyPair.secretKey, ep);

    auto result = HomSum(ctxtData, cc, QBFVINIT, Q, Bigq, keyPair,
                         numSlotsCKKS, numSlots, depth, levelsAvailableBeforeBootstrap);

    DecryptCKKS(result, cc, keyPair, numSlots, "Computed SUM");

    int64_t expectedSum = 0;
    for (size_t i = 0; i < data.size(); ++i) expectedSum += data[i];
    std::cerr << "  Expected SUM: " << expectedSum << std::endl;
}

// ===== COUNT Test =====
void TestCount(
    CryptoContextCKKSRNS::ContextType& cc,
    KeyPair<DCRTPoly>& keyPair,
    std::shared_ptr<M4DCRTParams>& ep,
    const BigInteger& QBFVInit,
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
    const std::string& testName) {

    std::cerr << "\n===== " << testName << " =====" << std::endl;
    std::cerr << "  Data: " << data << std::endl;

    auto x = data;
    if (x.size() < numSlots) x = Fill<int64_t>(x, numSlots);

    auto ctxtData = SchemeletRLWEMP::EncryptCoeff(x, QBFVINIT, PInput, keyPair.secretKey, ep);

    auto result = HomCount(ctxtData, cc, QBFVInit, PInput, POutput, Q, Bigq,
                           scaleTHI, order, numSlots, ringDim, dcrtBits,
                           ep, keyPair, numSlotsCKKS, depth,
                           levelsAvailableBeforeBootstrap, levelsComputation);

    DecryptCKKS(result, cc, keyPair, numSlots, "Computed COUNT");

    int64_t expectedCount = 0;
    for (size_t i = 0; i < data.size(); ++i) {
        if (data[i] != 0) expectedCount++;
    }
    std::cerr << "  Expected COUNT: " << expectedCount << std::endl;
}

// ===== AVG Test =====
void TestAvg(
    CryptoContextCKKSRNS::ContextType& cc,
    KeyPair<DCRTPoly>& keyPair,
    std::shared_ptr<M4DCRTParams>& ep,
    const BigInteger& QBFVInit,
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
    const std::string& testName) {

    std::cerr << "\n===== " << testName << " =====" << std::endl;
    std::cerr << "  Data: " << data << std::endl;

    auto x = data;
    if (x.size() < numSlots) x = Fill<int64_t>(x, numSlots);

    auto ctxtData = SchemeletRLWEMP::EncryptCoeff(x, QBFVINIT, PInput, keyPair.secretKey, ep);

    auto [sumResult, countResult] = HomAvg(ctxtData, cc, QBFVInit, PInput, POutput, Q, Bigq,
                                           scaleTHI, order, numSlots, ringDim, dcrtBits,
                                           ep, keyPair, numSlotsCKKS, depth,
                                           levelsAvailableBeforeBootstrap, levelsComputation);

    DecryptCKKS(sumResult, cc, keyPair, numSlots, "Computed SUM");
    DecryptCKKS(countResult, cc, keyPair, numSlots, "Computed COUNT");

    int64_t expectedSum = 0;
    int64_t expectedCount = 0;
    for (size_t i = 0; i < data.size(); ++i) {
        expectedSum += data[i];
        if (data[i] != 0) expectedCount++;
    }
    double expectedAvg = (expectedCount > 0) ? (double)expectedSum / expectedCount : 0;
    std::cerr << "  Expected SUM=" << expectedSum << ", COUNT=" << expectedCount
              << ", AVG=" << std::fixed << std::setprecision(2) << expectedAvg << std::endl;
}

// ===== MAX Test =====
void TestMax(
    CryptoContextCKKSRNS::ContextType& cc,
    KeyPair<DCRTPoly>& keyPair,
    std::shared_ptr<M4DCRTParams>& ep,
    const BigInteger& QBFVInit,
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
    const std::vector<int64_t>& data,
    const std::string& testName) {

    std::cerr << "\n===== " << testName << " =====" << std::endl;
    std::cerr << "  Data: " << data << std::endl;

    auto x = data;
    if (x.size() < numSlots) x = Fill<int64_t>(x, numSlots);

    auto ctxtData = SchemeletRLWEMP::EncryptCoeff(x, QBFVINIT, PInput, keyPair.secretKey, ep);

    auto result = HomMax(ctxtData, cc, QBFVInit, PInput, PDigit, Q, Bigq,
                         scaleTHI, scaleStepTHI, order, numSlots, ringDim, dcrtBits,
                         ep, keyPair, numSlotsCKKS, depth,
                         levelsAvailableBeforeBootstrap, levelsComputation);

    DecryptCKKS(result, cc, keyPair, numSlots, "Computed MAX");

    int64_t expectedMax = data.empty() ? 0 : *std::max_element(data.begin(), data.end());
    std::cerr << "  Expected MAX: " << expectedMax << std::endl;
}

// ===== MIN Test =====
void TestMin(
    CryptoContextCKKSRNS::ContextType& cc,
    KeyPair<DCRTPoly>& keyPair,
    std::shared_ptr<M4DCRTParams>& ep,
    const BigInteger& QBFVInit,
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
    const std::vector<int64_t>& data,
    const std::string& testName) {

    std::cerr << "\n===== " << testName << " =====" << std::endl;
    std::cerr << "  Data: " << data << std::endl;

    auto x = data;
    if (x.size() < numSlots) x = Fill<int64_t>(x, numSlots);

    auto ctxtData = SchemeletRLWEMP::EncryptCoeff(x, QBFVINIT, PInput, keyPair.secretKey, ep);

    auto result = HomMin(ctxtData, cc, QBFVInit, PInput, PDigit, Q, Bigq,
                         scaleTHI, scaleStepTHI, order, numSlots, ringDim, dcrtBits,
                         ep, keyPair, numSlotsCKKS, depth,
                         levelsAvailableBeforeBootstrap, levelsComputation);

    DecryptCKKS(result, cc, keyPair, numSlots, "Computed MIN");

    int64_t expectedMin = data.empty() ? 0 : *std::min_element(data.begin(), data.end());
    std::cerr << "  Expected MIN: " << expectedMin << std::endl;
}


int main() {
    std::cerr << "\n===== Aggregation Tests =====" << std::endl;
    std::cerr << "Testing SUM, COUNT, AVG, MAX, MIN" << std::endl;
    std::cerr << "==============================\n" << std::endl;

    // ===== Common Setup =====
    BigInteger PInput(BigInteger(4096));
    BigInteger POutput(BigInteger(4096));
    BigInteger PDigit(BigInteger(2));
    BigInteger Q(BigInteger(1) << 47);
    BigInteger Bigq(BigInteger(1) << 47);
    uint64_t scaleTHI     = 32;
    uint64_t scaleStepTHI = 32;
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
    } else {
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

    // FBT setup for EqualBootstrapping
    if (binaryLUT) {
        cc->EvalFBTSetup(coeffint, numSlotsCKKS, PInput, POutput, Bigq, keyPair.publicKey, {0, 0}, lvlb,
                         levelsAvailableAfterBootstrap, levelsComputation, order);
    } else {
        cc->EvalFBTSetup(coeffcomp, numSlotsCKKS, PInput, POutput, Bigq, keyPair.publicKey, {0, 0}, lvlb,
                         levelsAvailableAfterBootstrap, levelsComputation, order);
    }

    cc->EvalBootstrapKeyGen(keyPair.secretKey, numSlotsCKKS);

    // EvalSum keys (for SUM/COUNT/AVG)
    cc->EvalSumKeyGen(keyPair.secretKey);

    // Rotation keys (for MAX/MIN tournament)
    std::vector<int32_t> rotIndices;
    for (int32_t r = -static_cast<int32_t>(numSlots); r <= static_cast<int32_t>(numSlots); ++r) {
        if (r != 0) rotIndices.push_back(r);
    }
    cc->EvalAtIndexKeyGen(keyPair.secretKey, rotIndices);

    auto ep = SchemeletRLWEMP::GetElementParams(keyPair.secretKey, depth - (levelsAvailableBeforeBootstrap > 0));

    // ===== SUM Tests =====
    TestSum(cc, keyPair, QBFVINIT, PInput, Q, Bigq,
            numSlots, numSlotsCKKS, depth, levelsAvailableBeforeBootstrap,
            {10, 20, 30, 40, 50, 60, 70, 80},
            "SUM Test 1: simple increasing");

    TestSum(cc, keyPair, QBFVINIT, PInput, Q, Bigq,
            numSlots, numSlotsCKKS, depth, levelsAvailableBeforeBootstrap,
            {100, 200, 300, 400, 0, 0, 0, 0},
            "SUM Test 2: with zeros");

    TestSum(cc, keyPair, QBFVINIT, PInput, Q, Bigq,
            numSlots, numSlotsCKKS, depth, levelsAvailableBeforeBootstrap,
            {5, 5, 5, 5, 5, 5, 5, 5},
            "SUM Test 3: all same");

    // ===== COUNT Tests =====
    TestCount(cc, keyPair, ep, QBFVINIT, PInput, POutput, Q, Bigq,
              scaleTHI, order, numSlots, ringDim, dcrtBits,
              numSlotsCKKS, depth, levelsAvailableBeforeBootstrap, levelsComputation,
              {10, 20, 30, 40, 50, 60, 70, 80},
              "COUNT Test 1: all non-zero");

    TestCount(cc, keyPair, ep, QBFVINIT, PInput, POutput, Q, Bigq,
              scaleTHI, order, numSlots, ringDim, dcrtBits,
              numSlotsCKKS, depth, levelsAvailableBeforeBootstrap, levelsComputation,
              {10, 0, 30, 0, 50, 0, 70, 0},
              "COUNT Test 2: half zero");

    TestCount(cc, keyPair, ep, QBFVINIT, PInput, POutput, Q, Bigq,
              scaleTHI, order, numSlots, ringDim, dcrtBits,
              numSlotsCKKS, depth, levelsAvailableBeforeBootstrap, levelsComputation,
              {0, 0, 0, 0, 0, 0, 0, 0},
              "COUNT Test 3: all zero");

    // ===== AVG Tests =====
    TestAvg(cc, keyPair, ep, QBFVINIT, PInput, POutput, Q, Bigq,
            scaleTHI, order, numSlots, ringDim, dcrtBits,
            numSlotsCKKS, depth, levelsAvailableBeforeBootstrap, levelsComputation,
            {10, 20, 30, 40, 50, 60, 70, 80},
            "AVG Test 1: simple");

    TestAvg(cc, keyPair, ep, QBFVINIT, PInput, POutput, Q, Bigq,
            scaleTHI, order, numSlots, ringDim, dcrtBits,
            numSlotsCKKS, depth, levelsAvailableBeforeBootstrap, levelsComputation,
            {100, 0, 200, 0, 300, 0, 400, 0},
            "AVG Test 2: with zeros");

    // ===== MAX Tests =====
    TestMax(cc, keyPair, ep, QBFVINIT, PInput, PDigit, Q, Bigq,
            scaleTHI, scaleStepTHI, order, numSlots, ringDim, dcrtBits,
            numSlotsCKKS, depth, levelsAvailableBeforeBootstrap, levelsComputation,
            {10, 50, 30, 80, 20, 60, 40, 70},
            "MAX Test 1: varied values");

    TestMax(cc, keyPair, ep, QBFVINIT, PInput, PDigit, Q, Bigq,
            scaleTHI, scaleStepTHI, order, numSlots, ringDim, dcrtBits,
            numSlotsCKKS, depth, levelsAvailableBeforeBootstrap, levelsComputation,
            {100, 200, 300, 400, 500, 600, 700, 800},
            "MAX Test 2: increasing");

    TestMax(cc, keyPair, ep, QBFVINIT, PInput, PDigit, Q, Bigq,
            scaleTHI, scaleStepTHI, order, numSlots, ringDim, dcrtBits,
            numSlotsCKKS, depth, levelsAvailableBeforeBootstrap, levelsComputation,
            {5, 5, 5, 5, 5, 5, 5, 5},
            "MAX Test 3: all same");

    // ===== MIN Tests =====
    TestMin(cc, keyPair, ep, QBFVINIT, PInput, PDigit, Q, Bigq,
            scaleTHI, scaleStepTHI, order, numSlots, ringDim, dcrtBits,
            numSlotsCKKS, depth, levelsAvailableBeforeBootstrap, levelsComputation,
            {10, 50, 30, 80, 20, 60, 40, 70},
            "MIN Test 1: varied values");

    TestMin(cc, keyPair, ep, QBFVINIT, PInput, PDigit, Q, Bigq,
            scaleTHI, scaleStepTHI, order, numSlots, ringDim, dcrtBits,
            numSlotsCKKS, depth, levelsAvailableBeforeBootstrap, levelsComputation,
            {800, 700, 600, 500, 400, 300, 200, 100},
            "MIN Test 2: decreasing");

    TestMin(cc, keyPair, ep, QBFVINIT, PInput, PDigit, Q, Bigq,
            scaleTHI, scaleStepTHI, order, numSlots, ringDim, dcrtBits,
            numSlotsCKKS, depth, levelsAvailableBeforeBootstrap, levelsComputation,
            {42, 42, 42, 42, 42, 42, 42, 42},
            "MIN Test 3: all same");

    std::cerr << "\n===== All Aggregation Tests Completed! =====" << std::endl;

    return 0;
}
