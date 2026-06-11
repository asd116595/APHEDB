#define PROFILE
#include "equal_bootstrapping.h"

#include <iostream>

using namespace lbcrypto;

int main() {
    BigInteger PInput(BigInteger(4096));
    BigInteger POutput(BigInteger(4096));
    BigInteger Q(BigInteger(1) << 47);
    BigInteger Bigq(BigInteger(1) << 47);
    uint64_t scaleTHI     = 32;
    size_t order          = 1;
    uint32_t ringDim      = 4096;
    uint32_t numSlots     = 8;
    
    bool flagSP       = (numSlots <= ringDim / 2);  // sparse packing
    auto numSlotsCKKS = flagSP ? numSlots : numSlots / 2;

    uint32_t dcrtBits                       = Bigq.GetMSB() - 1;
    uint32_t firstMod                       = Bigq.GetMSB() - 1;
    uint32_t levelsAvailableAfterBootstrap  = 5;
    uint32_t levelsAvailableBeforeBootstrap = 0;
    uint32_t levelsComputation              = 5;
    uint32_t dnum                           = 3;
    SecretKeyDist secretKeyDist             = SPARSE_ENCAPSULATED;
    std::vector<uint32_t> lvlb              = {1,1};

    auto a = PInput.ConvertToInt<int64_t>();
    auto func = [a](int64_t x) -> int64_t {
        return (x % a) == 0 ? 1 : 0;
    };

    std::vector<int64_t> coeffint;
    std::vector<std::complex<double>> coeffcomp;
    bool binaryLUT = (PInput.ConvertToInt() == 2) && (order == 1);

    if (binaryLUT) {
        coeffint = {
            func(1),
            func(0) - func(1),
            0};
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

    if (binaryLUT) {
        cc->EvalFBTSetup(coeffint, numSlotsCKKS, PInput, POutput, Bigq, keyPair.publicKey, {0, 0}, lvlb,
                         levelsAvailableAfterBootstrap, levelsComputation, order);
    }
    else {
        cc->EvalFBTSetup(coeffcomp, numSlotsCKKS, PInput, POutput, Bigq, keyPair.publicKey, {0, 0}, lvlb,
                         levelsAvailableAfterBootstrap, levelsComputation, order);
    }

    cc->EvalBootstrapKeyGen(keyPair.secretKey, numSlotsCKKS);
    cc->EvalMultKeyGen(keyPair.secretKey);

    std::vector<int64_t> x_input1 = {100, 98, 0, 5, 16, 32, 128, 96};
    if (x_input1.size() < numSlots)
        x_input1 = Fill<int64_t>(x_input1, numSlots);

    auto ep = SchemeletRLWEMP::GetElementParams(keyPair.secretKey, depth - (levelsAvailableBeforeBootstrap > 0));

    auto ctxtBFV = SchemeletRLWEMP::EncryptCoeff(x_input1, QBFVINIT, PInput, keyPair.secretKey, ep);

    auto ctxt = EqualBootstrapping(ctxtBFV, cc, QBFVINIT, PInput, POutput, Q, Bigq, scaleTHI, order, numSlots,
                                   ringDim, dcrtBits, ep, keyPair, numSlotsCKKS, depth, levelsAvailableBeforeBootstrap,
                                   levelsComputation);
   
    return 0;
}
