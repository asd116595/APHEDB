#define PROFILE
#include "homlike.h"

#include <iostream>
#include <sstream>
#include <iomanip>

using namespace lbcrypto;

int main() {
    BigInteger PInput(BigInteger(4096));
    BigInteger POutput(BigInteger(4096));
    BigInteger Q(BigInteger(1) << 46);
    BigInteger Bigq(BigInteger(1) << 35);
    uint64_t scaleTHI = 1;
    size_t order = 1;
    uint32_t ringDim = 2048;
    uint32_t numSlots = 1024;   // 1024 slots
    uint32_t lmax = 8;           // 8 chars per name
    
    bool flagSP = (numSlots <= ringDim / 2);
    auto numSlotsCKKS = flagSP ? numSlots : numSlots / 2;

    uint32_t dcrtBits = Bigq.GetMSB() - 1;
    uint32_t firstMod = Bigq.GetMSB() - 1;
    uint32_t levelsAvailableAfterBootstrap = 5;
    uint32_t levelsAvailableBeforeBootstrap = 0;
    uint32_t levelsComputation = 5;
    uint32_t dnum = 3;
    SecretKeyDist secretKeyDist = SPARSE_ENCAPSULATED;
    std::vector<uint32_t> lvlb = {1, 1};

    auto b = PInput.ConvertToInt<int64_t>();
    auto funcMod = [b](int64_t x) -> int64_t {
        return (x % b);
    };

    std::vector<int64_t> coeffintMod;
    std::vector<std::complex<double>> coeffcompMod;
    bool binaryLUT = (PInput.ConvertToInt() == 2) && (order == 1);

    if (binaryLUT) {
        coeffintMod = {funcMod(1), funcMod(0) - funcMod(1), 0};
    } else {
        coeffcompMod = GetHermiteTrigCoefficients(funcMod, PInput.ConvertToInt(), order, scaleTHI);
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
        depth += FHECKKSRNS::GetFBTDepth(lvlb, coeffintMod, Bigq, order, secretKeyDist);
    else
        depth += FHECKKSRNS::GetFBTDepth(lvlb, coeffcompMod, Bigq, order, secretKeyDist);
    parameters.SetMultiplicativeDepth(depth);

    auto cc = GenCryptoContext(parameters);
    cc->Enable(PKE);
    cc->Enable(KEYSWITCH);
    cc->Enable(LEVELEDSHE);
    cc->Enable(ADVANCEDSHE);
    cc->Enable(FHE);

    auto keyPair = cc->KeyGen();
    std::cout << "CKKS ring dim: " << cc->GetRingDimension() << ", depth: " << depth << std::endl;

    cc->EvalMultKeyGen(keyPair.secretKey);

    if (binaryLUT)
        cc->EvalFBTSetup(coeffintMod, numSlotsCKKS, Bigq, PInput, Bigq, keyPair.publicKey, {0, 0}, lvlb,
                         levelsAvailableAfterBootstrap, levelsComputation, order);
    else
        cc->EvalFBTSetup(coeffcompMod, numSlotsCKKS, Bigq, PInput, Bigq, keyPair.publicKey, {0, 0}, lvlb,
                         levelsAvailableAfterBootstrap, levelsComputation, order);

    cc->EvalBootstrapKeyGen(keyPair.secretKey, numSlotsCKKS);
    cc->EvalAtIndexKeyGen(keyPair.secretKey, std::vector<int32_t>({-2}));

    // ===== Test Data: 128 names, each 8 chars =====
    std::vector<std::string> names;
    for (size_t i = 0; i < 128; ++i) {
        std::ostringstream oss;
        oss << "Name" << std::setw(3) << std::setfill('0') << i;
        std::string name = oss.str();
        while (name.length() < 8) name += ' ';
        names.push_back(name);
    }
    names[0] = "Alice    ";  // First name is Alice (8 chars)
    names[8] = "Bob      ";  // 9th name is Bob (8 chars)

    // Combine all names (128 * 8 = 1024 chars)
    std::string combined;
    for (size_t i = 0; i < numSlots / lmax && i < names.size(); ++i) {
        combined += names[i];
    }
    while (combined.size() < numSlots) combined += ' ';

    std::vector<int64_t> asciitext;
    for (char c : combined) {
        asciitext.push_back(static_cast<int>(c));
    }

    std::cout << numSlots << " chars, " << (numSlots / lmax) << " names" << std::endl;

    auto ep = SchemeletRLWEMP::GetElementParams(
        keyPair.secretKey, depth - (levelsAvailableBeforeBootstrap > 0));

    auto ctxtStr = SchemeletRLWEMP::EncryptCoeff(
        asciitext, QBFVINIT, PInput, keyPair.secretKey, ep);

    std::string core = "Alice";

    // ===== Case 1: Direct =====
    std::cout << "\n=== Case 1: Direct ===" << std::endl;
    RunCase1_Direct(ctxtStr, core, lmax, cc, QBFVINIT, PInput, POutput, Q, Bigq,
                    scaleTHI, order, numSlots, ringDim, dcrtBits, ep, keyPair,
                    numSlotsCKKS, depth, levelsAvailableBeforeBootstrap, levelsComputation);

    // ===== Case 2: _ in front =====
    std::cout << "\n=== Case 2: _ in front ===" << std::endl;
    RunCase2_UnderFront(ctxtStr, core, lmax, cc, QBFVINIT, PInput, POutput, Q, Bigq,
                        scaleTHI, order, numSlots, ringDim, dcrtBits, ep, keyPair,
                        numSlotsCKKS, depth, levelsAvailableBeforeBootstrap, levelsComputation);

    // ===== Case 3: _ in middle =====
    std::cout << "\n=== Case 3: _ in middle ===" << std::endl;
    RunCase3_UnderMid(ctxtStr, core, lmax, cc, QBFVINIT, PInput, POutput, Q, Bigq,
                      scaleTHI, order, numSlots, ringDim, dcrtBits, ep, keyPair,
                      numSlotsCKKS, depth, levelsAvailableBeforeBootstrap, levelsComputation);

    // ===== Case 4: _ in back =====
    std::cout << "\n=== Case 4: _ in back ===" << std::endl;
    RunCase4_UnderBack(ctxtStr, core, lmax, cc, QBFVINIT, PInput, POutput, Q, Bigq,
                       scaleTHI, order, numSlots, ringDim, dcrtBits, ep, keyPair,
                       numSlotsCKKS, depth, levelsAvailableBeforeBootstrap, levelsComputation);

    // ===== Case 5: % in front =====
    std::cout << "\n=== Case 5: % in front ===" << std::endl;
    RunCase5_PercentFront(ctxtStr, core, lmax, cc, QBFVINIT, PInput, POutput, Q, Bigq,
                          scaleTHI, order, numSlots, ringDim, dcrtBits, ep, keyPair,
                          numSlotsCKKS, depth, levelsAvailableBeforeBootstrap, levelsComputation);

    // ===== Case 6: % in middle =====
    std::cout << "\n=== Case 6: % in middle ===" << std::endl;
    RunCase6_PercentMid(ctxtStr, core, lmax, cc, QBFVINIT, PInput, POutput, Q, Bigq,
                        scaleTHI, order, numSlots, ringDim, dcrtBits, ep, keyPair,
                        numSlotsCKKS, depth, levelsAvailableBeforeBootstrap, levelsComputation);

    // ===== Case 7: % in back =====
    std::cout << "\n=== Case 7: % in back ===" << std::endl;
    RunCase7_PercentBack(ctxtStr, core, lmax, cc, QBFVINIT, PInput, POutput, Q, Bigq,
                         scaleTHI, order, numSlots, ringDim, dcrtBits, ep, keyPair,
                         numSlotsCKKS, depth, levelsAvailableBeforeBootstrap, levelsComputation);

    std::cout << "\nAll tests completed!" << std::endl;
    return 0;
}
