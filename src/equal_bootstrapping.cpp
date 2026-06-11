#define PROFILE
#include "equal_bootstrapping.h"

using namespace lbcrypto;

const BigInteger QBFVINIT(BigInteger(1) << 60);
const BigInteger QBFVINITLARGE(BigInteger(1) << 80);

lbcrypto::Ciphertext<lbcrypto::DCRTPoly> EqualBootstrapping(
    std::vector<lbcrypto::Poly>& ctxtBFV, 
    lbcrypto::CryptoContextCKKSRNS::ContextType& cc,
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
    std::shared_ptr<lbcrypto::M4DCRTParams>& ep, 
    lbcrypto::KeyPair<lbcrypto::DCRTPoly>& keyPair, 
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
        coeffint = {
            func(1),
            func(0) - func(1),
            0};  // those are coefficients for [1, cos^2(pi x)], not [1, cos(2pi x)] as in the general case.
    }
    else {
        coeffcomp = GetHermiteTrigCoefficients(func, PInput.ConvertToInt(), order, scaleTHI);  // divided by 2
    }

    SchemeletRLWEMP::ModSwitch(ctxtBFV, Q, QBFVInit);

    auto encryptedDigit = ctxtBFV;
    encryptedDigit[0].SwitchModulus(Bigq, 1, 0, 0);
    encryptedDigit[1].SwitchModulus(Bigq, 1, 0, 0);

    auto ctxt = SchemeletRLWEMP::ConvertRLWEToCKKS(
        *cc, encryptedDigit, keyPair.publicKey, Bigq, numSlotsCKKS,
        depth - (levelsAvailableBeforeBootstrap > 0));

    Ciphertext<DCRTPoly> ctxtAfterFBTNoDecoding;
    if (binaryLUT)
        ctxtAfterFBTNoDecoding = cc->EvalFBTNoDecoding(ctxt, coeffint, PInput.GetMSB() - 1, ep->GetModulus(), order);
    else
        ctxtAfterFBTNoDecoding = cc->EvalFBTNoDecoding(ctxt, coeffcomp, PInput.GetMSB() - 1, ep->GetModulus(), order);

    auto ctxtAfterFBTDecoding = cc->EvalHomDecoding(ctxtAfterFBTNoDecoding, scaleTHI, levelsComputation);
    
    // Print modulus and depth of eq_check
    

    
    std::cerr << "=== eq_check info ===" << std::endl;
    std::cerr << "Depth/Level: " << ctxtAfterFBTDecoding->GetLevel() << std::endl;
    std::cerr << "=====================" << std::endl;
    auto polys = SchemeletRLWEMP::ConvertCKKSToRLWE(ctxtAfterFBTDecoding, Q);

    auto computed = SchemeletRLWEMP::DecryptCoeff(polys, Q, POutput, keyPair.secretKey, ep, numSlotsCKKS, numSlots);

    std::cerr << "First " << numSlots << " elements of the obtained output % PInput: [";
    std::copy_n(computed.begin(), numSlots, std::ostream_iterator<int64_t>(std::cerr, " "));
    std::cerr << "]" << std::endl;
    
    


    return ctxtAfterFBTNoDecoding;
}
