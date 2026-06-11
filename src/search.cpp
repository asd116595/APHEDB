#include "search.h"
#include "equal_bootstrapping.h"

using namespace lbcrypto;

lbcrypto::Ciphertext<lbcrypto::DCRTPoly> Search(
    std::vector<lbcrypto::Poly>& ctxtStr,
    std::vector<std::vector<lbcrypto::Poly>> ctxtstr,
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

    const size_t numRotations = ctxtstr.size();
    lbcrypto::Ciphertext<lbcrypto::DCRTPoly> result;
    // 对每个旋转版本计算差值并检查相等性
    for (size_t i = 0; i < numRotations; ++i) {
        // 创建差值密文: ctxtStr - ctxtstr[i]
        std::vector<lbcrypto::Poly> diff = ctxtStr;
        diff[0] -= ctxtstr[i][0];
        diff[1] -= ctxtstr[i][1];

        
        auto eq_check = EqualBootstrapping(diff, cc, QBFVInit, PInput, POutput, Q, Bigq, 
            scaleTHI, order, numSlots, ringDim, dcrtBits, ep, keyPair, 
            numSlotsCKKS, depth, levelsAvailableBeforeBootstrap, levelsComputation);
        lbcrypto::Ciphertext<lbcrypto::DCRTPoly> rotated_result;
        lbcrypto::Ciphertext<lbcrypto::DCRTPoly> rotated_sum;
        for (size_t j = 0; j < numRotations; ++j) {
            
            if (j == 0) {
                rotated_sum = eq_check;
                rotated_result = eq_check;
            } else {
                rotated_result = cc->EvalRotate(rotated_result, 1);
                // 确保模数链对齐后再做加法
                rotated_sum = cc->EvalAdd(rotated_sum, rotated_result);
            }

        }

         auto ctxt1 = cc->EvalHomDecoding(rotated_sum, scaleTHI, levelsComputation);
    
        auto polys1 = SchemeletRLWEMP::ConvertCKKSToRLWE(ctxt1, Q);

        auto computed1 = SchemeletRLWEMP::DecryptCoeff(polys1, Q, POutput, keyPair.secretKey, ep, numSlotsCKKS, numSlots);

        std::cerr << "rotated_sum " << numSlots << " elements of the obtained output % PInput: [";
        std::copy_n(computed1.begin(), numSlots, std::ostream_iterator<int64_t>(std::cerr, " "));
        std::cerr << "]" << std::endl;

        std::vector<double> mask(numRotations, 0);
        mask[i] = 1; 
        if (mask.size() < numSlots)
            mask = Fill<double>(mask, numSlots);

        std::cout<<mask<<std::endl;

        Plaintext ptxt = cc->MakeCKKSPackedPlaintext(mask, 1, rotated_sum->GetLevel(), nullptr, numSlotsCKKS);
        auto c1 = cc->Encrypt(keyPair.publicKey, ptxt);
        auto masked_eq_check = cc->EvalMult(rotated_sum, c1);
        cc->ModReduceInPlace(masked_eq_check);

       


        if (i == 0) {
            result = masked_eq_check;
        } else {
            // 确保模数链对齐后再做加法
            result = cc->EvalAdd(result, masked_eq_check);
        }
        std::cerr << "Depth/Level: " << result->GetLevel() << std::endl;
    }
    std::cerr << "=== Search result info ===" << std::endl;
    
    auto ctxtAfter = cc->EvalHomDecoding(result, scaleTHI, levelsComputation-1);
    
    auto polys = SchemeletRLWEMP::ConvertCKKSToRLWE(ctxtAfter, Q);

    auto computed = SchemeletRLWEMP::DecryptCoeff(polys, Q, POutput, keyPair.secretKey, ep, numSlotsCKKS, numSlots);

    std::cerr << "First " << numSlots << " elements of the obtained output % PInput: [";
    std::copy_n(computed.begin(), numSlots, std::ostream_iterator<int64_t>(std::cerr, " "));
    std::cerr << "]" << std::endl;

    return result;
}
