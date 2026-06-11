#define PROFILE
#include "binfhecontext.h"

#include "math/hermite.h"
#include "openfhe.h"
#include "schemelet/rlwe-mp.h"

#include <functional>
#include <chrono>
#include <cmath>

using namespace lbcrypto;

const BigInteger QBFVINIT(BigInteger(1) << 60);
const BigInteger QBFVINITLARGE(BigInteger(1) << 80);

// 修改后的ArbitraryLUT函数
lbcrypto::Ciphertext<lbcrypto::DCRTPoly> ArbitraryLUT(
    std::vector<lbcrypto::Poly>& ctxtBFV, 
    lbcrypto::CryptoContextCKKSRNS::ContextType& cc,
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
    std::shared_ptr<lbcrypto::M4DCRTParams>& ep, 
    lbcrypto::KeyPair<lbcrypto::DCRTPoly>& keyPair, 
    uint32_t numSlotsCKKS, 
    uint32_t depth, 
    uint32_t levelsAvailableBeforeBootstrap, 
    uint32_t levelsComputation) {
    
    /* 1. 定义LUT函数 */
    auto a = PInput.ConvertToInt<int64_t>();
    auto func = [a](int64_t x) -> int64_t {
        return (x % a) == 0 ? 1 : 0;
    };
    
    /* 2. 计算系数 */
    std::vector<int64_t> coeffint;
    std::vector<std::complex<double>> coeffcomp;
    bool binaryLUT = (PInput.ConvertToInt() == 2) && (order == 1);
    
    if (binaryLUT) {
        coeffint = {
            func(1),
            func(0) - func(1)
        };
    }
    else {
        coeffcomp = GetHermiteTrigCoefficients(func, PInput.ConvertToInt(), order, scaleTHI);
    }
    
    /* 3. 设置LUT */
    std::vector<uint32_t> lvlb = {3, 3};
    uint32_t levelsAvailableAfterBootstrap = 0;
    
    if (binaryLUT)
        cc->EvalFBTSetup(coeffint, numSlotsCKKS, PInput, PDigit, Bigq, keyPair.publicKey, {0, 0}, lvlb,
                         levelsAvailableAfterBootstrap, 0, order);
    else
        cc->EvalFBTSetup(coeffcomp, numSlotsCKKS, PInput, PDigit, Bigq, keyPair.publicKey, {0, 0}, lvlb,
                         levelsAvailableAfterBootstrap, 0, order);
    
    /* 4. 转换RLWE到CKKS */
    auto ctxt = SchemeletRLWEMP::ConvertRLWEToCKKS(*cc, ctxtBFV, keyPair.publicKey, Bigq, numSlotsCKKS,
                                                   depth - (levelsAvailableBeforeBootstrap > 0));
    
    /* 5. 应用LUT */
    Ciphertext<DCRTPoly> ctxtAfterFBT;
    if (binaryLUT)
        ctxtAfterFBT = cc->EvalFBT(ctxt, coeffint, PInput.GetMSB() - 1, ep->GetModulus(), scaleTHI, 0, order);
    else
        ctxtAfterFBT = cc->EvalFBT(ctxt, coeffcomp, PInput.GetMSB() - 1, ep->GetModulus(), scaleTHI, 0, order);
    
    /* 6. 返回结果 */
    return ctxtAfterFBT;
}

int main() {
    // 参数设置
    BigInteger PInput(BigInteger(4096));
    BigInteger PDigit(BigInteger(2));
    BigInteger Q(BigInteger(1) << 33);
    BigInteger Bigq(BigInteger(1) << 33);
    uint64_t scaleTHI = 1;
    uint64_t scaleStepTHI = 1;
    size_t order = 1;
    uint32_t numSlots = 8;
    uint32_t ringDim = 4096;
    uint32_t dcrtBits = Bigq.GetMSB() - 1;
    uint32_t levelsAvailableBeforeBootstrap = 0;
    uint32_t levelsComputation = 1;
    
    /* 1. 确定打包方式 */
    bool flagSP = (numSlots <= ringDim / 2);
    auto numSlotsCKKS = flagSP ? numSlots : numSlots / 2;
    
    /* 2. 准备要加密的数据 */
    std::vector<int64_t> plaintextData = {
        (PInput.ConvertToInt<int64_t>() / 2), 
        (PInput.ConvertToInt<int64_t>() / 2) + 1, 
        0, 3, 16, 33, 64,
        (PInput.ConvertToInt<int64_t>() - 1)
    };
    
    // 填充到指定大小
    if (plaintextData.size() < numSlots) {
        while (plaintextData.size() < numSlots) {
            plaintextData.push_back(plaintextData.back());
        }
    }
    
    std::cerr << "Plaintext data to encrypt: ";
    std::copy_n(plaintextData.begin(), 8, std::ostream_iterator<int64_t>(std::cerr, " "));
    std::cerr << std::endl;
    
    /* 3. 设置加密参数 */
    SecretKeyDist secretKeyDist = SPARSE_ENCAPSULATED;
    
    CCParams<CryptoContextCKKSRNS> parameters;
    parameters.SetSecretKeyDist(secretKeyDist);
    parameters.SetSecurityLevel(HEStd_NotSet);
    parameters.SetScalingModSize(dcrtBits);
    parameters.SetScalingTechnique(FIXEDMANUAL);
    parameters.SetFirstModSize(dcrtBits);
    parameters.SetNumLargeDigits(3);
    parameters.SetBatchSize(numSlotsCKKS);
    parameters.SetRingDim(ringDim);
    
    uint32_t depth = levelsComputation;
    parameters.SetMultiplicativeDepth(depth);
    
    auto cc = GenCryptoContext(parameters);
    cc->Enable(PKE);
    cc->Enable(KEYSWITCH);
    cc->Enable(LEVELEDSHE);
    cc->Enable(ADVANCEDSHE);
    cc->Enable(FHE);
    
    std::cout << "CKKS scheme is using ring dimension " << cc->GetRingDimension() 
              << " and a multiplicative depth of " << depth << std::endl << std::endl;
    
    /* 4. 生成密钥 */
    auto keyPair = cc->KeyGen();
    cc->EvalMultKeyGen(keyPair.secretKey);
    cc->EvalBootstrapKeyGen(keyPair.secretKey, numSlotsCKKS);
    
    /* 5. 获取RLWE参数 */
    auto ep = SchemeletRLWEMP::GetElementParams(keyPair.secretKey, depth - (levelsAvailableBeforeBootstrap > 0));
    
    /* 6. 加密数据 */
    auto ctxtBFV = SchemeletRLWEMP::EncryptCoeff(plaintextData, QBFVINIT, PInput, keyPair.secretKey, ep);
    SchemeletRLWEMP::ModSwitch(ctxtBFV, Q, QBFVINIT);
    
    /* 7. 调用ArbitraryLUT */
    auto ctxtAfterFBT = ArbitraryLUT(ctxtBFV, cc, QBFVINIT, PInput, PDigit, Q, Bigq, 
                                     scaleTHI, scaleStepTHI, order, numSlots, ringDim, 
                                     dcrtBits, ep, keyPair, numSlotsCKKS, depth, 
                                     levelsAvailableBeforeBootstrap, levelsComputation);
    
    /* 8. 转换回RLWE并解密 */
    auto polys = SchemeletRLWEMP::ConvertCKKSToRLWE(ctxtAfterFBT, Q);
    auto result = SchemeletRLWEMP::DecryptCoeff(polys, Q, PDigit, keyPair.secretKey, ep, numSlotsCKKS, numSlots);
    
    /* 9. 输出结果 */
    std::cerr << "Decrypted result (first 8 elements): ";
    std::copy_n(result.begin(), 8, std::ostream_iterator<int64_t>(std::cerr, " "));
    std::cerr << std::endl;
    
    /* 10. 验证结果 */
    auto func = [a = PInput.ConvertToInt<int64_t>()](int64_t x) -> int64_t {
        return (x % a) == 0 ? 1 : 0;
    };
    
    std::vector<int64_t> expected;
    for (size_t i = 0; i < std::min(plaintextData.size(), result.size()); i++) {
        expected.push_back(func(plaintextData[i]));
    }
    
    std::cerr << "Expected result: ";
    std::copy_n(expected.begin(), 8, std::ostream_iterator<int64_t>(std::cerr, " "));
    std::cerr << std::endl;
    
    // 计算误差
    std::vector<int64_t> errors;
    for (size_t i = 0; i < result.size() && i < expected.size(); i++) {
        errors.push_back(std::abs(result[i] - expected[i]));
    }
    auto max_error = std::max_element(errors.begin(), errors.end());
    std::cerr << "Max absolute error: " << *max_error << std::endl;
    
    return 0;
}