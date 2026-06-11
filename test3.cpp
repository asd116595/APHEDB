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




lbcrypto::Ciphertext<lbcrypto::DCRTPoly> SignCipher(std::vector<lbcrypto::Poly> &ctxtBFV, lbcrypto::CryptoContextCKKSRNS::ContextType &cc,
    BigInteger QBFVInit, BigInteger PInput, BigInteger PDigit, BigInteger Q, BigInteger Bigq, uint64_t scaleTHI, uint64_t scaleStepTHI, size_t order, uint32_t numSlots, uint32_t ringDim, uint32_t dcrtBits,
    std::shared_ptr<lbcrypto::M4DCRTParams> &ep, lbcrypto::KeyPair<lbcrypto::DCRTPoly> &keyPair, uint32_t numSlotsCKKS, uint32_t depth, uint32_t levelsAvailableBeforeBootstrap, uint32_t levelsComputation
);


int main() {

    //MultiValueBootstrapping(QBFVINIT, BigInteger(256), BigInteger(256), (BigInteger(1) << 47), (BigInteger(1) << 47),32, 1, 256, 2048, 1);
    std::cerr << "\n\n Homomorphically evaluate the sign." << std::endl << std::endl;
    BigInteger PInput(BigInteger(4096));
    BigInteger PDigit(BigInteger(2));
    BigInteger Q(BigInteger(1) << 46);
    BigInteger Bigq(BigInteger(1) << 35);
    uint64_t scaleTHI     = 1;
    uint64_t scaleStepTHI = 1;
    size_t order          = 1;
    uint32_t ringDim      = 2048;
    uint32_t numSlots     = 8;
    
    

    bool flagSP       = (numSlots <= ringDim / 2);  // sparse packing
    auto numSlotsCKKS = flagSP ? numSlots : numSlots / 2;

    //auto numSlotsCKKS = numSlots;

    uint32_t dcrtBits                       = Bigq.GetMSB() - 1;
    uint32_t firstMod                       = Bigq.GetMSB() - 1;
    uint32_t levelsAvailableAfterBootstrap  = 5;
    uint32_t levelsAvailableBeforeBootstrap = 0;
    uint32_t levelsComputation              = 5;
    uint32_t dnum                           = 3;
    SecretKeyDist secretKeyDist             = SPARSE_ENCAPSULATED;
    std::vector<uint32_t> lvlb              = {1,1};//{3, 3};


    auto a = PInput.ConvertToInt<int64_t>();
    auto b = PDigit.ConvertToInt<int64_t>();

    auto funcMod = [b](int64_t x) -> int64_t {
        return (x % b);
    };
    auto funcStep = [a, b](int64_t x) -> int64_t {
        return (x % a) >= (b / 2);
    };

    std::vector<int64_t> coeffintMod;
    std::vector<std::complex<double>> coeffcompMod;
    std::vector<std::complex<double>> coeffcompStep;

    bool binaryLUT = (PDigit.ConvertToInt() == 2) && (order == 1);


    if (binaryLUT) {
        coeffintMod = {funcMod(1), funcMod(0) - funcMod(1)};
    }
    else {
        coeffcompMod  = GetHermiteTrigCoefficients(funcMod, PDigit.ConvertToInt(), order, scaleTHI);  // divided by 2
        coeffcompStep = GetHermiteTrigCoefficients(funcStep, PDigit.ConvertToInt(), order,
                                                   scaleStepTHI);  // divided by 2
    }

    // ===== 2. CryptoContext =====
    CCParams<CryptoContextCKKSRNS> parameters;
    parameters.SetSecretKeyDist(secretKeyDist);
    parameters.SetSecurityLevel(HEStd_NotSet);
    parameters.SetScalingModSize(dcrtBits);
    parameters.SetScalingTechnique(FIXEDMANUAL);
    parameters.SetFirstModSize(firstMod);
    parameters.SetNumLargeDigits(dnum);
    std::cerr << "\n\n Homomorphically evaluate the sign." << std::endl << std::endl;
    parameters.SetBatchSize(numSlotsCKKS);
    parameters.SetRingDim(ringDim);

    uint32_t depth = levelsAvailableAfterBootstrap+levelsComputation;
    //uint32_t depth = levelsAvailableAfterBootstrap;

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

    // ===== 3. KeyGen =====
    auto keyPair = cc->KeyGen();

    std::cout << "CKKS scheme is using ring dimension " << cc->GetRingDimension() << " and a multiplicative depth of "
              << depth << std::endl
              << std::endl;

    cc->EvalMultKeyGen(keyPair.secretKey);

    if (binaryLUT)
        cc->EvalFBTSetup(coeffintMod, numSlotsCKKS, PDigit, PInput, Bigq, keyPair.publicKey, {0, 0}, lvlb,
                         levelsAvailableAfterBootstrap, levelsComputation, order);
    else
        cc->EvalFBTSetup(coeffcompMod, numSlotsCKKS, PDigit, PInput, Bigq, keyPair.publicKey, {0, 0}, lvlb,
                         levelsAvailableAfterBootstrap, levelsComputation, order);

    cc->EvalBootstrapKeyGen(keyPair.secretKey, numSlotsCKKS);

    cc->EvalAtIndexKeyGen(keyPair.secretKey, std::vector<int32_t>({-2}));

    // ===== 4. Encrypt 输入 =====
 
    std::cerr << "\n\n Homomorphically evaluate the sign." << std::endl << std::endl;

    std::vector<int64_t> x_input1 = {100,98,0,5,16,32,128,96};
    if (x_input1.size() < numSlots)
        x_input1 = Fill<int64_t>(x_input1, numSlots);
    std::vector<int64_t> x_input2 = {97};
    if (x_input2.size() < numSlots)
        x_input2 = Fill<int64_t>(x_input2, numSlots);
    std::vector<int64_t> x_input3 = {18};
    if (x_input3.size() < numSlots)
        x_input3 = Fill<int64_t>(x_input3, numSlots);



    //cc->EvalHomDecoding(c1, scaleTHI, 0);
    auto ep = SchemeletRLWEMP::GetElementParams(keyPair.secretKey, depth - (levelsAvailableBeforeBootstrap > 0));


    //bool flagBR = (lvlb[0] != 1 || lvlb[1] != 1);

    auto ctxtBFV1 = SchemeletRLWEMP::EncryptCoeff(x_input1, QBFVINIT, PInput, keyPair.secretKey, ep);

    auto ckks11 = SchemeletRLWEMP::ConvertRLWEToCKKS(*cc, ctxtBFV1, keyPair.publicKey, QBFVINIT, numSlotsCKKS,
                                                       depth-(levelsAvailableBeforeBootstrap > 0));
    std::cout<<"depth:"<<depth-(levelsAvailableBeforeBootstrap > 0)<<"\n";

    auto ctxtBFV2 = SchemeletRLWEMP::EncryptCoeff(x_input2, QBFVINIT, PInput, keyPair.secretKey, ep);
    auto ctxtBFV3 = SchemeletRLWEMP::EncryptCoeff(x_input3, QBFVINIT, PInput, keyPair.secretKey, ep);
   
    std::vector<lbcrypto::Poly> ctxtBFV4 = { ctxtBFV2[0] - ctxtBFV1[0], ctxtBFV2[1] - ctxtBFV1[1]};
    std::vector<lbcrypto::Poly> ctxtBFV5 = { ctxtBFV1[0] - ctxtBFV3[0], ctxtBFV1[1] - ctxtBFV3[1]};

    auto start1 = std::chrono::high_resolution_clock::now();

    auto ctxtSign11 = SignCipher(ctxtBFV4, cc, QBFVINIT, PInput, PDigit, Q, Bigq, scaleTHI, scaleStepTHI, order, numSlots, ringDim, dcrtBits, ep, keyPair, numSlotsCKKS, depth, levelsAvailableBeforeBootstrap, levelsComputation);
    

    
    auto ctxtSign12 = SignCipher(ctxtBFV5, cc, QBFVINIT, PInput, PDigit, Q, Bigq, scaleTHI, scaleStepTHI, order, numSlots, ringDim, dcrtBits, ep, keyPair, numSlotsCKKS, depth, levelsAvailableBeforeBootstrap, levelsComputation);
    
    std::vector<double> x_one = {1};
    if (x_one.size() < numSlots)
        x_one = Fill<double>(x_one, numSlots);

    Plaintext ptxtone = cc->MakeCKKSPackedPlaintext(x_one, 1, depth-6-levelsComputation, nullptr, numSlotsCKKS);

    
        
    auto cone = cc->Encrypt(keyPair.publicKey, ptxtone);
        //auto ctxtAfterFBT2 = cc->EvalSub(c3,ctxtStepResult);
    auto ctxtAfterFBT1 = cc->EvalSub(cone,ctxtSign11);
    auto ctxtAfterFBT2 = cc->EvalSub(cone,ctxtSign12);

    auto ctxtAfterFBT3 = cc->EvalMult(ctxtAfterFBT1,ctxtAfterFBT2);

    cc->ModReduceInPlace(ctxtAfterFBT3);
        
    auto ctxtAfterFBT4 = cc->EvalHomDecoding(ctxtAfterFBT1, scaleTHI, levelsComputation);

    auto ctxtResult1 = SchemeletRLWEMP::ConvertCKKSToRLWE(ctxtAfterFBT4, Q);



    auto computed1 =
        SchemeletRLWEMP::DecryptCoeff(ctxtResult1, Q, 2, keyPair.secretKey, ep, numSlotsCKKS, numSlots);


    std::cerr << "First 8 elements of the last sign: [";
    std::copy_n(computed1.begin(), numSlots, std::ostream_iterator<int64_t>(std::cerr, " "));
    std::cerr << "]" << std::endl;

    auto end1 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> elapsed1 = end1 - start1;

    std::cout << "程序运行时间1：" << elapsed1.count() << " ms" << std::endl;

    std::vector<lbcrypto::Poly> ctxtBFV6 = { ctxtBFV3[0] - ctxtBFV1[0], ctxtBFV3[1] - ctxtBFV1[1]};
    std::vector<lbcrypto::Poly> ctxtBFV7 = { ctxtBFV1[0] - ctxtBFV2[0], ctxtBFV1[1] - ctxtBFV2[1]};

    
    return 0;
}



lbcrypto::Ciphertext<lbcrypto::DCRTPoly> SignCipher(std::vector<lbcrypto::Poly> &ctxtBFV, lbcrypto::CryptoContextCKKSRNS::ContextType &cc,
    BigInteger QBFVInit, BigInteger PInput, BigInteger PDigit, BigInteger Q, BigInteger Bigq, uint64_t scaleTHI, uint64_t scaleStepTHI, size_t order, uint32_t numSlots, uint32_t ringDim, uint32_t dcrtBits,
    std::shared_ptr<lbcrypto::M4DCRTParams> &ep, lbcrypto::KeyPair<lbcrypto::DCRTPoly> &keyPair, uint32_t numSlotsCKKS, uint32_t depth, uint32_t levelsAvailableBeforeBootstrap, uint32_t levelsComputation
)
    {
    
    //auto cc = ctxtBFV->GetCryptoContext();
    auto encryptedOrg = ctxtBFV;
    auto a = PInput.ConvertToInt<int64_t>();
    auto b = PDigit.ConvertToInt<int64_t>();
    auto funcMod = [b](int64_t x) -> int64_t {
        return (x % b);
    };
    auto funcStep = [a, b](int64_t x) -> int64_t {
        return (x % a) >= (b / 2);
    };


    std::vector<int64_t> coeffintMod;
    std::vector<std::complex<double>> coeffcompMod;
    std::vector<std::complex<double>> coeffcompStep;

    bool binaryLUT = (PDigit.ConvertToInt() == 2) && (order == 1);


    if (binaryLUT) {
        coeffintMod = {funcMod(1), funcMod(0) - funcMod(1)};
    }
    else {
        coeffcompMod  = GetHermiteTrigCoefficients(funcMod, PDigit.ConvertToInt(), order, scaleTHI);  // divided by 2
        coeffcompStep = GetHermiteTrigCoefficients(funcStep, PDigit.ConvertToInt(), order,
                                                   scaleStepTHI);  // divided by 2
    }

    SchemeletRLWEMP::ModSwitch(ctxtBFV, Q, QBFVInit);
    uint32_t QBFVBits = Q.GetMSB() - 1;
    std::cout<<"QBFVBits:"<<QBFVBits<<"\n";
    /* 8. Set up the sign loop parameters. */
    std::vector<int64_t> coeffint;
    std::vector<std::complex<double>> coeffcomp;
    if (binaryLUT)
        coeffint = coeffintMod;
    else
        coeffcomp = coeffcompMod;

    const bool checkeq2       = PDigit.ConvertToInt() == 2;
    const bool checkgt2       = PDigit.ConvertToInt() > 2;
    const uint32_t pDigitBits = PDigit.GetMSB() - 1;
    std::cout<<"pDigitBits:"<<pDigitBits<<"\n";

    BigInteger QNew;
    BigInteger pOrig = PInput;
    BigInteger Qorig = Q; // remember original Q so we can restore before returning

    bool step                = false;
    bool go                  = QBFVBits > dcrtBits;
    //size_t levelsToDrop      = 0;
    uint32_t postScalingBits = 0;

    /* 9. Start the sign loop. For arbitrary digit size, pNew > 2, the last iteration needs
     * to evaluate step pNew not mod pNew.
     * Currently this only works when log(pNew) divides log(p).
    */
    Ciphertext<DCRTPoly> ctxtStepResult;
    Ciphertext<DCRTPoly> ctxtAfterFBT1;

    while (go) {
        auto encryptedDigit = ctxtBFV;

        /* 9.1. Apply mod Bigq to extract the digit and convert it from RLWE to CKKS. */
        encryptedDigit[0].SwitchModulus(Bigq, 1, 0, 0);
        encryptedDigit[1].SwitchModulus(Bigq, 1, 0, 0);

        auto ctxt = SchemeletRLWEMP::ConvertRLWEToCKKS(*cc, encryptedDigit, keyPair.publicKey, Bigq, numSlotsCKKS,
                                                       depth - (levelsAvailableBeforeBootstrap > 0));

        Ciphertext<DCRTPoly> ctxtAfterFBT;
        //Ciphertext<DCRTPoly> ctxtAfterFBT1;
        
        bool isLastIteration = step || (checkeq2 && (QBFVBits - pDigitBits <= dcrtBits));
        if (isLastIteration) {
            //std::cout<<"ctxt"<<ctxt.G;
            if (binaryLUT)
                ctxtAfterFBT = cc->EvalFBTNoDecoding(ctxt, coeffint, pDigitBits, ep->GetModulus(), order);
            else
                ctxtAfterFBT = cc->EvalFBTNoDecoding(ctxt, coeffcomp, pDigitBits, ep->GetModulus(), order);
            ctxtStepResult = ctxtAfterFBT;

            //ctxtStepResult = cc->EvalMult(ctxtAfterFBT, scaleTHI * (1 << postScalingBits));
            auto ctxtAfterFBT3 = cc->EvalHomDecoding(ctxtStepResult, scaleTHI * (1 << postScalingBits), levelsComputation);
            std::vector<lbcrypto::Poly> ctxtResult;
            ctxtResult = SchemeletRLWEMP::ConvertCKKSToRLWE(ctxtAfterFBT3, Q);
            auto computed11 =
                SchemeletRLWEMP::DecryptCoeff(ctxtResult, Q, PInput, keyPair.secretKey, ep, numSlotsCKKS, numSlots);
            std::cerr << "First 8 elements of the obtained digit in last iteration: [";
            std::copy_n(computed11.begin(), 8, std::ostream_iterator<int64_t>(std::cerr, " "));
            std::cerr << "]" << std::endl;  

            break;
            //ctxtAfterFBT1 = cc->EvalMult(ctxtAfterFBT, scaleTHI * (1 << postScalingBits));


        } else {
            if (binaryLUT)
                ctxtAfterFBT = cc->EvalFBT(ctxt, coeffint, pDigitBits, ep->GetModulus(), scaleTHI * (1 << postScalingBits),
                                           levelsComputation, order);
            else
                ctxtAfterFBT = cc->EvalFBT(ctxt, coeffcomp, pDigitBits, ep->GetModulus(), scaleTHI * (1 << postScalingBits),
                                           levelsComputation, order);
            auto polys = SchemeletRLWEMP::ConvertCKKSToRLWE(ctxtAfterFBT, Q);
            ctxtBFV[0] = ctxtBFV[0] - polys[0];
            ctxtBFV[1] = ctxtBFV[1] - polys[1];

            QNew       = Q >> pDigitBits;
            ctxtBFV[0] = ctxtBFV[0].MultiplyAndRound(QNew, Q);
            ctxtBFV[0].SwitchModulus(QNew, 1, 0, 0);
            ctxtBFV[1] = ctxtBFV[1].MultiplyAndRound(QNew, Q);
            ctxtBFV[1].SwitchModulus(QNew, 1, 0, 0);
            Q >>= pDigitBits;
            PInput >>= pDigitBits;
            QBFVBits -= pDigitBits;
            postScalingBits += pDigitBits;
            
        }

        go = QBFVBits > dcrtBits;
        if (checkgt2 && !go && !step) {
            if (!binaryLUT)
                coeffcomp = coeffcompStep;
            scaleTHI = scaleStepTHI;
            step     = true;
            go       = true;

        }
    }
        
    

        //ctxtStepResult = cc->EvalRotate(ctxtStepResult, -2);
/*
        std::vector<double> x_input11 = {2,3,2,2,2,2,2,6};
        std::vector<double> x_input22 = {1,0,1,0,1,1,1,0};
        std::vector<double> x_input33 = {1,1,1,1,1,1,1,1};

        Plaintext ptxt1 = cc->MakeCKKSPackedPlaintext(x_input11, 1, depth-6-levelsComputation, nullptr, numSlotsCKKS);
        Plaintext ptxt2 = cc->MakeCKKSPackedPlaintext(x_input22, 1, depth-6-levelsComputation, nullptr, numSlotsCKKS);
        Plaintext ptxt3 = cc->MakeCKKSPackedPlaintext(x_input33, 1, depth-6-levelsComputation, nullptr, numSlotsCKKS);

        auto c1 = cc->Encrypt(keyPair.publicKey, ptxt1);
        
        auto c3 = cc->Encrypt(keyPair.publicKey, ptxt3);
        //auto ctxtAfterFBT2 = cc->EvalSub(c3,ctxtStepResult);
        auto ctxtAfterFBT2 = cc->EvalMult(ctxtStepResult, c1);

        cc->ModReduceInPlace(ctxtAfterFBT2);
        
        auto ctxtAfterFBT4 = cc->EvalHomDecoding(ctxtAfterFBT2, scaleTHI, levelsComputation-1);
        auto ctxtResult1 = SchemeletRLWEMP::ConvertCKKSToRLWE(ctxtAfterFBT4, Qorig);
        auto computed1 =
            SchemeletRLWEMP::DecryptCoeff(ctxtResult1, Qorig, pOrig, keyPair.secretKey, ep, numSlotsCKKS, numSlots);

        std::cerr << "First 8 elements of the last sign: [";
        std::copy_n(computed1.begin(), 8, std::ostream_iterator<int64_t>(std::cerr, " "));
        std::cerr << "]" << std::endl;
*/



    return ctxtStepResult;
    //return ctxtAfterFBT2;
}



