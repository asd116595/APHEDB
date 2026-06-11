#define PROFILE
#include "comparsion.h"

using namespace lbcrypto;

lbcrypto::Ciphertext<lbcrypto::DCRTPoly> SignCipher(
    std::vector<lbcrypto::Poly> ctxtBFV,   // by value - we modify a local copy
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
        coeffcompMod  = GetHermiteTrigCoefficients(funcMod, PDigit.ConvertToInt(), order, scaleTHI);
        coeffcompStep = GetHermiteTrigCoefficients(funcStep, PDigit.ConvertToInt(), order, scaleStepTHI);
    }

    SchemeletRLWEMP::ModSwitch(ctxtBFV, Q, QBFVInit);
    uint32_t QBFVBits = Q.GetMSB() - 1;

    std::vector<int64_t> coeffint;
    std::vector<std::complex<double>> coeffcomp;
    if (binaryLUT)
        coeffint = coeffintMod;
    else
        coeffcomp = coeffcompMod;

    const bool checkeq2       = PDigit.ConvertToInt() == 2;
    const bool checkgt2       = PDigit.ConvertToInt() > 2;
    const uint32_t pDigitBits = PDigit.GetMSB() - 1;

    BigInteger QNew;
    BigInteger Qorig = Q;

    bool step                = false;
    bool go                  = QBFVBits > dcrtBits;
    uint32_t postScalingBits = 0;

    Ciphertext<DCRTPoly> ctxtStepResult;

    while (go) {
        auto encryptedDigit = ctxtBFV;

        encryptedDigit[0].SwitchModulus(Bigq, 1, 0, 0);
        encryptedDigit[1].SwitchModulus(Bigq, 1, 0, 0);

        auto ctxt = SchemeletRLWEMP::ConvertRLWEToCKKS(*cc, encryptedDigit, keyPair.publicKey, Bigq, numSlotsCKKS,
                                                       depth - (levelsAvailableBeforeBootstrap > 0));

        Ciphertext<DCRTPoly> ctxtAfterFBT;

        bool isLastIteration = step || (checkeq2 && (QBFVBits - pDigitBits <= dcrtBits));
        if (isLastIteration) {
            if (binaryLUT)
                ctxtAfterFBT = cc->EvalFBTNoDecoding(ctxt, coeffint, pDigitBits, ep->GetModulus(), order);
            else
                ctxtAfterFBT = cc->EvalFBTNoDecoding(ctxt, coeffcomp, pDigitBits, ep->GetModulus(), order);
            ctxtStepResult = ctxtAfterFBT;

            // Debug: decrypt to verify
            auto ctxtAfterFBT3 = cc->EvalHomDecoding(ctxtStepResult, scaleTHI * (1 << postScalingBits), levelsComputation);
            std::vector<lbcrypto::Poly> ctxtResult;
            ctxtResult = SchemeletRLWEMP::ConvertCKKSToRLWE(ctxtAfterFBT3, Q);
            auto computed11 =
                SchemeletRLWEMP::DecryptCoeff(ctxtResult, Q, PInput, keyPair.secretKey, ep, numSlotsCKKS, numSlots);
            std::cerr << "SignCipher last iteration - first 8 elements: [";
            std::copy_n(computed11.begin(), 8, std::ostream_iterator<int64_t>(std::cerr, " "));
            std::cerr << "]" << std::endl;

            break;
        }
        else {
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

    return ctxtStepResult;
}


// ===== Comparison Wrapper Functions =====
// These all call SignCipher internally.
// sign(diff) = 1 if diff > 0, else 0

lbcrypto::Ciphertext<lbcrypto::DCRTPoly> GreaterThan(
    const std::vector<lbcrypto::Poly>& ctxtX,
    const std::vector<lbcrypto::Poly>& ctxtY,
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


    std::vector<lbcrypto::Poly> diff = {
        ctxtY[0] - ctxtX[0],
        ctxtY[1] - ctxtX[1]
    };

    return SignCipher(diff, cc, QBFVInit, PInput, PDigit, Q, Bigq,
                      scaleTHI, scaleStepTHI, order, numSlots, ringDim, dcrtBits,
                      ep, keyPair, numSlotsCKKS, depth,
                      levelsAvailableBeforeBootstrap, levelsComputation);
}


lbcrypto::Ciphertext<lbcrypto::DCRTPoly> GreaterEqual(
    const std::vector<lbcrypto::Poly>& ctxtX,
    const std::vector<lbcrypto::Poly>& ctxtY,
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


    std::vector<lbcrypto::Poly> diff = {
        ctxtX[0] - ctxtY[0],
        ctxtX[1] - ctxtY[1]
    };

    auto signResult = SignCipher(diff, cc, QBFVInit, PInput, PDigit, Q, Bigq,
                                  scaleTHI, scaleStepTHI, order, numSlots, ringDim, dcrtBits,
                                  ep, keyPair, numSlotsCKKS, depth,
                                  levelsAvailableBeforeBootstrap, levelsComputation);

    // Compute 1 - sign(y - x)
    // Use a fresh encryption of 1 at the same level as signResult
    std::vector<double> x_one = {1};
    if (x_one.size() < numSlots)
        x_one = Fill<double>(x_one, numSlots);

    Plaintext ptxtone = cc->MakeCKKSPackedPlaintext(x_one, 1, depth-6-levelsComputation, nullptr, numSlotsCKKS);
    auto cone = cc->Encrypt(keyPair.publicKey, ptxtone);

    auto ctxtAfterFBT4 = cc->EvalHomDecoding(cc->EvalSub(cone, signResult), scaleTHI, levelsComputation);
    auto ctxtResult1 = SchemeletRLWEMP::ConvertCKKSToRLWE(ctxtAfterFBT4, Q);
    auto computed1 =
            SchemeletRLWEMP::DecryptCoeff(ctxtResult1, Q, PInput, keyPair.secretKey, ep, numSlotsCKKS, numSlots);

        std::cerr << "First 8 elements of the last sign: [";
        std::copy_n(computed1.begin(), 8, std::ostream_iterator<int64_t>(std::cerr, " "));
        std::cerr << "]" << std::endl;


    return cc->EvalSub(cone, signResult);
}


lbcrypto::Ciphertext<lbcrypto::DCRTPoly> LessThan(
    const std::vector<lbcrypto::Poly>& ctxtX,
    const std::vector<lbcrypto::Poly>& ctxtY,
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

    // (i.e., x < y)
    std::vector<lbcrypto::Poly> diff = {
        ctxtX[0] - ctxtY[0],
        ctxtX[1] - ctxtY[1]
    };

    return SignCipher(diff, cc, QBFVInit, PInput, PDigit, Q, Bigq,
                      scaleTHI, scaleStepTHI, order, numSlots, ringDim, dcrtBits,
                      ep, keyPair, numSlotsCKKS, depth,
                      levelsAvailableBeforeBootstrap, levelsComputation);
}


lbcrypto::Ciphertext<lbcrypto::DCRTPoly> LessEqual(
    const std::vector<lbcrypto::Poly>& ctxtX,
    const std::vector<lbcrypto::Poly>& ctxtY,
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

    // sign(x - y) = 1 if x > y
    // 1 - sign(x - y) = 1 if x <= y
    std::vector<lbcrypto::Poly> diff = {
        ctxtY[0] - ctxtX[0],
        ctxtY[1] - ctxtX[1]
    };

    auto signResult = SignCipher(diff, cc, QBFVInit, PInput, PDigit, Q, Bigq,
                                  scaleTHI, scaleStepTHI, order, numSlots, ringDim, dcrtBits,
                                  ep, keyPair, numSlotsCKKS, depth,
                                  levelsAvailableBeforeBootstrap, levelsComputation);

    
    std::vector<double> x_one = {1};
    if (x_one.size() < numSlots)
        x_one = Fill<double>(x_one, numSlots);

    Plaintext ptxtone = cc->MakeCKKSPackedPlaintext(x_one, 1, depth-6-levelsComputation, nullptr, numSlotsCKKS);
    auto cone = cc->Encrypt(keyPair.publicKey, ptxtone);

    return cc->EvalSub(cone, signResult);
}
