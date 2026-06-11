#define PROFILE
#include "range.h"
#include <chrono>

using namespace lbcrypto;

void BetweenSignCore(
    std::vector<Poly>& ctxtBFV1,
    std::vector<Poly>& ctxtBFV2,
    CryptoContextCKKSRNS::ContextType& cc,
    BigInteger QBFVInit,
    BigInteger PInput,
    BigInteger PDigit,
    BigInteger& Qorig,
    BigInteger& pOrig,
    BigInteger Q,
    BigInteger Bigq,
    uint64_t scaleTHI,
    uint64_t scaleStepTHI,
    size_t order,
    uint32_t numSlots,
    uint32_t ringDim,
    uint32_t dcrtBits,
    std::shared_ptr<M4DCRTParams>& ep,
    KeyPair<DCRTPoly>& keyPair,
    uint32_t numSlotsCKKS,
    uint32_t depth,
    uint32_t levelsAvailableBeforeBootstrap,
    uint32_t levelsComputation,
    Ciphertext<DCRTPoly>& ctxtAfterFBT1,
    Ciphertext<DCRTPoly>& ctxtAfterFBT2) {

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

    SchemeletRLWEMP::ModSwitch(ctxtBFV1, Q, QBFVInit);
    SchemeletRLWEMP::ModSwitch(ctxtBFV2, Q, QBFVInit);
    uint32_t QBFVBits = Q.GetMSB() - 1;
    std::cout << "QBFVBits:" << QBFVBits << "\n";

    std::vector<int64_t> coeffint;
    std::vector<std::complex<double>> coeffcomp;
    if (binaryLUT)
        coeffint = coeffintMod;
    else
        coeffcomp = coeffcompMod;

    const bool checkeq2       = PDigit.ConvertToInt() == 2;
    const bool checkgt2       = PDigit.ConvertToInt() > 2;
    const uint32_t pDigitBits = PDigit.GetMSB() - 1;

    Qorig = Q;
    pOrig = PInput;

    BigInteger QNew;

    bool step                = false;
    bool go                  = QBFVBits > dcrtBits;
    uint32_t postScalingBits = 0;

    while (go) {
        auto encryptedDigit1 = ctxtBFV1;
        auto encryptedDigit2 = ctxtBFV2;

        encryptedDigit1[0].SwitchModulus(Bigq, 1, 0, 0);
        encryptedDigit1[1].SwitchModulus(Bigq, 1, 0, 0);

        encryptedDigit2[0].SwitchModulus(Bigq, 1, 0, 0);
        encryptedDigit2[1].SwitchModulus(Bigq, 1, 0, 0);

        auto ctxt1 = SchemeletRLWEMP::ConvertRLWEToCKKS(*cc, encryptedDigit1, keyPair.publicKey, Bigq, numSlotsCKKS,
                                                       depth - (levelsAvailableBeforeBootstrap > 0));
        auto ctxt2 = SchemeletRLWEMP::ConvertRLWEToCKKS(*cc, encryptedDigit2, keyPair.publicKey, Bigq, numSlotsCKKS,
                                                       depth - (levelsAvailableBeforeBootstrap > 0));

        bool isLastIteration = step || (checkeq2 && (QBFVBits - pDigitBits <= dcrtBits));
        if (isLastIteration) {
            if (binaryLUT) {
                ctxtAfterFBT1 = cc->EvalFBTNoDecoding(ctxt1, coeffint, pDigitBits, ep->GetModulus(), order);
                ctxtAfterFBT2 = cc->EvalFBTNoDecoding(ctxt2, coeffint, pDigitBits, ep->GetModulus(), order);
            } else {
                ctxtAfterFBT1 = cc->EvalFBTNoDecoding(ctxt1, coeffcomp, pDigitBits, ep->GetModulus(), order);
                ctxtAfterFBT2 = cc->EvalFBTNoDecoding(ctxt2, coeffcomp, pDigitBits, ep->GetModulus(), order);
            }
            break;
        }
        else {
            if (binaryLUT) {
                ctxtAfterFBT1 = cc->EvalFBT(ctxt1, coeffint, pDigitBits, ep->GetModulus(), scaleTHI * (1 << postScalingBits),
                                           levelsComputation, order);
                ctxtAfterFBT2 = cc->EvalFBT(ctxt2, coeffint, pDigitBits, ep->GetModulus(), scaleTHI * (1 << postScalingBits),
                                           levelsComputation, order);
            } else {
                ctxtAfterFBT1 = cc->EvalFBT(ctxt1, coeffcomp, pDigitBits, ep->GetModulus(), scaleTHI * (1 << postScalingBits),
                                           levelsComputation, order);
                ctxtAfterFBT2 = cc->EvalFBT(ctxt2, coeffcomp, pDigitBits, ep->GetModulus(), scaleTHI * (1 << postScalingBits),
                                           levelsComputation, order);
            }

            auto polys1 = SchemeletRLWEMP::ConvertCKKSToRLWE(ctxtAfterFBT1, Q);
            auto polys2 = SchemeletRLWEMP::ConvertCKKSToRLWE(ctxtAfterFBT2, Q);
            ctxtBFV1[0] = ctxtBFV1[0] - polys1[0];
            ctxtBFV1[1] = ctxtBFV1[1] - polys1[1];

            ctxtBFV2[0] = ctxtBFV2[0] - polys2[0];
            ctxtBFV2[1] = ctxtBFV2[1] - polys2[1];

            QNew       = Q >> pDigitBits;
            ctxtBFV1[0] = ctxtBFV1[0].MultiplyAndRound(QNew, Q);
            ctxtBFV1[0].SwitchModulus(QNew, 1, 0, 0);
            ctxtBFV1[1] = ctxtBFV1[1].MultiplyAndRound(QNew, Q);
            ctxtBFV1[1].SwitchModulus(QNew, 1, 0, 0);

            ctxtBFV2[0] = ctxtBFV2[0].MultiplyAndRound(QNew, Q);
            ctxtBFV2[0].SwitchModulus(QNew, 1, 0, 0);
            ctxtBFV2[1] = ctxtBFV2[1].MultiplyAndRound(QNew, Q);
            ctxtBFV2[1].SwitchModulus(QNew, 1, 0, 0);

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

}


// ===== a < x < b (strict on both sides) =====
// diff1 = b - x, sign(b-x) = 1 when x < b
// diff2 = x - a, sign(x-a) = 1 when x > a
// result = 1 - (sign1 + sign2)
Ciphertext<DCRTPoly> HomRangeStrict(
    const std::vector<Poly>& ctxtX,
    const std::vector<Poly>& ctxtA,
    const std::vector<Poly>& ctxtB,
    CryptoContextCKKSRNS::ContextType& cc,
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
    std::shared_ptr<M4DCRTParams>& ep,
    KeyPair<DCRTPoly>& keyPair,
    uint32_t numSlotsCKKS,
    uint32_t depth,
    uint32_t levelsAvailableBeforeBootstrap,
    uint32_t levelsComputation) {

    std::cerr << "\n=== HomRangeStrict (a < x < b) ===" << std::endl;

    auto start = std::chrono::high_resolution_clock::now();

    std::vector<Poly> diff1 = { ctxtX[0] - ctxtB[0], ctxtX[1] - ctxtB[1] };
    std::vector<Poly> diff2 = { ctxtA[0] - ctxtX[0], ctxtA[1] - ctxtX[1] };

    BigInteger Qorig, pOrig;
    Ciphertext<DCRTPoly> ctxtAfterFBT1, ctxtAfterFBT2;
    BetweenSignCore(diff1, diff2, cc, QBFVInit, PInput, PDigit, Qorig, pOrig, Q, Bigq,
                    scaleTHI, scaleStepTHI, order, numSlots, ringDim, dcrtBits,
                    ep, keyPair, numSlotsCKKS, depth,
                    levelsAvailableBeforeBootstrap, levelsComputation,
                    ctxtAfterFBT1, ctxtAfterFBT2);

    // result = 1 - (sign1 + sign2)
    std::vector<double> x_input11 = {1};
    if (x_input11.size() < numSlots)
        x_input11 = Fill<double>(x_input11, numSlots);

    Plaintext ptxt1 = cc->MakeCKKSPackedPlaintext(x_input11, 1, depth-6-levelsComputation, nullptr, numSlotsCKKS);

    auto ctxtAfterFBT5 = cc->EvalAdd(cc->EvalSub(ptxt1, ctxtAfterFBT1), cc->EvalSub(ptxt1, ctxtAfterFBT2));
    auto ctxtAfterFBT7 = cc->EvalSub(ptxt1, ctxtAfterFBT5);

    auto ctxtAfterFBT6 = cc->EvalHomDecoding(ctxtAfterFBT7, scaleTHI, levelsComputation);
    auto ctxtResult = SchemeletRLWEMP::ConvertCKKSToRLWE(ctxtAfterFBT6, Qorig);
    auto computed =
        SchemeletRLWEMP::DecryptCoeff(ctxtResult, Qorig, pOrig, keyPair.secretKey, ep, numSlotsCKKS, numSlots);

    std::cerr << "HomRangeStrict result: [";
    std::copy_n(computed.begin(), numSlots, std::ostream_iterator<int64_t>(std::cerr, " "));
    std::cerr << "]" << std::endl;

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> elapsed = end - start;
    std::cout << "HomRangeStrict time: " << elapsed.count() << " ms" << std::endl;

    return {};
}


// ===== a <= x < b (lower inclusive) =====
// diff1 = b - x, sign(b-x) = 1 when x < b
// diff2 = a - x, sign(a-x) = 1 when x < a
// result = 1 - (sign1 + sign2)
Ciphertext<DCRTPoly> HomRangeLowerInc(
    const std::vector<Poly>& ctxtX,
    const std::vector<Poly>& ctxtA,
    const std::vector<Poly>& ctxtB,
    CryptoContextCKKSRNS::ContextType& cc,
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
    std::shared_ptr<M4DCRTParams>& ep,
    KeyPair<DCRTPoly>& keyPair,
    uint32_t numSlotsCKKS,
    uint32_t depth,
    uint32_t levelsAvailableBeforeBootstrap,
    uint32_t levelsComputation) {

    std::cerr << "\n=== HomRangeLowerInc (a <= x < b) ===" << std::endl;

    auto start = std::chrono::high_resolution_clock::now();

    std::vector<Poly> diff1 = { ctxtX[0] - ctxtB[0], ctxtX[1] - ctxtB[1] };
    std::vector<Poly> diff2 = { ctxtX[0] - ctxtA[0], ctxtX[1] - ctxtA[1] };

    BigInteger Qorig, pOrig;
    Ciphertext<DCRTPoly> ctxtAfterFBT1, ctxtAfterFBT2;
    BetweenSignCore(diff1, diff2, cc, QBFVInit, PInput, PDigit, Qorig, pOrig, Q, Bigq,
                    scaleTHI, scaleStepTHI, order, numSlots, ringDim, dcrtBits,
                    ep, keyPair, numSlotsCKKS, depth,
                    levelsAvailableBeforeBootstrap, levelsComputation,
                    ctxtAfterFBT1, ctxtAfterFBT2);

    
    std::vector<double> x_input11 = {1};
    if (x_input11.size() < numSlots)
        x_input11 = Fill<double>(x_input11, numSlots);

    Plaintext ptxt1 = cc->MakeCKKSPackedPlaintext(x_input11, 1, depth-6-levelsComputation, nullptr, numSlotsCKKS);

    auto ctxtAfterFBT5 = cc->EvalAdd(cc->EvalSub(ptxt1, ctxtAfterFBT1), ctxtAfterFBT2);
    auto ctxtAfterFBT7 = cc->EvalSub(ptxt1, ctxtAfterFBT5);

    auto ctxtAfterFBT6 = cc->EvalHomDecoding(ctxtAfterFBT7, scaleTHI, levelsComputation);
    auto ctxtResult = SchemeletRLWEMP::ConvertCKKSToRLWE(ctxtAfterFBT6, Qorig);
    auto computed =
        SchemeletRLWEMP::DecryptCoeff(ctxtResult, Qorig, pOrig, keyPair.secretKey, ep, numSlotsCKKS, numSlots);

    std::cerr << "HomRangeLowerInc result: [";
    std::copy_n(computed.begin(), numSlots, std::ostream_iterator<int64_t>(std::cerr, " "));
    std::cerr << "]" << std::endl;

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> elapsed = end - start;
    std::cout << "HomRangeLowerInc time: " << elapsed.count() << " ms" << std::endl;

    return {};
}


// ===== a < x <= b (upper inclusive) =====
// diff1 = x - b, sign(x-b) = 1 when x > b
// diff2 = x - a, sign(x-a) = 1 when x > a
// result = 1 - (sign1 + sign2)
Ciphertext<DCRTPoly> HomRangeUpperInc(
    const std::vector<Poly>& ctxtX,
    const std::vector<Poly>& ctxtA,
    const std::vector<Poly>& ctxtB,
    CryptoContextCKKSRNS::ContextType& cc,
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
    std::shared_ptr<M4DCRTParams>& ep,
    KeyPair<DCRTPoly>& keyPair,
    uint32_t numSlotsCKKS,
    uint32_t depth,
    uint32_t levelsAvailableBeforeBootstrap,
    uint32_t levelsComputation) {

    std::cerr << "\n=== HomRangeUpperInc (a < x <= b) ===" << std::endl;

    auto start = std::chrono::high_resolution_clock::now();

    std::vector<Poly> diff1 = { ctxtB[0] - ctxtX[0], ctxtB[1] - ctxtX[1] };
    std::vector<Poly> diff2 = { ctxtA[0] - ctxtX[0], ctxtA[1] - ctxtX[1] };

    BigInteger Qorig, pOrig;
    Ciphertext<DCRTPoly> ctxtAfterFBT1, ctxtAfterFBT2;
    BetweenSignCore(diff1, diff2, cc, QBFVInit, PInput, PDigit, Qorig, pOrig, Q, Bigq,
                    scaleTHI, scaleStepTHI, order, numSlots, ringDim, dcrtBits,
                    ep, keyPair, numSlotsCKKS, depth,
                    levelsAvailableBeforeBootstrap, levelsComputation,
                    ctxtAfterFBT1, ctxtAfterFBT2);

    // result = 1 - (sign1 + sign2)
    std::vector<double> x_input11 = {1};
    if (x_input11.size() < numSlots)
        x_input11 = Fill<double>(x_input11, numSlots);

    Plaintext ptxt1 = cc->MakeCKKSPackedPlaintext(x_input11, 1, depth-6-levelsComputation, nullptr, numSlotsCKKS);

    auto ctxtAfterFBT5 = cc->EvalAdd(ctxtAfterFBT1, cc->EvalSub(ptxt1, ctxtAfterFBT2));
    auto ctxtAfterFBT7 = cc->EvalSub(ptxt1, ctxtAfterFBT5);

    auto ctxtAfterFBT6 = cc->EvalHomDecoding(ctxtAfterFBT7, scaleTHI, levelsComputation);
    auto ctxtResult = SchemeletRLWEMP::ConvertCKKSToRLWE(ctxtAfterFBT6, Qorig);
    auto computed =
        SchemeletRLWEMP::DecryptCoeff(ctxtResult, Qorig, pOrig, keyPair.secretKey, ep, numSlotsCKKS, numSlots);

    std::cerr << "HomRangeUpperInc result: [";
    std::copy_n(computed.begin(), numSlots, std::ostream_iterator<int64_t>(std::cerr, " "));
    std::cerr << "]" << std::endl;

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> elapsed = end - start;
    std::cout << "HomRangeUpperInc time: " << elapsed.count() << " ms" << std::endl;

    return {};
}


// ===== a <= x <= b (both inclusive) =====
// diff1 = a - x, sign(a-x) = 1 when x < a
// diff2 = x - b, sign(x-b) = 1 when x > b
// result = 1 - (sign1 + sign2)
Ciphertext<DCRTPoly> HomRangeBothInc(
    const std::vector<Poly>& ctxtX,
    const std::vector<Poly>& ctxtA,
    const std::vector<Poly>& ctxtB,
    CryptoContextCKKSRNS::ContextType& cc,
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
    std::shared_ptr<M4DCRTParams>& ep,
    KeyPair<DCRTPoly>& keyPair,
    uint32_t numSlotsCKKS,
    uint32_t depth,
    uint32_t levelsAvailableBeforeBootstrap,
    uint32_t levelsComputation) {

    std::cerr << "\n=== HomRangeBothInc (a <= x <= b) ===" << std::endl;

    auto start = std::chrono::high_resolution_clock::now();

    std::vector<Poly> diff1 = { ctxtX[0] - ctxtA[0], ctxtX[1] - ctxtA[1] };
    std::vector<Poly> diff2 = { ctxtB[0] - ctxtX[0], ctxtB[1] - ctxtX[1] };

    BigInteger Qorig, pOrig;
    Ciphertext<DCRTPoly> ctxtAfterFBT1, ctxtAfterFBT2;
    BetweenSignCore(diff1, diff2, cc, QBFVInit, PInput, PDigit, Qorig, pOrig, Q, Bigq,
                    scaleTHI, scaleStepTHI, order, numSlots, ringDim, dcrtBits,
                    ep, keyPair, numSlotsCKKS, depth,
                    levelsAvailableBeforeBootstrap, levelsComputation,
                    ctxtAfterFBT1, ctxtAfterFBT2);

    // result = 1 - (sign1 + sign2)
    std::vector<double> x_input11 = {1};
    if (x_input11.size() < numSlots)
        x_input11 = Fill<double>(x_input11, numSlots);

    Plaintext ptxt1 = cc->MakeCKKSPackedPlaintext(x_input11, 1, depth-6-levelsComputation, nullptr, numSlotsCKKS);

    auto ctxtAfterFBT5 = cc->EvalAdd(ctxtAfterFBT1, ctxtAfterFBT2);
    auto ctxtAfterFBT7 = cc->EvalSub(ptxt1, ctxtAfterFBT5);

    auto ctxtAfterFBT6 = cc->EvalHomDecoding(ctxtAfterFBT7, scaleTHI, levelsComputation);
    auto ctxtResult = SchemeletRLWEMP::ConvertCKKSToRLWE(ctxtAfterFBT6, Qorig);
    auto computed =
        SchemeletRLWEMP::DecryptCoeff(ctxtResult, Qorig, pOrig, keyPair.secretKey, ep, numSlotsCKKS, numSlots);

    std::cerr << "HomRangeBothInc result: [";
    std::copy_n(computed.begin(), numSlots, std::ostream_iterator<int64_t>(std::cerr, " "));
    std::cerr << "]" << std::endl;

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> elapsed = end - start;
    std::cout << "HomRangeBothInc time: " << elapsed.count() << " ms" << std::endl;

    return {};
}
