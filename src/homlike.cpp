#include "homlike.h"
#include <chrono>

using namespace lbcrypto;

namespace {

// ---------------------------------------------------------------
// RunOneCase: 调用一次 EqualBootstrapping + 旋转相加
// ---------------------------------------------------------------
Ciphertext<DCRTPoly> RunOneCase(
    const std::vector<Poly>& ctxtStr,
    const std::vector<int64_t>& patternVec,
    uint32_t lmax,
    CryptoContextCKKSRNS::ContextType& cc,
    BigInteger QBFVInit, BigInteger PInput, BigInteger POutput,
    BigInteger Q, BigInteger Bigq, uint64_t scaleTHI,
    size_t order, uint32_t numSlots, uint32_t ringDim, uint32_t dcrtBits,
    std::shared_ptr<M4DCRTParams>& ep, KeyPair<DCRTPoly>& keyPair,
    uint32_t numSlotsCKKS, uint32_t depth,
    uint32_t levelsAvailableBeforeBootstrap, uint32_t levelsComputation) {

    auto asciipattern = patternVec;
    if (asciipattern.size() < numSlots)
        asciipattern = Fill<int64_t>(asciipattern, numSlots);

    auto ctxtPattern = SchemeletRLWEMP::EncryptCoeff(
        asciipattern, QBFVInit, PInput, keyPair.secretKey, ep);

    std::vector<Poly> diff = ctxtStr;
    diff[0] -= ctxtPattern[0];
    diff[1] -= ctxtPattern[1];

    auto eq_check = EqualBootstrapping(diff, cc, QBFVInit, PInput, POutput, Q, Bigq,
                                       scaleTHI, order, numSlots, ringDim, dcrtBits,
                                       ep, keyPair, numSlotsCKKS, depth,
                                       levelsAvailableBeforeBootstrap, levelsComputation);

    Ciphertext<DCRTPoly> rotated_result;
    Ciphertext<DCRTPoly> rotated_sum;

    for (size_t j = 0; j < lmax; ++j) {
        if (j == 0) {
            rotated_sum = eq_check;
            rotated_result = eq_check;
        } else {
            rotated_result = cc->EvalRotate(rotated_result, 1);
            rotated_sum = cc->EvalAdd(rotated_sum, rotated_result);
        }
    }

    std::vector<double> mask(lmax, 0);
    mask[0] = 1;
    if (mask.size() < numSlots)
        mask = Fill<double>(mask, numSlots);

    Plaintext ptxt = cc->MakeCKKSPackedPlaintext(
        mask, 1, rotated_sum->GetLevel(), nullptr, numSlotsCKKS);
    auto c1 = cc->Encrypt(keyPair.publicKey, ptxt);
    auto masked_eq_check = cc->EvalMult(rotated_sum, c1);
    cc->ModReduceInPlace(masked_eq_check);

    return masked_eq_check;
}

// ---------------------------------------------------------------
// Build pattern vector, pad to lmax
// ---------------------------------------------------------------
std::vector<int64_t> BuildPatternVec(const std::string& str, size_t lmax) {
    std::vector<int64_t> vec;
    for (char c : str) {
        if (c == '_') vec.push_back(0);
        else if (c != '%') vec.push_back(static_cast<int64_t>(c));
    }
    while (vec.size() < lmax) vec.push_back(' ');  // 空格填充
    return vec;
}

} // anonymous namespace

// ===============================================================
// Case 1: Direct match (no wildcards) – 1 sub-case
// ===============================================================
Ciphertext<DCRTPoly> RunCase1_Direct(
    const std::vector<Poly>& ctxtStr,
    const std::string& core,
    uint32_t lmax,
    CryptoContextCKKSRNS::ContextType& cc,
    BigInteger QBFVInit, BigInteger PInput, BigInteger POutput,
    BigInteger Q, BigInteger Bigq, uint64_t scaleTHI,
    size_t order, uint32_t numSlots, uint32_t ringDim, uint32_t dcrtBits,
    std::shared_ptr<M4DCRTParams>& ep, KeyPair<DCRTPoly>& keyPair,
    uint32_t numSlotsCKKS, uint32_t depth,
    uint32_t levelsAvailableBeforeBootstrap, uint32_t levelsComputation) {

    std::cerr << "\n===== RunCase1_Direct =====" << std::endl;
    std::cerr << "  Core: \"" << core << "\" expectedScore=" << lmax << std::endl;

    auto case_start = std::chrono::high_resolution_clock::now();

    auto patternVec = BuildPatternVec(core, lmax);
    auto result = RunOneCase(ctxtStr, patternVec, lmax,
                            cc, QBFVInit, PInput, POutput, Q, Bigq, scaleTHI,
                            order, numSlots, ringDim, dcrtBits,
                            ep, keyPair, numSlotsCKKS, depth,
                            levelsAvailableBeforeBootstrap, levelsComputation);

    if (result) {
        auto ctxtAfter = cc->EvalHomDecoding(result, scaleTHI, levelsComputation - 1);
        auto polys = SchemeletRLWEMP::ConvertCKKSToRLWE(ctxtAfter, Q);
        auto computed = SchemeletRLWEMP::DecryptCoeff(
            polys, Q, POutput, keyPair.secretKey, ep, numSlotsCKKS, numSlots);
        std::cerr << "Result: ";
        std::copy_n(computed.begin(), 16, std::ostream_iterator<int64_t>(std::cerr, " "));
        std::cerr << std::endl;
    }
    
    auto case_end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> case_dt = case_end - case_start;
    std::cerr << "  Timing: RunCase1_Direct = " << case_dt.count() << " ms" << std::endl;

    return result;
}

// ===============================================================
// Case 2: _ in front – 1 sub-case
// ===============================================================
Ciphertext<DCRTPoly> RunCase2_UnderFront(
    const std::vector<Poly>& ctxtStr,
    const std::string& core,
    uint32_t lmax,
    CryptoContextCKKSRNS::ContextType& cc,
    BigInteger QBFVInit, BigInteger PInput, BigInteger POutput,
    BigInteger Q, BigInteger Bigq, uint64_t scaleTHI,
    size_t order, uint32_t numSlots, uint32_t ringDim, uint32_t dcrtBits,
    std::shared_ptr<M4DCRTParams>& ep, KeyPair<DCRTPoly>& keyPair,
    uint32_t numSlotsCKKS, uint32_t depth,
    uint32_t levelsAvailableBeforeBootstrap, uint32_t levelsComputation) {

    std::cerr << "\n===== RunCase2_UnderFront =====" << std::endl;
    std::cerr << "  Core: \"" << core << "\" expectedScore=" << lmax << std::endl;

    auto case_start = std::chrono::high_resolution_clock::now();

    std::string raw = "_" + core;
    auto patternVec = BuildPatternVec(raw, lmax);
    auto result = RunOneCase(ctxtStr, patternVec, lmax,
                            cc, QBFVInit, PInput, POutput, Q, Bigq, scaleTHI,
                            order, numSlots, ringDim, dcrtBits,
                            ep, keyPair, numSlotsCKKS, depth,
                            levelsAvailableBeforeBootstrap, levelsComputation);

    if (result) {
        auto ctxtAfter = cc->EvalHomDecoding(result, scaleTHI, levelsComputation - 1);
        auto polys = SchemeletRLWEMP::ConvertCKKSToRLWE(ctxtAfter, Q);
        auto computed = SchemeletRLWEMP::DecryptCoeff(
            polys, Q, POutput, keyPair.secretKey, ep, numSlotsCKKS, numSlots);
        std::cerr << "Result: ";
        std::copy_n(computed.begin(), 16, std::ostream_iterator<int64_t>(std::cerr, " "));
        std::cerr << std::endl;
    }
    auto case_end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> case_dt = case_end - case_start;
    std::cerr << "  Timing: RunCase2_UnderFront = " << case_dt.count() << " ms" << std::endl;

    return result;
}

// ===============================================================
// Case 3: _ in middle – 1 sub-case
// ===============================================================
Ciphertext<DCRTPoly> RunCase3_UnderMid(
    const std::vector<Poly>& ctxtStr,
    const std::string& core,
    uint32_t lmax,
    CryptoContextCKKSRNS::ContextType& cc,
    BigInteger QBFVInit, BigInteger PInput, BigInteger POutput,
    BigInteger Q, BigInteger Bigq, uint64_t scaleTHI,
    size_t order, uint32_t numSlots, uint32_t ringDim, uint32_t dcrtBits,
    std::shared_ptr<M4DCRTParams>& ep, KeyPair<DCRTPoly>& keyPair,
    uint32_t numSlotsCKKS, uint32_t depth,
    uint32_t levelsAvailableBeforeBootstrap, uint32_t levelsComputation) {

    std::cerr << "\n===== RunCase3_UnderMid =====" << std::endl;
    std::cerr << "  Core: \"" << core << "\" expectedScore=" << lmax << std::endl;

    auto case_start = std::chrono::high_resolution_clock::now();

    std::string raw = core;
    raw.insert(raw.begin() + 1, '_');
    auto patternVec = BuildPatternVec(raw, lmax);
    auto result = RunOneCase(ctxtStr, patternVec, lmax,
                            cc, QBFVInit, PInput, POutput, Q, Bigq, scaleTHI,
                            order, numSlots, ringDim, dcrtBits,
                            ep, keyPair, numSlotsCKKS, depth,
                            levelsAvailableBeforeBootstrap, levelsComputation);

    if (result) {
        auto ctxtAfter = cc->EvalHomDecoding(result, scaleTHI, levelsComputation - 1);
        auto polys = SchemeletRLWEMP::ConvertCKKSToRLWE(ctxtAfter, Q);
        auto computed = SchemeletRLWEMP::DecryptCoeff(
            polys, Q, POutput, keyPair.secretKey, ep, numSlotsCKKS, numSlots);
        std::cerr << "Result: ";
        std::copy_n(computed.begin(), 16, std::ostream_iterator<int64_t>(std::cerr, " "));
        std::cerr << std::endl;
    }
    auto case_end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> case_dt = case_end - case_start;
    std::cerr << "  Timing: RunCase3_UnderMid = " << case_dt.count() << " ms" << std::endl;

    return result;
}

// ===============================================================
// Case 4: _ in back – 1 sub-case
// ===============================================================
Ciphertext<DCRTPoly> RunCase4_UnderBack(
    const std::vector<Poly>& ctxtStr,
    const std::string& core,
    uint32_t lmax,
    CryptoContextCKKSRNS::ContextType& cc,
    BigInteger QBFVInit, BigInteger PInput, BigInteger POutput,
    BigInteger Q, BigInteger Bigq, uint64_t scaleTHI,
    size_t order, uint32_t numSlots, uint32_t ringDim, uint32_t dcrtBits,
    std::shared_ptr<M4DCRTParams>& ep, KeyPair<DCRTPoly>& keyPair,
    uint32_t numSlotsCKKS, uint32_t depth,
    uint32_t levelsAvailableBeforeBootstrap, uint32_t levelsComputation) {

    std::cerr << "\n===== RunCase4_UnderBack =====" << std::endl;
    std::cerr << "  Core: \"" << core << "\" expectedScore=" << lmax << std::endl;

    auto case_start = std::chrono::high_resolution_clock::now();

    std::string raw = core + "_";
    auto patternVec = BuildPatternVec(raw, lmax);
    auto result = RunOneCase(ctxtStr, patternVec, lmax,
                            cc, QBFVInit, PInput, POutput, Q, Bigq, scaleTHI,
                            order, numSlots, ringDim, dcrtBits,
                            ep, keyPair, numSlotsCKKS, depth,
                            levelsAvailableBeforeBootstrap, levelsComputation);

    if (result) {
        auto ctxtAfter = cc->EvalHomDecoding(result, scaleTHI, levelsComputation - 1);
        auto polys = SchemeletRLWEMP::ConvertCKKSToRLWE(ctxtAfter, Q);
        auto computed = SchemeletRLWEMP::DecryptCoeff(
            polys, Q, POutput, keyPair.secretKey, ep, numSlotsCKKS, numSlots);
        std::cerr << "Result: ";
        std::copy_n(computed.begin(), 16, std::ostream_iterator<int64_t>(std::cerr, " "));
        std::cerr << std::endl;
    }
    auto case_end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> case_dt = case_end - case_start;
    std::cerr << "  Timing: RunCase4_UnderBack = " << case_dt.count() << " ms" << std::endl;

    return result;
}

// ===============================================================
// Case 5: % in front – multiple sub-cases, results ADDED
// ===============================================================
Ciphertext<DCRTPoly> RunCase5_PercentFront(
    const std::vector<Poly>& ctxtStr,
    const std::string& core,
    uint32_t lmax,
    CryptoContextCKKSRNS::ContextType& cc,
    BigInteger QBFVInit, BigInteger PInput, BigInteger POutput,
    BigInteger Q, BigInteger Bigq, uint64_t scaleTHI,
    size_t order, uint32_t numSlots, uint32_t ringDim, uint32_t dcrtBits,
    std::shared_ptr<M4DCRTParams>& ep, KeyPair<DCRTPoly>& keyPair,
    uint32_t numSlotsCKKS, uint32_t depth,
    uint32_t levelsAvailableBeforeBootstrap, uint32_t levelsComputation) {

    std::cerr << "\n===== RunCase5_PercentFront =====" << std::endl;
    std::cerr << "  Core: \"" << core << "\" expectedScore=" << lmax << std::endl;

    auto case_start = std::chrono::high_resolution_clock::now();
    
    size_t coreLen = core.size();
    size_t numUnderscores = lmax - coreLen;
    std::cerr << "  Core: \"" << core << "\" numUnderscores=" << numUnderscores << std::endl;

    Ciphertext<DCRTPoly> result = nullptr;

    for (size_t u = 0; u < numUnderscores; ++u) {
        std::string raw;
        for (size_t i = 0; i < u; ++i) raw += '_';
        raw += core;
        auto patternVec = BuildPatternVec(raw, lmax);
        std::cerr << "  u=" << u << " raw=\"" << raw << "\" expectedScore=" << (lmax - u) << std::endl;

        auto subResult = RunOneCase(ctxtStr, patternVec, lmax,
                                   cc, QBFVInit, PInput, POutput, Q, Bigq, scaleTHI,
                                   order, numSlots, ringDim, dcrtBits,
                                   ep, keyPair, numSlotsCKKS, depth,
                                   levelsAvailableBeforeBootstrap, levelsComputation);
        if (result == nullptr) result = subResult;
        else result = cc->EvalAdd(result, subResult);
    }

    if (result) {
        auto ctxtAfter = cc->EvalHomDecoding(result, scaleTHI, levelsComputation - 1);
        auto polys = SchemeletRLWEMP::ConvertCKKSToRLWE(ctxtAfter, Q);
        auto computed = SchemeletRLWEMP::DecryptCoeff(
            polys, Q, POutput, keyPair.secretKey, ep, numSlotsCKKS, numSlots);
        std::cerr << "Result: ";
        std::copy_n(computed.begin(), 16, std::ostream_iterator<int64_t>(std::cerr, " "));
        std::cerr << std::endl;
    }
    auto case_end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> case_dt = case_end - case_start;
    std::cerr << "  Timing: RunCase5_PercentFront = " << case_dt.count() << " ms" << std::endl;
    return result;
}

// ===============================================================
// Case 6: % in middle – multiple sub-cases, results ADDED
// ===============================================================
Ciphertext<DCRTPoly> RunCase6_PercentMid(
    const std::vector<Poly>& ctxtStr,
    const std::string& core,
    uint32_t lmax,
    CryptoContextCKKSRNS::ContextType& cc,
    BigInteger QBFVInit, BigInteger PInput, BigInteger POutput,
    BigInteger Q, BigInteger Bigq, uint64_t scaleTHI,
    size_t order, uint32_t numSlots, uint32_t ringDim, uint32_t dcrtBits,
    std::shared_ptr<M4DCRTParams>& ep, KeyPair<DCRTPoly>& keyPair,
    uint32_t numSlotsCKKS, uint32_t depth,
    uint32_t levelsAvailableBeforeBootstrap, uint32_t levelsComputation) {

    std::cerr << "\n===== RunCase6_PercentMid =====" << std::endl;
    std::cerr << "  Core: \"" << core << "\" expectedScore=" << lmax << std::endl;

    auto case_start = std::chrono::high_resolution_clock::now();

    size_t coreLen = core.size();
    size_t numUnderscores = lmax - coreLen;
    std::cerr << "  Core: \"" << core << "\" numUnderscores=" << numUnderscores << std::endl;

    Ciphertext<DCRTPoly> result = nullptr;

    for (size_t u = 1; u <= numUnderscores; ++u) {
        std::string raw;
        raw += core[0];
        for (size_t i = 0; i < u; ++i) raw += '_';
        raw += core.substr(1);
        auto patternVec = BuildPatternVec(raw, lmax);
        std::cerr << "  u=" << u << " raw=\"" << raw << "\" expectedScore=" << (lmax - u + 1) << std::endl;

        auto subResult = RunOneCase(ctxtStr, patternVec, lmax,
                                   cc, QBFVInit, PInput, POutput, Q, Bigq, scaleTHI,
                                   order, numSlots, ringDim, dcrtBits,
                                   ep, keyPair, numSlotsCKKS, depth,
                                   levelsAvailableBeforeBootstrap, levelsComputation);
        if (result == nullptr) result = subResult;
        else result = cc->EvalAdd(result, subResult);
    }

    if (result) {
        auto ctxtAfter = cc->EvalHomDecoding(result, scaleTHI, levelsComputation - 1);
        auto polys = SchemeletRLWEMP::ConvertCKKSToRLWE(ctxtAfter, Q);
        auto computed = SchemeletRLWEMP::DecryptCoeff(
            polys, Q, POutput, keyPair.secretKey, ep, numSlotsCKKS, numSlots);
        std::cerr << "Result: ";
        std::copy_n(computed.begin(), 16, std::ostream_iterator<int64_t>(std::cerr, " "));
        std::cerr << std::endl;
    }
    auto case_end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> case_dt = case_end - case_start;
    std::cerr << "  Timing: RunCase6_PercentMid = " << case_dt.count() << " ms" << std::endl;
    return result;
}

// ===============================================================
// Case 7: % in back – multiple sub-cases, results ADDED
// ===============================================================
Ciphertext<DCRTPoly> RunCase7_PercentBack(
    const std::vector<Poly>& ctxtStr,
    const std::string& core,
    uint32_t lmax,
    CryptoContextCKKSRNS::ContextType& cc,
    BigInteger QBFVInit, BigInteger PInput, BigInteger POutput,
    BigInteger Q, BigInteger Bigq, uint64_t scaleTHI,
    size_t order, uint32_t numSlots, uint32_t ringDim, uint32_t dcrtBits,
    std::shared_ptr<M4DCRTParams>& ep, KeyPair<DCRTPoly>& keyPair,
    uint32_t numSlotsCKKS, uint32_t depth,
    uint32_t levelsAvailableBeforeBootstrap, uint32_t levelsComputation) {

    std::cerr << "\n===== RunCase7_PercentBack =====" << std::endl;
    std::cerr << "  Core: \"" << core << "\" expectedScore=" << lmax << std::endl;

    auto case_start = std::chrono::high_resolution_clock::now();

    size_t coreLen = core.size();
    size_t numUnderscores = lmax - coreLen;
    std::cerr << "  Core: \"" << core << "\" numUnderscores=" << numUnderscores << std::endl;

    Ciphertext<DCRTPoly> result = nullptr;

    for (size_t u = 1; u <= numUnderscores; ++u) {
        std::string raw = core;
        for (size_t i = 0; i < u; ++i) raw += '_';
        auto patternVec = BuildPatternVec(raw, lmax);
        std::cerr << "  u=" << u << " raw=\"" << raw << "\" expectedScore=" << (lmax - u + 1) << std::endl;

        auto subResult = RunOneCase(ctxtStr, patternVec, lmax,
                                   cc, QBFVInit, PInput, POutput, Q, Bigq, scaleTHI,
                                   order, numSlots, ringDim, dcrtBits,
                                   ep, keyPair, numSlotsCKKS, depth,
                                   levelsAvailableBeforeBootstrap, levelsComputation);
        if (result == nullptr) result = subResult;
        else result = cc->EvalAdd(result, subResult);
    }

    if (result) {
        auto ctxtAfter = cc->EvalHomDecoding(result, scaleTHI, levelsComputation - 1);
        auto polys = SchemeletRLWEMP::ConvertCKKSToRLWE(ctxtAfter, Q);
        auto computed = SchemeletRLWEMP::DecryptCoeff(
            polys, Q, POutput, keyPair.secretKey, ep, numSlotsCKKS, numSlots);
        std::cerr << "Result: ";
        std::copy_n(computed.begin(), 16, std::ostream_iterator<int64_t>(std::cerr, " "));
        std::cerr << std::endl;
    }
    auto case_end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> case_dt = case_end - case_start;
    std::cerr << "  Timing: RunCase7_PercentBack = " << case_dt.count() << " ms" << std::endl;
    return result;
}
