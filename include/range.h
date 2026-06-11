#ifndef RANGE_H
#define RANGE_H

#include "binfhecontext.h"
#include "math/hermite.h"
#include "openfhe.h"
#ifndef SCHEMELET_RLWE_MP_GUARD
#define SCHEMELET_RLWE_MP_GUARD
#include "schemelet/rlwe-mp.h"
#endif

#include <functional>
#include <chrono>
#include <cmath>

using namespace lbcrypto;

/**
 * @brief BetweenSignCore - core sign computation for two RLWE ciphertexts
 *
 * Computes sign(diff1) and sign(diff2) simultaneously using the digit-extraction
 * + FBT pipeline (same as SignCipher but for two inputs at once).
 *
 * Logic (copied from BewteenSign in testbetween.cpp):
 *   sign1 = sign(ctxtBFV1) = 1 if ctxtBFV1 > 0, else 0
 *   sign2 = sign(ctxtBFV2) = 1 if ctxtBFV2 > 0, else 0
 *
 * Returns {sign1, sign2} in CKKS domain (EvalFBTNoDecoding output) via output params.
 * The caller is responsible for computing result = 1 - (sign1 + sign2), decoding,
 * and timing.
 *
 * NOTE: ctxtBFV1 and ctxtBFV2 are MODIFIED (ModSwitch + digit extraction).
 * Caller should pass copies if originals need to be preserved.
 */
void BetweenSignCore(
    std::vector<lbcrypto::Poly>& ctxtBFV1,
    std::vector<lbcrypto::Poly>& ctxtBFV2,
    lbcrypto::CryptoContextCKKSRNS::ContextType& cc,
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
    std::shared_ptr<lbcrypto::M4DCRTParams>& ep,
    lbcrypto::KeyPair<lbcrypto::DCRTPoly>& keyPair,
    uint32_t numSlotsCKKS,
    uint32_t depth,
    uint32_t levelsAvailableBeforeBootstrap,
    uint32_t levelsComputation,
    lbcrypto::Ciphertext<lbcrypto::DCRTPoly>& ctxtAfterFBT1,
    lbcrypto::Ciphertext<lbcrypto::DCRTPoly>& ctxtAfterFBT2);

/**
 * @brief HomRangeStrict - checks if a < x < b (strict on both sides)
 *
 * diff1 = b - x, diff2 = x - a
 * result = sign(b-x) + sign(x-a) - 1
 */
lbcrypto::Ciphertext<lbcrypto::DCRTPoly> HomRangeStrict(
    const std::vector<lbcrypto::Poly>& ctxtX,
    const std::vector<lbcrypto::Poly>& ctxtA,
    const std::vector<lbcrypto::Poly>& ctxtB,
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
    uint32_t levelsComputation);

/**
 * @brief HomRangeLowerInc - checks if a <= x < b (lower inclusive)
 *
 * diff1 = b - x, diff2 = a - x
 * sign1 = sign(b-x), sign2 = 1 - sign(a-x)
 * result = sign1 + sign2 - 1
 */
lbcrypto::Ciphertext<lbcrypto::DCRTPoly> HomRangeLowerInc(
    const std::vector<lbcrypto::Poly>& ctxtX,
    const std::vector<lbcrypto::Poly>& ctxtA,
    const std::vector<lbcrypto::Poly>& ctxtB,
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
    uint32_t levelsComputation);

/**
 * @brief HomRangeUpperInc - checks if a < x <= b (upper inclusive)
 *
 * diff1 = x - b, diff2 = x - a
 * sign1 = 1 - sign(x-b), sign2 = sign(x-a)
 * result = sign1 + sign2 - 1
 */
lbcrypto::Ciphertext<lbcrypto::DCRTPoly> HomRangeUpperInc(
    const std::vector<lbcrypto::Poly>& ctxtX,
    const std::vector<lbcrypto::Poly>& ctxtA,
    const std::vector<lbcrypto::Poly>& ctxtB,
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
    uint32_t levelsComputation);

/**
 * @brief HomRangeBothInc - checks if a <= x <= b (both inclusive)
 *
 * diff1 = a - x, diff2 = x - b
 * sign1 = 1 - sign(a-x), sign2 = 1 - sign(x-b)
 * result = sign1 + sign2 - 1
 */
lbcrypto::Ciphertext<lbcrypto::DCRTPoly> HomRangeBothInc(
    const std::vector<lbcrypto::Poly>& ctxtX,
    const std::vector<lbcrypto::Poly>& ctxtA,
    const std::vector<lbcrypto::Poly>& ctxtB,
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
    uint32_t levelsComputation);

#endif // RANGE_H
