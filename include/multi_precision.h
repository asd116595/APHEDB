#ifndef MULTI_PRECISION_H
#define MULTI_PRECISION_H

#include "binfhecontext.h"
#include "math/hermite.h"
#include "openfhe.h"
#ifndef SCHEMELET_RLWE_MP_GUARD
#define SCHEMELET_RLWE_MP_GUARD
#include "schemelet/rlwe-mp.h"
#endif

#include "comparsion.h"

#include <functional>
#include <chrono>
#include <cmath>
#include <iterator>

using namespace lbcrypto;

/**
 * @brief MultiPrecisionSign - compares two multi-precision ciphertext arrays
 *
 * Each number is represented as an array of RLWE ciphertexts, where each
 * ciphertext encrypts one "digit" (e.g., 8 bits). The first element is the
 * most significant digit.
 *
 * Algorithm:
 *   For i from 0 to numDigits-1 (MSB to LSB):
 *     - Compute greater_i = sign(X[i] - Y[i])  (1 if X[i] > Y[i])
 *     - Compute equal_i   = (X[i] == Y[i])     (1 if equal)
 *     - Compute result[i] = greater_i + equal_i * result[i+1]
 *   Final: result[0] = 1 if X > Y, else 0
 *
 * This is done recursively from LSB to MSB.
 *
 * @param ctxtX     Array of RLWE ciphertexts for X (MSB first)
 * @param ctxtY     Array of RLWE ciphertexts for Y (MSB first)
 * @param numDigits Number of digits in each array
 * @param cc        CKKS crypto context
 * @param QBFVInit  Initial modulus for BFV
 * @param PInput    Input plaintext modulus
 * @param PDigit    Digit modulus for sign evaluation
 * @param Q         RLWE modulus (original copy, will be modified)
 * @param Bigq      CKKS scaling factor
 * @param scaleTHI  Scaling factor for homomorphic decoding
 * @param scaleStepTHI Scaling factor for step function
 * @param order     Order of Hermite interpolation
 * @param numSlots  Number of slots
 * @param ringDim   Ring dimension
 * @param dcrtBits  DCRT bits
 * @param ep        Element parameters
 * @param keyPair   Key pair
 * @param numSlotsCKKS Number of CKKS slots
 * @param depth     Multiplicative depth
 * @param levelsAvailableBeforeBootstrap Levels before bootstrap
 * @param levelsComputation Levels for computation
 * @return Ciphertext<DCRTPoly> Encrypted 1 if X > Y, else encrypted 0
 */
lbcrypto::Ciphertext<lbcrypto::DCRTPoly> MultiPrecisionSign(
    const std::vector<std::vector<lbcrypto::Poly>>& ctxtX,
    const std::vector<std::vector<lbcrypto::Poly>>& ctxtY,
    uint32_t numDigits,
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
 * @brief MultiPrecisionEqual - checks if two multi-precision ciphertext arrays are equal
 *
 * @return Ciphertext<DCRTPoly> Encrypted 1 if X == Y, else encrypted 0
 */
lbcrypto::Ciphertext<lbcrypto::DCRTPoly> MultiPrecisionEqual(
    const std::vector<std::vector<lbcrypto::Poly>>& ctxtX,
    const std::vector<std::vector<lbcrypto::Poly>>& ctxtY,
    uint32_t numDigits,
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
 * @brief MultiPrecisionGreaterThan - checks if X > Y (multi-precision)
 */
lbcrypto::Ciphertext<lbcrypto::DCRTPoly> MultiPrecisionGreaterThan(
    const std::vector<std::vector<lbcrypto::Poly>>& ctxtX,
    const std::vector<std::vector<lbcrypto::Poly>>& ctxtY,
    uint32_t numDigits,
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
 * @brief MultiPrecisionGreaterEqual - checks if X >= Y (multi-precision)
 */
lbcrypto::Ciphertext<lbcrypto::DCRTPoly> MultiPrecisionGreaterEqual(
    const std::vector<std::vector<lbcrypto::Poly>>& ctxtX,
    const std::vector<std::vector<lbcrypto::Poly>>& ctxtY,
    uint32_t numDigits,
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
 * @brief MultiPrecisionLessThan - checks if X < Y (multi-precision)
 */
lbcrypto::Ciphertext<lbcrypto::DCRTPoly> MultiPrecisionLessThan(
    const std::vector<std::vector<lbcrypto::Poly>>& ctxtX,
    const std::vector<std::vector<lbcrypto::Poly>>& ctxtY,
    uint32_t numDigits,
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
 * @brief MultiPrecisionLessEqual - checks if X <= Y (multi-precision)
 */
lbcrypto::Ciphertext<lbcrypto::DCRTPoly> MultiPrecisionLessEqual(
    const std::vector<std::vector<lbcrypto::Poly>>& ctxtX,
    const std::vector<std::vector<lbcrypto::Poly>>& ctxtY,
    uint32_t numDigits,
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
 * @brief EqualBootstrapping - checks if a ciphertext value equals zero
 *
 * Uses FBT with equal-LUT: func(x) = 1 if x == 0, else 0
 * Returns the result as a CKKS ciphertext (EvalFBTNoDecoding output).
 *
 * @param ctxtBFV   Input RLWE ciphertext
 * @param cc        CKKS crypto context
 * @param QBFVInit  Initial modulus
 * @param PInput    Input plaintext modulus
 * @param POutput   Output plaintext modulus (typically same as PInput)
 * @param Q         RLWE modulus
 * @param Bigq      CKKS scaling factor
 * @param scaleTHI  Scaling factor
 * @param order     Order of Hermite interpolation
 * @param numSlots  Number of slots
 * @param ringDim   Ring dimension
 * @param dcrtBits  DCRT bits
 * @param ep        Element parameters
 * @param keyPair   Key pair
 * @param numSlotsCKKS Number of CKKS slots
 * @param depth     Multiplicative depth
 * @param levelsAvailableBeforeBootstrap Levels before bootstrap
 * @param levelsComputation Levels for computation
 * @return Ciphertext<DCRTPoly> Encrypted 1 if value == 0, else encrypted 0
 */
lbcrypto::Ciphertext<lbcrypto::DCRTPoly> EqualBootstrapping(
    std::vector<lbcrypto::Poly> ctxtBFV,
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
    uint32_t levelsComputation);

#endif // MULTI_PRECISION_H
