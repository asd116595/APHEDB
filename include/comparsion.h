#ifndef COMPARSION_H
#define COMPARSION_H

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
 * @brief SignCipher function - performs sign evaluation on ciphertext
 * 
 * Homomorphically evaluates the sign of a ciphertext by extracting digits
 * and applying functional bootstrapping.
 * Returns sign(diff) where sign(x) = 1 if x > 0, else 0.
 * 
 * @param ctxtBFV Input RLWE ciphertext (difference of two encrypted values), passed by value
 * @param cc CKKS crypto context
 * @param QBFVInit Initial modulus for BFV
 * @param PInput Input plaintext modulus
 * @param PDigit Digit modulus for the sign evaluation
 * @param Q RLWE modulus
 * @param Bigq CKKS scaling factor
 * @param scaleTHI Scaling factor for homomorphic decoding
 * @param scaleStepTHI Scaling factor for step function in homomorphic decoding
 * @param order Order of Hermite interpolation
 * @param numSlots Number of slots
 * @param ringDim Ring dimension
 * @param dcrtBits DCRT bits
 * @param ep Element parameters
 * @param keyPair Key pair
 * @param numSlotsCKKS Number of CKKS slots
 * @param depth Multiplicative depth
 * @param levelsAvailableBeforeBootstrap Levels available before bootstrap
 * @param levelsComputation Levels for computation
 * @return Ciphertext<DCRTPoly> Result ciphertext representing the sign
 */
lbcrypto::Ciphertext<lbcrypto::DCRTPoly> SignCipher(
    std::vector<lbcrypto::Poly> ctxtBFV,   // passed by VALUE to avoid modifying caller's data
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
 * @brief GreaterThan - checks if x > y homomorphically
 * 
 * Computes sign(x - y): returns encrypted 1 if x > y, else encrypted 0.
 * 
 * @param ctxtX Encrypted x
 * @param ctxtY Encrypted y
 * @param cc CKKS crypto context
 * @param QBFVInit Initial modulus for BFV
 * @param PInput Input plaintext modulus
 * @param PDigit Digit modulus for the sign evaluation
 * @param Q RLWE modulus
 * @param Bigq CKKS scaling factor
 * @param scaleTHI Scaling factor for homomorphic decoding
 * @param scaleStepTHI Scaling factor for step function in homomorphic decoding
 * @param order Order of Hermite interpolation
 * @param numSlots Number of slots
 * @param ringDim Ring dimension
 * @param dcrtBits DCRT bits
 * @param ep Element parameters
 * @param keyPair Key pair
 * @param numSlotsCKKS Number of CKKS slots
 * @param depth Multiplicative depth
 * @param levelsAvailableBeforeBootstrap Levels available before bootstrap
 * @param levelsComputation Levels for computation
 * @return Ciphertext<DCRTPoly> Encrypted 1 if x > y, else encrypted 0
 */
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
    uint32_t levelsComputation);

/**
 * @brief GreaterEqual - checks if x >= y homomorphically
 * 
 * Computes 1 - sign(y - x): returns encrypted 1 if x >= y, else encrypted 0.
 * 
 * @param ctxtX Encrypted x
 * @param ctxtY Encrypted y
 * @param cc CKKS crypto context
 * @param QBFVInit Initial modulus for BFV
 * @param PInput Input plaintext modulus
 * @param PDigit Digit modulus for the sign evaluation
 * @param Q RLWE modulus
 * @param Bigq CKKS scaling factor
 * @param scaleTHI Scaling factor for homomorphic decoding
 * @param scaleStepTHI Scaling factor for step function in homomorphic decoding
 * @param order Order of Hermite interpolation
 * @param numSlots Number of slots
 * @param ringDim Ring dimension
 * @param dcrtBits DCRT bits
 * @param ep Element parameters
 * @param keyPair Key pair
 * @param numSlotsCKKS Number of CKKS slots
 * @param depth Multiplicative depth
 * @param levelsAvailableBeforeBootstrap Levels available before bootstrap
 * @param levelsComputation Levels for computation
 * @return Ciphertext<DCRTPoly> Encrypted 1 if x >= y, else encrypted 0
 */
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
    uint32_t levelsComputation);

/**
 * @brief LessThan - checks if x < y homomorphically
 * 
 * Computes sign(y - x): returns encrypted 1 if x < y, else encrypted 0.
 * 
 * @param ctxtX Encrypted x
 * @param ctxtY Encrypted y
 * @param cc CKKS crypto context
 * @param QBFVInit Initial modulus for BFV
 * @param PInput Input plaintext modulus
 * @param PDigit Digit modulus for the sign evaluation
 * @param Q RLWE modulus
 * @param Bigq CKKS scaling factor
 * @param scaleTHI Scaling factor for homomorphic decoding
 * @param scaleStepTHI Scaling factor for step function in homomorphic decoding
 * @param order Order of Hermite interpolation
 * @param numSlots Number of slots
 * @param ringDim Ring dimension
 * @param dcrtBits DCRT bits
 * @param ep Element parameters
 * @param keyPair Key pair
 * @param numSlotsCKKS Number of CKKS slots
 * @param depth Multiplicative depth
 * @param levelsAvailableBeforeBootstrap Levels available before bootstrap
 * @param levelsComputation Levels for computation
 * @return Ciphertext<DCRTPoly> Encrypted 1 if x < y, else encrypted 0
 */
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
    uint32_t levelsComputation);

/**
 * @brief LessEqual - checks if x <= y homomorphically
 * 
 * Computes 1 - sign(x - y): returns encrypted 1 if x <= y, else encrypted 0.
 * 
 * @param ctxtX Encrypted x
 * @param ctxtY Encrypted y
 * @param cc CKKS crypto context
 * @param QBFVInit Initial modulus for BFV
 * @param PInput Input plaintext modulus
 * @param PDigit Digit modulus for the sign evaluation
 * @param Q RLWE modulus
 * @param Bigq CKKS scaling factor
 * @param scaleTHI Scaling factor for homomorphic decoding
 * @param scaleStepTHI Scaling factor for step function in homomorphic decoding
 * @param order Order of Hermite interpolation
 * @param numSlots Number of slots
 * @param ringDim Ring dimension
 * @param dcrtBits DCRT bits
 * @param ep Element parameters
 * @param keyPair Key pair
 * @param numSlotsCKKS Number of CKKS slots
 * @param depth Multiplicative depth
 * @param levelsAvailableBeforeBootstrap Levels available before bootstrap
 * @param levelsComputation Levels for computation
 * @return Ciphertext<DCRTPoly> Encrypted 1 if x <= y, else encrypted 0
 */
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
    uint32_t levelsComputation);

#endif // COMPARSION_H
