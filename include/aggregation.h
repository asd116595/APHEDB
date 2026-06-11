#ifndef AGGREGATION_H
#define AGGREGATION_H

// NOTE: comparsion.h and equal_bootstrapping.h both include
// "schemelet/rlwe-mp.h" which has NO include guards.
// Use a manual guard to prevent double-inclusion.
#ifndef SCHEMELET_RLWE_MP_GUARD
#define SCHEMELET_RLWE_MP_GUARD
#include "schemelet/rlwe-mp.h"
#endif

#include "comparsion.h"
#include "equal_bootstrapping.h"

#include <vector>

using namespace lbcrypto;

/**
 * @brief Homomorphic SUM: sums all slot values in CKKS domain
 *
 * Converts RLWE to CKKS, then uses EvalSum to add all slots together.
 * Result: each slot contains the total sum.
 *
 * For many values (e.g. 1000), ringDim=4096 supports up to 4096 coefficients.
 *
 * @param ctxtData Input RLWE ciphertext (coefficient-packed)
 * @param cc CKKS crypto context
 * @param QBFVInit Initial modulus for BFV
 * @param Q RLWE modulus
 * @param Bigq CKKS scaling factor
 * @param keyPair Key pair
 * @param numSlotsCKKS Number of CKKS slots
 * @param numSlots Number of slots
 * @param depth Multiplicative depth
 * @param levelsAvailableBeforeBootstrap Levels available before bootstrap
 * @return Ciphertext<DCRTPoly> Result ciphertext (CKKS domain, sum in all slots)
 */
Ciphertext<DCRTPoly> HomSum(
    const std::vector<Poly>& ctxtData,
    CryptoContextCKKSRNS::ContextType& cc,
    BigInteger QBFVInit,
    BigInteger Q,
    BigInteger Bigq,
    KeyPair<DCRTPoly>& keyPair,
    uint32_t numSlotsCKKS,
    uint32_t numSlots,
    uint32_t depth,
    uint32_t levelsAvailableBeforeBootstrap);

/**
 * @brief Homomorphic COUNT: counts non-zero values across all slots
 *
 * Uses EqualBootstrapping to check each value != 0, then sums the results.
 * Returns CKKS ciphertext where each slot = count of non-zero values.
 *
 * For many values, ringDim=4096 supports up to 4096 coefficients per ciphertext.
 *
 * @param ctxtData Input RLWE ciphertext (coefficient-packed)
 * @param cc CKKS crypto context
 * @param QBFVInit Initial modulus for BFV
 * @param PInput Input plaintext modulus
 * @param POutput Output plaintext modulus
 * @param Q RLWE modulus
 * @param Bigq CKKS scaling factor
 * @param scaleTHI Scaling factor for homomorphic decoding
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
 * @return Ciphertext<DCRTPoly> Result ciphertext (CKKS domain, count in all slots)
 */
Ciphertext<DCRTPoly> HomCount(
    const std::vector<Poly>& ctxtData,
    CryptoContextCKKSRNS::ContextType& cc,
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
    std::shared_ptr<M4DCRTParams>& ep,
    KeyPair<DCRTPoly>& keyPair,
    uint32_t numSlotsCKKS,
    uint32_t depth,
    uint32_t levelsAvailableBeforeBootstrap,
    uint32_t levelsComputation);

/**
 * @brief Homomorphic AVG: computes encrypted SUM and encrypted COUNT
 *
 * Returns a pair: (sumCipher, countCipher).
 * Average = sum / count can be computed after decryption.
 *
 * @param ctxtData Input RLWE ciphertext (coefficient-packed)
 * @param cc CKKS crypto context
 * @param QBFVInit Initial modulus for BFV
 * @param PInput Input plaintext modulus
 * @param POutput Output plaintext modulus
 * @param Q RLWE modulus
 * @param Bigq CKKS scaling factor
 * @param scaleTHI Scaling factor for homomorphic decoding
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
 * @return std::pair<Ciphertext, Ciphertext> (sumResult, countResult)
 */
std::pair<Ciphertext<DCRTPoly>, Ciphertext<DCRTPoly>> HomAvg(
    const std::vector<Poly>& ctxtData,
    CryptoContextCKKSRNS::ContextType& cc,
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
    std::shared_ptr<M4DCRTParams>& ep,
    KeyPair<DCRTPoly>& keyPair,
    uint32_t numSlotsCKKS,
    uint32_t depth,
    uint32_t levelsAvailableBeforeBootstrap,
    uint32_t levelsComputation);

/**
 * @brief Homomorphic MAX: finds maximum value across all slots
 *
 * Tournament-style reduction in RLWE coefficient domain:
 * Uses RLWE-level subtraction and SignCipher-based comparison.
 *
 * Strategy (works for up to ringDim values):
 *   For each tournament step s = 1, 2, 4, 8, ...:
 *     - Split: extract pairs using coefficient rotation (EvalAutomorphism)
 *     - Compare: compute diff in RLWE, call GreaterThan
 *     - Select: max = sel * a + (1-sel) * b in CKKS domain
 *     - Convert result back to RLWE for next step
 *
 * Result: RLWE ciphertext with the maximum value in all positions.
 *
 * @param ctxtData Input RLWE ciphertext (coefficient-packed)
 * @param cc CKKS crypto context
 * @param QBFVInit Initial modulus for BFV
 * @param PInput Input plaintext modulus
 * @param PDigit Digit modulus for sign evaluation
 * @param Q RLWE modulus
 * @param Bigq CKKS scaling factor
 * @param scaleTHI Scaling factor for homomorphic decoding
 * @param scaleStepTHI Scaling factor for step function
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
 * @return Ciphertext<DCRTPoly> Result ciphertext (RLWE domain, max in all slots)
 */
Ciphertext<DCRTPoly> HomMax(
    const std::vector<Poly>& ctxtData,
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
    uint32_t levelsComputation);

/**
 * @brief Homomorphic MIN: finds minimum value across all slots
 *
 * Same tournament approach as HomMax, but uses LessThan as selector.
 *
 * @param ctxtData Input RLWE ciphertext (coefficient-packed)
 * @param cc CKKS crypto context
 * @param QBFVInit Initial modulus for BFV
 * @param PInput Input plaintext modulus
 * @param PDigit Digit modulus for sign evaluation
 * @param Q RLWE modulus
 * @param Bigq CKKS scaling factor
 * @param scaleTHI Scaling factor for homomorphic decoding
 * @param scaleStepTHI Scaling factor for step function
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
 * @return Ciphertext<DCRTPoly> Result ciphertext (RLWE domain, min in all slots)
 */
Ciphertext<DCRTPoly> HomMin(
    const std::vector<Poly>& ctxtData,
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
    uint32_t levelsComputation);

#endif // AGGREGATION_H
