#ifndef GROUP_EQUAL_H
#define GROUP_EQUAL_H

#include "equal_bootstrapping.h"

#include <vector>
#include <string>
#include <utility>

using namespace lbcrypto;

/**
 * @brief GroupEqual - performs equality comparison within groups and returns 0/1 ciphertext string
 * 
 * Given encrypted data and a grouping scheme, this function:
 * 1. Splits data into groups (each group can be mapped to values 1/2/3)
 * 2. For each group, performs EqualBootstrapping to check if all elements in the group are equal
 * 3. Returns a ciphertext representing a binary string (0/1 for each position)
 * 
 * @param ctxtData Input RLWE ciphertext containing the data
 * @param groupAssignments Vector mapping each slot index to a group ID (0, 1, 2, ...)
 *        e.g., {0,0,1,1,2,2,0,0} means slots 0,1,6,7 belong to group 0; 2,3 to group 1; 4,5 to group 2
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
 * @return Ciphertext<DCRTPoly> Result ciphertext: each slot contains encrypted 1 if that slot's value
 *         equals the group's reference value, encrypted 0 otherwise
 */
lbcrypto::Ciphertext<lbcrypto::DCRTPoly> GroupEqual(
    const std::vector<lbcrypto::Poly>& ctxtData,
    const std::vector<int64_t>& groupAssignments,
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

/**
 * @brief GroupEqualMultiGroup - compares each slot's value against multiple group reference values
 * 
 * For each group, a reference value is extracted. Then each slot is compared against
 * all group reference values. Returns a 2D result where result[i][j] = 1 if slot i's value
 * equals group j's reference value, else 0.
 * 
 * This is useful for: given data split into groups, determine which group each value belongs to.
 * 
 * @param ctxtData Input RLWE ciphertext containing the data
 * @param groupAssignments Vector mapping each slot index to a group ID
 * @param numGroups Total number of groups
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
 * @return std::vector<Ciphertext<DCRTPoly>> Vector of ciphertexts, one per group.
 *         result[j][i] = encrypted 1 if slot i's value equals group j's reference value
 */
std::vector<lbcrypto::Ciphertext<lbcrypto::DCRTPoly>> GroupEqualMultiGroup(
    const std::vector<lbcrypto::Poly>& ctxtData,
    const std::vector<int64_t>& groupAssignments,
    int64_t numGroups,
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

#endif // GROUP_EQUAL_H
