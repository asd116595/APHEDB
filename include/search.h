#ifndef SEARCH_H
#define SEARCH_H

#include "equal_bootstrapping.h"

#include <vector>

using namespace lbcrypto;

/**
 * @brief Search function - performs search operation using EqualBootstrapping
 * 
 * @param ctxtStr Input RLWE ciphertext
 * @param ctxtstr Input CKKS ciphertext
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
 * @return Ciphertext<DCRTPoly> Result ciphertext
 */
lbcrypto::Ciphertext<lbcrypto::DCRTPoly> Search(
    std::vector<lbcrypto::Poly>& ctxtStr,
    std::vector<std::vector<lbcrypto::Poly>> ctxtstr,
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

#endif // SEARCH_H
