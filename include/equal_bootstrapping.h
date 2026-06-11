#ifndef EQUAL_BOOTSTRAPPING_H
#define EQUAL_BOOTSTRAPPING_H

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

extern const BigInteger QBFVINIT;
extern const BigInteger QBFVINITLARGE;

/**
 * @brief EqualBootstrapping function - performs equality check bootstrapping
 * 
 * @param ctxtBFV Input RLWE ciphertext
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
lbcrypto::Ciphertext<lbcrypto::DCRTPoly> EqualBootstrapping(
    std::vector<lbcrypto::Poly>& ctxtBFV, 
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

#endif // EQUAL_BOOTSTRAPPING_H
