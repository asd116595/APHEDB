#ifndef HOMLIKE_H
#define HOMLIKE_H

#include "equal_bootstrapping.h"

#include <vector>
#include <string>

using namespace lbcrypto;

/**
 * @brief Case 1: Direct match (no wildcards) – 1 sub-case
 */
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
    uint32_t levelsAvailableBeforeBootstrap, uint32_t levelsComputation);

/**
 * @brief Case 2: _ in front – 1 sub-case
 */
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
    uint32_t levelsAvailableBeforeBootstrap, uint32_t levelsComputation);

/**
 * @brief Case 3: _ in middle – 1 sub-case
 */
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
    uint32_t levelsAvailableBeforeBootstrap, uint32_t levelsComputation);

/**
 * @brief Case 4: _ in back – 1 sub-case
 */
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
    uint32_t levelsAvailableBeforeBootstrap, uint32_t levelsComputation);

/**
 * @brief Case 5: % in front – multiple sub-cases, results ADDED
 */
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
    uint32_t levelsAvailableBeforeBootstrap, uint32_t levelsComputation);

/**
 * @brief Case 6: % in middle – multiple sub-cases, results ADDED
 */
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
    uint32_t levelsAvailableBeforeBootstrap, uint32_t levelsComputation);

/**
 * @brief Case 7: % in back – multiple sub-cases, results ADDED
 */
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
    uint32_t levelsAvailableBeforeBootstrap, uint32_t levelsComputation);

#endif // HOMLIKE_H
