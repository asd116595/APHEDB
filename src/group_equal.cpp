#define PROFILE
#include "group_equal.h"

using namespace lbcrypto;

Ciphertext<DCRTPoly> GroupEqual(
    const std::vector<Poly>& ctxtData,
    const std::vector<int64_t>& groupAssignments,
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
    uint32_t levelsComputation) {

    // Find unique groups (in order of first appearance)
    std::vector<int64_t> groups;
    for (auto g : groupAssignments) {
        bool found = false;
        for (auto existing : groups) {
            if (existing == g) { found = true; break; }
        }
        if (!found) groups.push_back(g);
    }

    std::cerr << "\n=== GroupEqual: processing " << groups.size() << " groups ===" << std::endl;

    // Convert RLWE data to CKKS for rotation operations
    // level = depth - (levelsAvailableBeforeBootstrap > 0)
    auto ctxtCKKS = SchemeletRLWEMP::ConvertRLWEToCKKS(
        *cc, ctxtData, keyPair.publicKey, QBFVInit, numSlotsCKKS,
        depth - (levelsAvailableBeforeBootstrap > 0));

    // Initialize result as encrypted 0
    std::vector<double> zeroVec(numSlots, 0);
    if (zeroVec.size() < numSlots) zeroVec = Fill<double>(zeroVec, numSlots);
    Plaintext ptxtZero = cc->MakeCKKSPackedPlaintext(zeroVec, 1, ctxtCKKS->GetLevel(), nullptr, numSlotsCKKS);
    auto result = cc->Encrypt(keyPair.publicKey, ptxtZero);

    for (size_t gIdx = 0; gIdx < groups.size(); ++gIdx) {
        int64_t groupId = groups[gIdx];

        // Find reference slot for this group
        int32_t refSlot = -1;
        for (size_t i = 0; i < groupAssignments.size(); ++i) {
            if (groupAssignments[i] == groupId) {
                refSlot = static_cast<int32_t>(i);
                break;
            }
        }
        if (refSlot < 0) continue;

        std::cerr << "  Group " << groupId << ": refSlot=" << refSlot << std::endl;

        // Create mask for this group: 1 for slots in this group, 0 for others
        std::vector<double> groupMask(numSlots, 0);
        for (size_t i = 0; i < groupAssignments.size() && i < numSlots; ++i) {
            if (groupAssignments[i] == groupId) {
                groupMask[i] = 1;
            }
        }

        // ===== Efficient broadcast using EvalSum =====
        // Step 1: Rotate to bring refSlot value to position 0
        Ciphertext<DCRTPoly> ctxtRotated;
        if (refSlot == 0) {
            ctxtRotated = ctxtCKKS;
        } else {
            ctxtRotated = cc->EvalRotate(ctxtCKKS, -refSlot);
        }

        // Step 2: Mask to keep only slot 0's value, zero out others
        std::vector<double> pos0Mask(numSlots, 0);
        pos0Mask[0] = 1;
        Plaintext ptxtPos0Mask = cc->MakeCKKSPackedPlaintext(pos0Mask, 1, ctxtRotated->GetLevel(), nullptr, numSlotsCKKS);
        auto ctxtRefAt0 = cc->EvalMult(ctxtRotated, ptxtPos0Mask);
        cc->ModReduceInPlace(ctxtRefAt0);

        // Step 3: EvalSum broadcasts the value at slot 0 to all slots
        // Since only slot 0 is non-zero, EvalSum puts ref_value into every slot
        // EvalSum uses O(log n) rotations, much more efficient than for-loop
        auto ctxtRefBroadcast = cc->EvalSum(ctxtRefAt0, numSlots);

        // Step 4: Rotate back so ref value aligns with original positions
        Ciphertext<DCRTPoly> ctxtRefAligned;
        if (refSlot == 0) {
            ctxtRefAligned = ctxtRefBroadcast;
        } else {
            ctxtRefAligned = cc->EvalRotate(ctxtRefBroadcast, refSlot);
        }

        // Step 5: diff = data - ref
        auto ctxtDiffCKKS = cc->EvalSub(ctxtCKKS, ctxtRefAligned);

        // Step 6: Convert diff CKKS to RLWE
        auto ctxtDiffRLWE = SchemeletRLWEMP::ConvertCKKSToRLWE(ctxtDiffCKKS, Q);

        // Step 7: EqualBootstrapping checks if each coefficient == 0
        auto eqResult = EqualBootstrapping(ctxtDiffRLWE, cc, QBFVInit, PInput, POutput, Q, Bigq,
                                            scaleTHI, order, numSlots, ringDim, dcrtBits,
                                            ep, keyPair, numSlotsCKKS, depth,
                                            levelsAvailableBeforeBootstrap, levelsComputation);

        // Step 8: Mask result with group mask
        Plaintext ptxtMask = cc->MakeCKKSPackedPlaintext(groupMask, 1, eqResult->GetLevel(), nullptr, numSlotsCKKS);
        auto ctxtMasked = cc->EvalMult(eqResult, ptxtMask);
        cc->ModReduceInPlace(ctxtMasked);

        // Accumulate into result
        result = cc->EvalAdd(result, ctxtMasked);
    }

    return result;
}


std::vector<Ciphertext<DCRTPoly>> GroupEqualMultiGroup(
    const std::vector<Poly>& ctxtData,
    const std::vector<int64_t>& groupAssignments,
    int64_t numGroups,
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
    uint32_t levelsComputation) {

    std::cerr << "\n=== GroupEqualMultiGroup: " << numGroups << " groups ===" << std::endl;

    std::vector<Ciphertext<DCRTPoly>> results;

    // Convert RLWE data to CKKS
    auto ctxtCKKS = SchemeletRLWEMP::ConvertRLWEToCKKS(
        *cc, ctxtData, keyPair.publicKey, QBFVInit, numSlotsCKKS,
        depth - (levelsAvailableBeforeBootstrap > 0));

    for (int64_t g = 0; g < numGroups; ++g) {
        // Find reference slot for group g
        int32_t refSlot = -1;
        for (size_t i = 0; i < groupAssignments.size(); ++i) {
            if (groupAssignments[i] == g) {
                refSlot = static_cast<int32_t>(i);
                break;
            }
        }
        if (refSlot < 0) {
            std::cerr << "  Group " << g << ": no members, skipping" << std::endl;
            continue;
        }

        std::cerr << "  Group " << g << ": refSlot=" << refSlot << std::endl;

        // Rotate to bring refSlot to position 0
        Ciphertext<DCRTPoly> ctxtRotated;
        if (refSlot == 0) {
            ctxtRotated = ctxtCKKS;
        } else {
            ctxtRotated = cc->EvalRotate(ctxtCKKS, -refSlot);
        }

        // Mask to keep only slot 0's value
        std::vector<double> pos0Mask(numSlots, 0);
        pos0Mask[0] = 1;
        Plaintext ptxtPos0Mask = cc->MakeCKKSPackedPlaintext(pos0Mask, 1, ctxtRotated->GetLevel(), nullptr, numSlotsCKKS);
        auto ctxtRefAt0 = cc->EvalMult(ctxtRotated, ptxtPos0Mask);
        cc->ModReduceInPlace(ctxtRefAt0);

        // EvalSum broadcasts slot 0's value to all slots
        auto ctxtRefBroadcast = cc->EvalSum(ctxtRefAt0, numSlots);

        // Rotate back
        Ciphertext<DCRTPoly> ctxtRefAligned;
        if (refSlot == 0) {
            ctxtRefAligned = ctxtRefBroadcast;
        } else {
            ctxtRefAligned = cc->EvalRotate(ctxtRefBroadcast, refSlot);
        }

        // diff = data - ref
        auto ctxtDiffCKKS = cc->EvalSub(ctxtCKKS, ctxtRefAligned);

        // Convert to RLWE and call EqualBootstrapping
        auto ctxtDiffRLWE = SchemeletRLWEMP::ConvertCKKSToRLWE(ctxtDiffCKKS, Q);

        auto eqResult = EqualBootstrapping(ctxtDiffRLWE, cc, QBFVInit, PInput, POutput, Q, Bigq,
                                            scaleTHI, order, numSlots, ringDim, dcrtBits,
                                            ep, keyPair, numSlotsCKKS, depth,
                                            levelsAvailableBeforeBootstrap, levelsComputation);

        results.push_back(eqResult);
    }

    return results;
}
