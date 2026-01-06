#pragma once
#include "gates.h"

/// @file
/// @brief Construction helpers for controlled-NOT (C^n X) style gates.
///
/**
 * @details
 * This header defines templated helpers to construct controlled-X matrices
 * where a single target qubit is flipped when all control qubits are 1.
 * It also declares common aliases: `CX` / `CNOT`, `CCX` / `TOFFOLI`, `CCCX`, etc.
 */

namespace Gates
{
    /**
     * @brief  Constructs the C^n X (multi-control X) matrix for QBitCount qubits.
     *
     * @tparam QBitCount  Number of qubits for which the resulting matrix is sized (controls + target).
     * @return A unitary matrix of size 2^QBitCount × 2^QBitCount implementing the multi-controlled-X.
     *
     * @remarks
     *  - The construction starts from the identity and swaps the two basis states
     *    corresponding to all-control qubits being 1 with the target qubit being 0/1.
     *  - For QBitCount == 2 this is the standard CNOT. For QBitCount == 3 it is Toffoli.
     */
    template<dimension_t QBitCount>
    constexpr matrix_t<ConstexprMath::pow2(QBitCount), ConstexprMath::pow2(QBitCount)>
    make_CnX_matrix() noexcept
    {
        constexpr dimension_t Dim = ConstexprMath::pow2(QBitCount);

        // Start from the identity matrix of the given qubit count
        matrix_t<Dim, Dim> M = identity_matrix<QBitCount>();

        // The two basis indices to swap are the last two states:
        //   i0 = |11...10>  (second-to-last)
        //   i1 = |11...11>  (last)
        constexpr dimension_t i0 = Dim - 2; // |111...10>
        constexpr dimension_t i1 = Dim - 1; // |111...11>

        // Replace the identity entries in the bottom-right corner with a 2x2 X
        // (i.e., zero the diagonal entries for these two states)
        M[i0][i0] = cplx_t::fromReal(0.0);
        M[i1][i1] = cplx_t::fromReal(0.0);

        // Place the X swap in the corner (flip |...10> <-> |...11>)
        M[i0][i1] = cplx_t::fromReal(1.0);
        M[i1][i0] = cplx_t::fromReal(1.0);

        return M;
    }

    // Common named instances for convenience
    constexpr auto CX = make_CnX_matrix<dimension_t(2)>();
    constexpr auto CNOT = CX;
    constexpr auto CCX = make_CnX_matrix<dimension_t(3)>();
    constexpr auto TOFFOLI = CCX;
    constexpr auto CCCX = make_CnX_matrix<dimension_t(4)>();
}