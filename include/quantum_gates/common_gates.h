#pragma once
#include "core_types.h"
#include "constexprmath/constexpr_core_functions.h"


/// @file
/// @brief Common single- and two-qubit gate matrices and helpers.
///
/**
 * @details
 * Provides common quantum gate matrices
 * and the `identity_matrix<QBitCount>()` helper for arbitrary qubit counts.
 */
namespace QCC
{
    namespace Gates
    {
        /// @brief Mathematical constants used for gate definitions
        constexpr float_t sqrt2 = 1.41421356237309505;
        constexpr float_t inv_sqrt2 = 1.0 / sqrt2;

        /**
         * @brief Produce an identity matrix for `QBitCount` qubits (2^QBitCount � 2^QBitCount).
         * @tparam QBitCount  Number of qubits.
         * @return A diagonal identity matrix with ones on the main diagonal.
         */
        template<Ket::dimension_t QBitCount>
        constexpr Ket::matrix_t<ConstexprMath::pow2(QBitCount)> identityMatrix() noexcept
        {
            constexpr dimension_t Dim = ConstexprMath::pow2(QBitCount);

            // Zero-initialize and set the diagonal entries to 1.0
            matrix_t<Dim> Identity = {};
            for (dimension_t i = 0; i < Dim; ++i)
            {
                Identity[i][i] = cplx_t::fromReal(1.0);
            }
            return Identity;
        }

        // Single-qubit identity
        constexpr Ket::matrix_t<2> I = identityMatrix<1U>();

        // Hadamard gate: H = (1/sqrt(2)) * [[1, 1], [1, -1]]
        constexpr Ket::matrix_t<2> H = { {
            { cplx_t(1.0 / sqrt2, 0.0), cplx_t(1.0 / sqrt2, 0.0) },
            { cplx_t(1.0 / sqrt2, 0.0), cplx_t(-1.0 / sqrt2, 0.0) }
        } };

        // Pauli-X (NOT)
        constexpr Ket::matrix_t<2> X = { {
            { cplx_t(0.0, 0.0), cplx_t(1.0, 0.0) },
            { cplx_t(1.0, 0.0), cplx_t(0.0, 0.0) }
        } };

        // Pauli-Y
        constexpr Ket::matrix_t<2> Y = { {
            { cplx_t(0.0, 0.0), cplx_t(0.0, -1.0) },
            { cplx_t(0.0, 1.0), cplx_t(0.0, 0.0) }
        } };

        // Pauli-Z
        constexpr Ket::matrix_t<2> Z = { {
            { cplx_t(1.0, 0.0), cplx_t(0.0, 0.0) },
            { cplx_t(0.0, 0.0), cplx_t(-1.0, 0.0) }
        } };

        // 2-qubit SWAP gate
        constexpr Ket::matrix_t<4> SWAP = { {
            { cplx_t(1.0, 0.0), cplx_t(0.0, 0.0), cplx_t(0.0, 0.0), cplx_t(0.0, 0.0) },
            { cplx_t(0.0, 0.0), cplx_t(0.0, 0.0), cplx_t(1.0, 0.0), cplx_t(0.0, 0.0) },
            { cplx_t(0.0, 0.0), cplx_t(1.0, 0.0), cplx_t(0.0, 0.0), cplx_t(0.0, 0.0) },
            { cplx_t(0.0, 0.0), cplx_t(0.0, 0.0), cplx_t(0.0, 0.0), cplx_t(1.0, 0.0) }
        } };

        // Control-first CNOT gate (CX)
        constexpr Ket::matrix_t<4> CX = { {
            { cplx_t(1.0, 0.0), cplx_t(0.0, 0.0), cplx_t(0.0, 0.0), cplx_t(0.0, 0.0) },
            { cplx_t(0.0, 0.0), cplx_t(0.0, 0.0), cplx_t(0.0, 0.0), cplx_t(1.0, 0.0) },
            { cplx_t(0.0, 0.0), cplx_t(0.0, 0.0), cplx_t(1.0, 0.0), cplx_t(0.0, 0.0) },
            { cplx_t(0.0, 0.0), cplx_t(1.0, 0.0), cplx_t(0.0, 0.0), cplx_t(0.0, 0.0) }
        } };
        constexpr auto CNOT = CX;

        // Control-first Toffoli gate (CCX)
        constexpr Ket::matrix_t<8> CCX = { {
            { cplx_t(1.0,0.0), cplx_t(0.0,0.0), cplx_t(0.0,0.0), cplx_t(0.0,0.0),
              cplx_t(0.0,0.0), cplx_t(0.0,0.0), cplx_t(0.0,0.0), cplx_t(0.0,0.0) },

            { cplx_t(0.0,0.0), cplx_t(1.0,0.0), cplx_t(0.0,0.0), cplx_t(0.0,0.0),
              cplx_t(0.0,0.0), cplx_t(0.0,0.0), cplx_t(0.0,0.0), cplx_t(0.0,0.0) },

            { cplx_t(0.0,0.0), cplx_t(0.0,0.0), cplx_t(1.0,0.0), cplx_t(0.0,0.0),
              cplx_t(0.0,0.0), cplx_t(0.0,0.0), cplx_t(0.0,0.0), cplx_t(0.0,0.0) },

            { cplx_t(0.0,0.0), cplx_t(0.0,0.0), cplx_t(0.0,0.0), cplx_t(1.0,0.0),
              cplx_t(0.0,0.0), cplx_t(0.0,0.0), cplx_t(0.0,0.0), cplx_t(0.0,0.0) },

            { cplx_t(0.0,0.0), cplx_t(0.0,0.0), cplx_t(0.0,0.0), cplx_t(0.0,0.0),
              cplx_t(1.0,0.0), cplx_t(0.0,0.0), cplx_t(0.0,0.0), cplx_t(0.0,0.0) },

            { cplx_t(0.0,0.0), cplx_t(0.0,0.0), cplx_t(0.0,0.0), cplx_t(0.0,0.0),
              cplx_t(0.0,0.0), cplx_t(1.0,0.0), cplx_t(0.0,0.0), cplx_t(0.0,0.0) },

            { cplx_t(0.0,0.0), cplx_t(0.0,0.0), cplx_t(0.0,0.0), cplx_t(0.0,0.0),
              cplx_t(0.0,0.0), cplx_t(0.0,0.0), cplx_t(0.0,0.0), cplx_t(1.0,0.0) },

            { cplx_t(0.0,0.0), cplx_t(0.0,0.0), cplx_t(0.0,0.0), cplx_t(0.0,0.0),
              cplx_t(0.0,0.0), cplx_t(0.0,0.0), cplx_t(1.0,0.0), cplx_t(0.0,0.0) }
        } };
        constexpr auto TOFFOLI = CCX;



    } // namespace Gates    
}