#pragma once
#include "common_gates.h"
#include "constexprmath/constexpr_trigon.h"


namespace Gates
{
    /**
     * @brief     Construct the inverse Quantum Fourier Transform (QFT†) dense matrix for QBitCount qubits.
     *
     * @tparam QBitCount  Number of qubits (matrix dimension = 2^QBitCount).
     * @return A unitary matrix implementing the inverse QFT:
     *         U[j,k] = (1/sqrt(N)) * exp(-2π i j k / N)
     *
     * @note This routine is constexpr-friendly and uses the project's constexpr trigonometric
     *       helpers for compile-time matrix construction.
     */
    template<dimension_t QBitCount>
    constexpr matrix_t<ConstexprMath::pow2(QBitCount), ConstexprMath::pow2(QBitCount)>
        make_IQFT_matrix() noexcept
    {
        constexpr dimension_t Dim = ConstexprMath::pow2(QBitCount);
        constexpr double InvSqrtDim = 1.0 / ConstexprMath::sqrt(static_cast<double>(Dim));

        matrix_t<Dim, Dim> IQFTMatrix{};

        for (dimension_t j = 0; j < Dim; ++j)
        {
            for (dimension_t k = 0; k < Dim; ++k)
            {
                // angle = 2π * j * k / N
                const double Angle =
                    2.0 * ConstexprMath::Pi * static_cast<double>(j * k)
                    / static_cast<double>(Dim);

                // For inverse QFT use the negative exponent:
                // exp(-i * angle) = cos(angle) - i * sin(angle)
                const double re = ConstexprMath::cos(Angle);
                const double im = -ConstexprMath::sin(Angle);

                IQFTMatrix[j][k] = cplx_t(
                    InvSqrtDim * re,
                    InvSqrtDim * im
                );
            }
        }

        return IQFTMatrix;
    }
}