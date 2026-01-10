#pragma once
#include "gates.h"
#include "constexpr_trigon.h"


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
        constexpr dimension_t N = ConstexprMath::pow2(QBitCount);
        constexpr double inv_sqrt_N = 1.0 / ConstexprMath::sqrt(static_cast<double>(N));

        matrix_t<N, N> M{};

        for (dimension_t j = 0; j < N; ++j)
        {
            for (dimension_t k = 0; k < N; ++k)
            {
                // angle = 2π * j * k / N
                const double angle =
                    2.0 * ConstexprMath::Pi * static_cast<double>(j * k)
                    / static_cast<double>(N);

                // For inverse QFT use the negative exponent:
                // exp(-i * angle) = cos(angle) - i * sin(angle)
                const double re = ConstexprMath::cos_constexpr(angle);
                const double im = -ConstexprMath::sin_constexpr(angle);

                M[j][k] = cplx_t(
                    inv_sqrt_N * re,
                    inv_sqrt_N * im
                );
            }
        }

        return M;
    }
}