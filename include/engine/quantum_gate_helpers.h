#pragma once
#include "types.h"


/// @file
/// @brief Helpers for applying gate matrices and validating matrix properties.
///
/**
 * @details
 * This header provides:
 *  - Type trait `is_gate_matrix` to identify gate matrix types.
 *  - `apply_unitary` : constexpr matrix-vector multiplication (used to apply a local gate).
 *  - `is_valid_square_matrix` : compile-time check for square matrices with power-of-two size.
 *  - `is_unitary` : constexpr runtime/checkable check that a matrix is unitary.
 */

/// @brief  Type trait to check if a type is a gate matrix:
/// square array of complex numbers where the dimension is a power of two.
template<typename T>
struct is_gate_matrix : std::false_type {};

template<dimension_t N>
struct is_gate_matrix<std::array<std::array<cplx_t, N>, N>> {
    static constexpr dimension_t dim = N;
	static constexpr bool value = ConstexprMath::is_power_of_two(N);
};

template<typename T>
inline constexpr bool is_gate_matrix_v = is_gate_matrix<T>::value;

/// @brief  Checks whether a provided square complex matrix is unitary.
/// @tparam Dim     Dimension of the square matrix (2^k).
/// @param mat      The matrix to check.
/// @return         True if mat * mat^† equals the identity, false otherwise.
///
/// @details
/// This function computes the conjugate transpose (Hermitian adjoint) of `mat`,
/// multiplies `mat` by its conjugate transpose and verifies that the product is
/// equal to the identity matrix within exact arithmetic. Because this routine
/// uses exact double checks for equality of re/im components, it is intended
/// for constexpr / compile-time constructed matrices used as gates.
template<auto M>
constexpr bool is_unitary()
{
	// Get the matrix dimension from the type trait
    constexpr dimension_t N = is_gate_matrix<std::remove_cvref_t<decltype(M)>>::dim;

	// Epsilon for floating-point comparison
    constexpr double epsilon = 1E-9;

	// Check unitarity: M * M^† == I
    for (dimension_t i = 0; i < N; i++)
        for (dimension_t j = 0; j < N; j++)
        {
            cplx_t sum{ 0,0 };
            for (dimension_t k = 0; k < N; k++)
                sum = sum + M[k][i].conj() * M[k][j];

            constexpr auto abs = [](double v) -> double
            {
                return (v < 0.0) ? -v : v;
            };

            if (i == j && (abs(sum.re - 1) > epsilon || abs(sum.im) > epsilon))
                return false;
            if (i != j && (abs(sum.re) > epsilon || abs(sum.im) > epsilon))
                return false;
        }
    return true;
}

/// @brief  Applies a unitary matrix to a state vector via matrix-vector multiplication.
template<dimension_t Dim>
constexpr state_vector_t<Dim>
apply_unitary(const matrix_t<Dim, Dim>& U,
    const state_vector_t<Dim>& v) noexcept
{
    // Result initialized to zero amplitudes
    state_vector_t<Dim> result{ cplx_t::zero() };

    // Standard matrix-vector product:
    // result[i] = sum_j U[i][j] * v[j]
    for (dimension_t i = 0; i < Dim; ++i)
    {
        // Accumulate the i-th output amplitude
        for (dimension_t j = 0; j < Dim; ++j)
        {
            // Multiply the matrix row by the vector and add
            result[i] += U[i][j] * v[j];
        }
    }

    return result;
}
