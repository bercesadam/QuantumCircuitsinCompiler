#pragma once
#include "types.h"


/// @file
/// @brief Helpers for applying gate matrices and validating matrix properties.
///
/**
 * @details
 * This header provides:
 *  - `apply_unitary` : constexpr matrix-vector multiplication (used to apply a local gate).
 *  - `is_valid_square_matrix` : compile-time check for square matrices with power-of-two size.
 *  - `is_unitary` : constexpr runtime/checkable check that a matrix is unitary.
 */

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

/// @brief  Checks if a matrix type is a valid square matrix whose dimensions are powers of two.
/// @tparam MatrixType  The matrix type to check. It is expected to expose `size()` and `at<0>().size()`.
/// @return     True if the matrix is square and both dimensions are non-zero powers of two, false otherwise.
///
/// @note This is a compile-time check (consteval). It validates shape constraints only; it does not inspect
///       element values.
template<dimension_t Dim>
constexpr bool is_valid_square_matrix(const matrix_t<Dim, Dim>& mat) noexcept
{
    constexpr std::size_t rows = mat.size();

    // If there are no rows or rows is not a power of two the type is invalid
    if constexpr ((rows == 0) || (!ConstexprMath::is_power_of_two(rows)))
    {
        return false;
    }
    else
    {
        // Inspect columns via the first row's container size
        constexpr std::size_t cols = mat[0].size();

        // Columns must be non-zero and a power of two
        if constexpr ((cols == 0) || !ConstexprMath::is_power_of_two(cols))
        {
            return false;
        }

        // Valid square iff rows == cols
        return rows == cols;
    }
}

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
template<dimension_t Dim>
constexpr bool is_unitary(const matrix_t<Dim, Dim>& mat) noexcept
{
    // Compute the conjugate transpose (mat^†)
    matrix_t<Dim, Dim> conj_transpose{};
    for (dimension_t i = 0; i < Dim; ++i) {
        for (dimension_t j = 0; j < Dim; ++j) {
            // element (j,i) of the conjugate transpose = conjugate of (i,j)
			//conj_transpose[j][i] = cplx_t(mat[i][j].re, -mat[i][j].im);
        }
    }   

    // Compute product = mat * mat^†
    matrix_t<Dim, Dim> product{};
    for (dimension_t i = 0; i < Dim; ++i) {
        for (dimension_t j = 0; j < Dim; ++j) {
            cplx_t sum = cplx_t::zero();

            // Accumulate the dot product of row i of mat and column j of conj_transpose
            for (dimension_t k = 0; k < Dim; ++k) {
                //sum += mat[i][k] * conj_transpose[k][j];
            }

            product[i][j] = sum;
        }
    }

    // Verify that the product equals the identity matrix exactly
    for (dimension_t i = 0; i < Dim; ++i) {
        for (dimension_t j = 0; j < Dim; ++j) {
            if (i == j) {
                // Diagonal must be 1 + 0i
                if (product[i][j].re != 1.0 || product[i][j].im != 0.0) {
                    return false;
                }
            } else {
                // Off-diagonals must be 0 + 0i
                if (product[i][j].re != 0.0 || product[i][j].im != 0.0) {
                    return false;
                }
            }
        }
    }
    return true;
}
