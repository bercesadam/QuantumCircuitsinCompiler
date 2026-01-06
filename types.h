#pragma once
#include <array>
#include "constexpr_complex.h"
#include "constexpr_pow2.h"

/// @file
/// @brief Core type aliases used across the project: sizes, complex number type, vectors and matrices.
///
/// @details
/// - `dimension_t` and `index_t` are the unsigned integral types used for sizes and indices.
/// - `cplx_t` is the project's constexpr-capable complex type.
/// - `state_vector_t<N>` is a std::array of `cplx_t` with N elements representing amplitudes.
/// - `matrix_t<R,C>` is a 2D std::array representing a matrix of complex amplitudes.
/// - `qbit_list_t<QBitCount>` is a fixed-size array of qubit indices used to specify affected qubits.
using dimension_t = std::size_t;
using index_t = std::size_t;

using cplx_t = ConstexprMath::Complex<double>;

/// @brief State vector with compile-time fixed size (array of complex amplitudes).
template<dimension_t StateCount>
using state_vector_t = std::array<cplx_t, StateCount>;

/// @brief Matrix type with compile-time fixed rows/columns.
/// @tparam Rows  Number of rows.
/// @tparam Cols  Number of columns.
template<dimension_t Rows, dimension_t Cols>
using matrix_t = std::array<std::array<cplx_t, Cols>, Rows>;

/// @brief Fixed-size list of qubit indices.
/// @tparam QBitCount  Number of qubits in the list.
template<index_t QBitCount>
using qbit_list_t = std::array<dimension_t, QBitCount>;
