#pragma once

#include "core_types.h"
#include "hamiltonian/hamiltonian.h"
#include "wavefunction/state_vector.h"


/// @file
/// @brief Crank–Nicolson solver for the time-dependent Schrödinger equation.
///
/// @details
/// This file implements a numerical time evolution scheme for quantum
/// state vectors based on the Crank–Nicolson method. The solver advances
/// a discretized wavefunction by solving the time-dependent Schrödinger
/// equation (TDSE)
///
///   iħ ∂ψ(t)/∂t = H ψ(t)
///
/// where ψ(t) is a state vector in a finite-dimensional Hilbert space and
/// H is the Hamiltonian operator of the system.
///
/// The Hamiltonian is provided explicitly as a matrix representation,
/// typically originating from a spatial discretization of the kinetic
/// and potential energy operators (e.g. finite-difference approximation
/// of the Laplacian plus a potential term).
///
/// Time integration is performed using the Crank–Nicolson scheme, which
/// corresponds to an implicit midpoint discretization in time:
///
///   ( I + i·Δt/(2ħ) · H ) · ψⁿ⁺¹ = ( I - i·Δt/(2ħ) · H ) · ψⁿ
///
/// This scheme is:
///  - second-order accurate in time,
///  - unconditionally stable,
///  - norm-preserving for Hermitian Hamiltonians,
///  - and therefore suitable for unitary quantum time evolution.
///
/// In this implementation, the Hamiltonian is assumed to be time-independent
/// (or piecewise constant in time). As a consequence, the Crank–Nicolson
/// system matrices are constructed once and reused for each time step.
///
/// If the Hamiltonian matrix is tridiagonal—as is the case for many
/// one-dimensional finite-difference discretizations—the resulting linear
/// system can be solved efficiently in O(N) time using the Thomas algorithm.


/// @brief Compact storage representation of a tridiagonal matrix.
/// as we have useful information only in the three non-zero diagonals.
template<dimension_t Dim>
struct tridiagonal_matrix_t
{
	/// @brief  Lower diagonal (a₁ ... a_{N−1})
	/// @details
	/// a[0] is unused.
	std::array<cplx_t, Dim> a{};

	/// @brief  Main diagonal (b₀ ... b_{N−1})
	std::array<cplx_t, Dim> b{};

	/// @brief  Upper diagonal (c₀ ... c_{N−2})
	/// @details
	/// c[Dim − 1] is unused.
	std::array<cplx_t, Dim> c{};
};

/// @brief  Helper function to construct the Crank–Nicolson system matrices A and B.
/// @tparam Dim     Dimension of the Hilbert space.
/// @param  hamiltonian  Hamiltonian operator of the system.
/// @param  dt           Time step size.
/// @param  A            Output matrix A = I + i·dt/(2ħ)·H.
/// @param  B            Output matrix B = I - i·dt/(2ħ)·H.
///
/// @details
/// This function builds the two matrices required by the Crank–Nicolson
/// time integration scheme. The matrices arise from the implicit midpoint
/// discretization of the time-dependent Schrödinger equation.
///
/// If the Hamiltonian matrix is tridiagonal, both A and B remain
/// tridiagonal, enabling efficient O(N) time stepping.
template<dimension_t Dim>
static constexpr void buildCrankNicolsonMatrices(const Hamiltonian<Dim>& hamiltonian, double dt,
	tridiagonal_matrix_t<Dim>& A, tridiagonal_matrix_t<Dim>& B) noexcept
{
	const auto& H = hamiltonian.getMatrix();
	const double hBar = hamiltonian.getConstants().hBar;

	// i * dt / (2ħ)
	const cplx_t Factor(0.0, dt / (2.0 * hBar));

	for (dimension_t i = 0; i < Dim; ++i)
	{
		// Build main diagonal
		A.b[i] = cplx_t::fromReal(1.0) + Factor * H[i][i];
		B.b[i] = cplx_t::fromReal(1.0) - Factor * H[i][i];

		//  Build lower diagonal
		if (i > 0)
		{
			A.a[i] = Factor * H[i][i - 1];
			B.a[i] = -Factor * H[i][i - 1];
		}

		//  Build upper diagonal
		if (i + 1 < Dim)
		{
			A.c[i] = Factor * H[i][i + 1];
			B.c[i] = -Factor * H[i][i + 1];
		}
	}
}

/// @brief  Helper function to compute the product of an instance of the helper tridiagonal matrix type
///         and a state vector.
/// @tparam Dim     Dimension of the vector space.
/// @param  M       Tridiagonal matrix.
/// @param  x       Input vector.
/// @return         Resulting vector M · x.
///
/// @details
/// This routine performs an efficient matrix–vector multiplication
/// exploiting the tridiagonal structure of the matrix. It is primarily
/// used to construct the right-hand side of the Crank–Nicolson system.
template<dimension_t Dim>
static constexpr StateVector<Dim>
multiplyTrigiagonal(const tridiagonal_matrix_t<Dim>& M,	const StateVector<Dim>& x) noexcept
{
	StateVector<Dim> Result{ cplx_t::zero() };

	for (dimension_t i = 0; i < Dim; ++i)
	{
		// Main diagonal contribution
		Result[i] += M.b[i] * x[i];

		// Lower diagonal contribution
		if (i > 0)
			Result[i] += M.a[i] * x[i - 1];

		// Upper diagonal contribution
		if (i + 1 < Dim)
			Result[i] += M.c[i] * x[i + 1];
	}

	return Result;
}

/// @brief  Solves a tridiagonal linear system using the Thomas algorithm.
/// @tparam Dim     Dimension of the linear system.
/// @param  M       Tridiagonal coefficient matrix.
/// @param  d       Right-hand-side vector.
/// @return         Solution vector x satisfying M · x = d.
///
/// @details
/// This function implements the Thomas algorithm, consisting of a
/// forward elimination phase followed by backward substitution.
/// The algorithm achieves linear time complexity by exploiting the
/// tridiagonal structure of the system.
///
/// The matrix is passed by value and modified internally.
template<dimension_t Dim>
constexpr StateVector<Dim> solveTridiagonal(tridiagonal_matrix_t<Dim> M, StateVector<Dim> psi) noexcept
{
	// --- FORWARD ELIMINATION ---
	for (dimension_t i = 1; i < Dim; ++i)
	{
		// Elimination multiplier
		const cplx_t w = M.a[i] / M.b[i - 1];

		// Update main diagonal
		M.b[i] = M.b[i] - w * M.c[i - 1];

		// Update right-hand side
		psi[i] = psi[i] - w * psi[i - 1];
	}

	// --- BACK SUBSTITUTION ---
	StateVector<Dim> Result{};

	Result[Dim - 1] = psi[Dim - 1] / M.b[Dim - 1];

	for (dimension_t i = Dim - 1; i-- > 0;)
	{
		Result[i] = (psi[i] - M.c[i] * Result[i + 1]) / M.b[i];
	}

	return Result;
}


/// @brief Callable object performing one Crank–Nicolson time step.
///
/// @details
/// The operator precomputes the Crank–Nicolson matrices A and B and
/// applies them to advance a quantum state vector by a single time step.
///
/// Usage:
///
///   CrankNicolsonTimeEvolutionOperator<Dim> evol(hamiltonian, dt);
///   psi = evol(psi);
template<dimension_t Dim>
class CrankNicolsonSolver
{
	// Precomputed matrices
	tridiagonal_matrix_t<Dim> m_A;
	tridiagonal_matrix_t<Dim> m_B;

public:
	/// @brief  Constructs the time evolution operator.
	/// @param  hamiltonian  Hamiltonian operator of the system.
	/// @param  dt           Time step size.
	///
	/// @details
	/// The constructor precomputes the Crank–Nicolson matrices A and B,
	/// which are reused for each time step.
	constexpr CrankNicolsonSolver(const Hamiltonian<Dim>& hamiltonian, double dt) noexcept
	{
		buildCrankNicolsonMatrices(hamiltonian, dt, m_A, m_B);
	}

	/// @brief  Advances the state vector by one time step.
	/// @param  psi     State vector at time step n.
	/// @return         State vector at time step n+1.
	///
	/// @details
	/// The function computes the right-hand side B · ψⁿ and then solves
	/// the linear system A · ψⁿ⁺¹ = RHS, resulting in unitary time evolution.
	constexpr StateVector<Dim>
		operator()(const StateVector<Dim>& psi) const noexcept
	{
		// RHS = B · ψⁿ
		auto rhs = multiplyTrigiagonal(m_B, psi);

		// Solve A · ψⁿ⁺¹ = RHS
		return solveTridiagonal(m_A, rhs);
	}
};
