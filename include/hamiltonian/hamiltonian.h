#pragma once
#include "core_types.h"
#include "softculomb_potential.h"


///@brief Constants used for contructing a Hamiltonian class instance
struct HamiltonConstants
{
	// Reduced Planck's constant, in J·s
	// but usually set to 1.0 for simplicity in quantum simulations
	double hBar = 1.0545718e-34;
	// Weight of the particle (in kg), default is electron mass
	// but can be set to 1.0 for simplicity
	double m = 9.10938356e-31;
	// Spatial discretization step (in meters), default is 
	// a small value like 1e-10 for atomic scale simulations
	double dx = 1e-10;
};

/// @brief Represents the Hamiltonian operator in 1D discreitized space
/// @tparam Dim The dimension of the Hamiltonian matrix
/// @details This class constructs the Hamiltonian matrix for a quantum system
/// 		based on the provided constants and potential function.
///			Realizes the following equation: 
/// ///			H = - (ħ² / 2m·Δx²) · (d²/dx²) + V(x)
template<dimension_t Dim>
class Hamiltonian
{
	HamiltonConstants m_constants;
	matrix_t<Dim> m_hamiltonianMatrix;

public:
	constexpr matrix_t<Dim> getMatrix() const noexcept
	{
		return m_hamiltonianMatrix;
	}

	constexpr HamiltonConstants getConstants() const noexcept
	{
		return m_constants;
	}

	// Design limitation: currently only one potential can be used to construct the Hamiltonian
	constexpr Hamiltonian(const HamiltonConstants& constants, const SoftCoulombPotential& potential) noexcept
		: m_constants(constants)
	{
		// Initialize Hamiltonian matrix with zeros
		m_hamiltonianMatrix = {};

		// α = ħ² / (2m·Δx²)
		const double Alpha = constants.hBar * constants.hBar / (2.0 * constants.m * constants.dx * constants.dx);
		for (dimension_t i = 0; i < Dim; ++i)
		{
			for (dimension_t j = 0; j < Dim; ++j)
			{
				if (i == j)
				{
					// Calculating position for the Potential callable: i * Δx					
					double Position = i * constants.dx;

					// Diagonal elements: Kinetic + Potential energy
					// Kinetic part: -2α
					// Potential part: V(position)
					// Total: -2α + V(position)
					m_hamiltonianMatrix[i][j] = cplx_t::fromReal(2.0 * Alpha + potential(Position));
				}
				else if (j == i + 1 || j == i - 1)
				{
					m_hamiltonianMatrix[i][j] = cplx_t::fromReal(-Alpha);
				}
				else
				{
					m_hamiltonianMatrix[i][j] = cplx_t::fromReal(0.0);
				}
			}
		}
	}

	// Delegating constructor to allow default zero potential as convenience
	constexpr Hamiltonian(const HamiltonConstants& constants) noexcept
		: Hamiltonian(constants, SoftCoulombPotential())
	{
	}
};
