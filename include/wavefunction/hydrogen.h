#pragma once
#include "core_types.h"
#include "state_vector.h"

namespace Ket
{
	/// @brief Represents a pair of quantum numbers (n, l) using the spectroscopic notation:
	/// - n : Principal quantum number (n >= 1), controls energy and radial extent
	/// - l : Orbital angular momentum quantum number (0 <= l <= n-1),
	///       controls angular structure and parity
	class QuantumNumber
	{
		/// Principal quantum number
		unsigned int m_n;
		/// Orbital angular momentum quantum number
		unsigned int m_l;

		/// @brief Private constructor enforcing valid factory usage.
		/// @param n Principal quantum number
		/// @param l Orbital angular momentum quantum number
		constexpr QuantumNumber(unsigned int n, unsigned int l)
			: m_n(n), m_l(l)
		{
		}

	public:
		/// @brief Returns the principal quantum number.
		constexpr unsigned int n() { return m_n; }

		/// @brief Returns the orbital angular momentum quantum number.
		constexpr unsigned int l() { return m_l; }

		/// @name Hydrogen orbital presets
		/// @{
		static constexpr QuantumNumber _1s() { return QuantumNumber{ 1, 0 }; }
		static constexpr QuantumNumber _2s() { return QuantumNumber{ 2, 0 }; }
		static constexpr QuantumNumber _2p() { return QuantumNumber{ 2, 1 }; }
		static constexpr QuantumNumber _3s() { return QuantumNumber{ 3, 0 }; }
		static constexpr QuantumNumber _3p() { return QuantumNumber{ 3, 1 }; }
		static constexpr QuantumNumber _3d() { return QuantumNumber{ 3, 2 }; }
		static constexpr QuantumNumber _4s() { return QuantumNumber{ 4, 0 }; }
		static constexpr QuantumNumber _4p() { return QuantumNumber{ 4, 1 }; }
		static constexpr QuantumNumber _4d() { return QuantumNumber{ 4, 2 }; }
		static constexpr QuantumNumber _4f() { return QuantumNumber{ 4, 3 }; }
		/// @}
	};


	/// @brief Generator functor for simplified 1D hydrogen-like orbital wavefunctions.
	///
	/// This class constructs a state vector that mimics the radial structure of
	/// hydrogen atom orbitals, projected onto a one-dimensional spatial grid.
	///
	/// The implementation is intentionally simplified:
	/// - Angular dependence is encoded via parity and polynomial factors
	/// - Radial nodes are approximated using a Laguerre-like product
	/// - The result is suitable for visualization, experimentation,
	///   and constexpr evaluation rather than exact spectroscopy
	///
	/// @tparam Dim Size of the discrete spatial grid
	template<dimension_t Dim>
	struct HydrogenOrbitals
	{
		/// @brief Generates a hydrogen-like orbital wavefunction.
		///
		/// @param q   Quantum number pair (n, l)
		/// @param a0  Effective Bohr radius (controls spatial scale)
		/// @param dx  Spatial discretization step
		/// @param x0  Position of the atomic center
		///
		/// @return Normalized quantum state vector representing the orbital
		constexpr StateVector<Dim>
			operator()(QuantumNumber q,
				float_t a0, float_t dx,
				float_t x0 = 0.0) const noexcept
		{
			StateVector<Dim> psi{};
			const unsigned int n = q.n();
			const unsigned int l = q.l();

			for (dimension_t i = 0; i < Dim; ++i)
			{
				// Spatial coordinate corresponding to grid index i
				const float_t x = i * dx;

				// Radial distance from the atomic center
				const float_t r = ConstexprMath::abs(x - x0);

				// --- Exponential decay term ---
				//
				// Models the asymptotic behavior of hydrogenic orbitals:
				// exp(-r / (n * a0))
				const float_t ExpPart =
					ConstexprMath::exp_taylor<30>(-r / (n * a0));

				// --- Centrifugal-like polynomial prefactor ---
				//
				// Represents the r^l behavior near the origin caused
				// by angular momentum
				float_t PolyPart = 1.0;
				for (unsigned k = 0; k < l; ++k)
					PolyPart *= (r / (n * a0));

				// --- Radial nodes (simplified Laguerre analogue) ---
				//
				// Introduces (n - l - 1) radial nodes using a product form.
				// This is not an exact Laguerre polynomial but preserves
				// qualitative nodal structure.
				float_t NodePart = 1.0;
				for (unsigned k = 1; k <= n - l - 1; ++k)
					NodePart *= (1.0 - r / (k * n * a0));

				// Combined radial wavefunction value
				float_t Value = ExpPart * PolyPart * NodePart;

				// --- Parity handling ---
				//
				// Odd-l orbitals have antisymmetric parity in 1D projection
				if (l % 2 == 1 && x < x0)
					Value = -Value;

				// Store as a purely real complex amplitude
				psi[i] = cplx_t::fromReal(Value);
			}

			// Normalize the wavefunction to unit probability
			psi.normalize();

			return psi;
		}
	};
}
