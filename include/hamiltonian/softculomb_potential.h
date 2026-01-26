#pragma once
#include "constexprmath/constexpr_core_functions.h"

namespace Ket
{
	/// @brief One-dimensional soft-Coulomb atomic potential.
	///
	/// This class implements a regularized (soft-core) Coulomb potential
	/// commonly used in 1D quantum simulations to avoid the singularity
	/// at the origin.
	///
	/// Mathematical form:
	/// V(x) = -Z / sqrt( (x - x0)^2 + a^2 )
	///
	/// Physical interpretation:
	/// - Z  : Effective nuclear charge
	/// - a  : Softening parameter that removes the 1/r divergence at x = x0
	/// - x0 : Position of the atomic center
	class SoftCoulombPotential
	{
		/// Effective nuclear charge (atomic number)
		float_t m_Z;
		/// Softening parameter controlling short-range regularization
		float_t m_a;
		/// Position of the atomic center
		float_t m_x0;

	public:
		/// @brief Constructs a soft-Coulomb potential.
		///
		/// @param Z   Effective nuclear charge
		/// @param a   Softening parameter (prevents singularity at the origin)
		/// @param x0  Position of the atomic center
		constexpr SoftCoulombPotential(float_t Z = 1.0,
			float_t a = 1e-10,
			float_t x0 = 0.0) noexcept
			: m_Z(Z), m_a(a), m_x0(x0)
		{
		}

		/// @brief Evaluates the potential at a given spatial position.
		///
		/// @param  x Spatial coordinate
		/// @return Potential energy V(x)
		constexpr float_t operator()(float_t x) const noexcept
		{
			// Displacement from the atomic center
			const float_t dx = x - m_x0;

			// Soft-Coulomb potential evaluation
			return -m_Z / ConstexprMath::sqrt(dx * dx + m_a * m_a);
		}
	};
}
