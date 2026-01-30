#pragma once
#include "constexprmath/constexpr_core_functions.h"

namespace KetCat
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

	/// @brief Radial soft-Coulomb potential for spherically symmetric systems.
	/// @details
	/// This class implements a radial soft-Coulomb potential commonly used
	/// in quantum simulations of atomic systems. The potential includes both
	/// the Coulomb attraction term and the centrifugal barrier term for
	/// angular momentum states.
	class SoftCoulombRadialPotential
	{
		double m_Zeff;     // Effective nuclear charge Z_eff
		double m_a;        // Softening parameter a > 0
		unsigned m_l;      // Orbital quantum number ℓ ≥ 0
		double m_hbar;     // ħ
		double m_mu;       // Reduced mass μ

	public:
		constexpr SoftCoulombRadialPotential(double Zeff = 1.0,
			double a = 2e-2,
			unsigned l = 0,
			double hbar = 1.0,
			double mu = 1.0) noexcept
			: m_Zeff(Zeff), m_a(a), m_l(l), m_hbar(hbar), m_mu(mu)
		{
		}

		/// @brief  Evaluate V(r) at radius r.
		/// @param  r   Radius r ≥ 0.
		/// @return     Potential value V(r).
		///
		/// @details
		/// Computation:
		///   r² + a² → r2s
		///   V(r) = −Z_eff / √(r2s)  +  ℓ(ℓ+1)·ħ² / (2μ·r2s)
		constexpr double operator()(double r) const noexcept
		{
			const double SquareRadiusWithSoftening = r * r + m_a * m_a;  // r² + a²
			const double CoulombAttraction = -m_Zeff / ConstexprMath::sqrt(SquareRadiusWithSoftening);
			const double CentrifugalBarrier = (static_cast<double>(m_l) * (m_l + 1) * m_hbar * m_hbar) / (2.0 * m_mu * SquareRadiusWithSoftening);
			return CoulombAttraction + CentrifugalBarrier;
		}
	};
}
