#pragma once
#include "constexprmath/constexpr_core_functions.h"


/// @brief 1D soft-Coulomb atompotenciál
/// V(x) = -Z / sqrt(x^2 + a^2)
class SoftCoulombPotential
{
	double m_Z;     // atomszám
	double m_a;     // softening parameter
	double m_x0;    // atom középpontja

public:
	constexpr SoftCoulombPotential(double Z = 1.0, double a = 1e-10, double x0 = 0.0) noexcept
		: m_Z(Z), m_a(a), m_x0(x0)
	{
	}

	constexpr double operator()(double x) const noexcept
	{
		const double dx = x - m_x0;
		return -m_Z / ConstexprMath::sqrt(dx * dx + m_a * m_a);
	}
};
