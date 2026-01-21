#pragma once
#include "core_types.h"
#include "state_vector.h"


class QuantumNumber
{
	unsigned int m_n;
	unsigned int m_l;

	constexpr QuantumNumber(unsigned int n, unsigned int l)
		: m_n(n), m_l(l)
	{}

public:
	constexpr unsigned int n() { return m_n; }
	constexpr unsigned int l() { return m_l; }

	static constexpr QuantumNumber _2s() { return QuantumNumber{ 2, 0 }; }
	static constexpr QuantumNumber _2p() { return QuantumNumber{ 2, 1 }; }
	static constexpr QuantumNumber _3s() { return QuantumNumber{ 3, 0 }; }
	static constexpr QuantumNumber _3p() { return QuantumNumber{ 3, 1 }; }
	static constexpr QuantumNumber _3d() { return QuantumNumber{ 3, 2 }; }
	static constexpr QuantumNumber _4s() { return QuantumNumber{ 4, 0 }; }
	static constexpr QuantumNumber _4p() { return QuantumNumber{ 4, 1 }; }
	static constexpr QuantumNumber _4d() { return QuantumNumber{ 4, 2 }; }
	static constexpr QuantumNumber _4f() { return QuantumNumber{ 4, 3 }; }
};


template<dimension_t Dim>
struct HydrogenOrbitals
{
	/// @param n   főkvantumszám (n >= 1)
	/// @param l   "impulzusmomentum" analóg (0 <= l <= n-1)
	/// @param a0  effektív Bohr-sugár
	/// @param dx  rácslépés
	/// @param x0  atom közepe
	constexpr StateVector<Dim>
	operator()(QuantumNumber q,
	           double a0, double dx,
	           double x0 = 0.0) const noexcept
	{
		StateVector<Dim> psi{};
		const unsigned int n = q.n();
		const unsigned int l = q.l();

		for (dimension_t i = 0; i < Dim; ++i)
		{
			const double x = i * dx;
			const double r = ConstexprMath::abs(x - x0);

			// --- exponenciális lecsengés ---
			const double ExpPart =
				ConstexprMath::exp_taylor<30>(-r / (n * a0));

			// --- "centrifugális" prefaktor ---
			double PolyPart = 1.0;
			for (unsigned k = 0; k < l; ++k)
				PolyPart *= (r / (n * a0));

			// --- csomópontok (Laguerre-analóg, egyszerűsített) ---
			double NodePart = 1.0;
			for (unsigned k = 1; k <= n - l - 1; ++k)
				NodePart *= (1.0 - r / (k * n * a0));

			double Value = ExpPart * PolyPart * NodePart;

			// --- paritás ---
			if (l % 2 == 1 && x < x0)
				Value = -Value;

			psi[i] = cplx_t::fromReal(Value);
		}

		psi.normalize();

		return psi;
	}
};
