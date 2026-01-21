#pragma once
#include "core_types.h"
#include "state_vector.h"


template<dimension_t Dim>
struct EigenState
{
	// n = 1,2,3,... (kvantumszám)
	constexpr StateVector<Dim>
		operator()(unsigned int n, double dx, double L) const noexcept
	{
		StateVector<Dim> psi{};

		for (dimension_t i = 0; i < Dim; ++i)
		{
			// fizikai pozíció (Dirichlet peremek között)
			const double x = (i + 1) * dx;

			// sajátállapot alakja
			const double value =
				ConstexprMath::sin(
					n * ConstexprMath::Pi * x / L
				);

			psi[i] = cplx_t(value, 0.0);
		}

		psi.normalize();

		return psi;
	}
};

