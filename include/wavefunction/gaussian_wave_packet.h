#pragma once
#include "core_types.h"
#include "state_vector.h"


/// @brief Functor to generate a Gaussian wave packet state vector.
/// @tparam Dim The dimension of the state vector to generate.
///	@details This functor generates a Gaussian wave packet characterized by its
///			center position (x0), central wave number (k0), and standard deviation (sigma).
///			Returns a state vector of dimension Dim with complex amplitudes.
template<dimension_t Dim>
struct GaussianWavePacket
{
	constexpr StateVector<Dim> operator()(double x0, double k0, double sigma, double dx) const noexcept
	{
		StateVector<Dim> wavePacket = {};

		for (dimension_t n = 0; n < Dim; ++n)
		{
			// Position corresponding to index n
			const double x = (n + 1) * dx;

			// Gaussian envelope calculation 
			// exp(-((x - x0)^2) / (4 * sigma^2))
			const double exponent = -((x - x0) * (x - x0)) / (4.0 * sigma * sigma);
			const double envelope = ConstexprMath::exp_taylor<20>(exponent);

			// Plane wave component calculation: cos(k0 * x) + i * sin(k0 * x)
			const double realPart = ConstexprMath::cos(k0 * x);
			const double imagPart = ConstexprMath::sin(k0 * x);

			// Combine envelope and plane wave to form the complex amplitude
			wavePacket[n] = cplx_t(envelope * realPart, envelope * imagPart);
		}

		wavePacket.normalize();

		return wavePacket;
	}
};
