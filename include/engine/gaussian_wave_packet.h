#pragma once
#include "types.h"

/// @brief Compute the exponential function using Taylor series expansion.
/// @param x The exponent value.
/// @param N The number of terms in the Taylor series (default is 20).
template <unsigned int Terms, std::floating_point FloatType>
constexpr FloatType exp_taylor(FloatType x)
{
	FloatType sum = 1.0;
	FloatType term = 1.0;
	for (unsigned n = 1; n <= Terms; ++n)
	{
		term *= x / n;
		sum += term;
	}
	return sum;
}

template <unsigned int Iterations, std::floating_point FloatType>
constexpr FloatType sqrt_newton_raphson(FloatType x)
{
	FloatType Guess = 1.0;
	for (unsigned int i = 0; i < Iterations; ++i)
	{
		Guess = 0.5 * (Guess + x / Guess);
	}
	return Guess;
}

/// @brief Functor to generate a Gaussian wave packet state vector.
/// @tparam Dim The dimension of the state vector to generate.
///	@details This functor generates a Gaussian wave packet characterized by its
///			center position (x0), central wave number (k0), and standard deviation (sigma).
///			Returns a state vector of dimension Dim with complex amplitudes.
template<dimension_t Dim>
struct GaussianWavePacket
{
	constexpr state_vector_t<Dim> operator()(double x0, double k0, double sigma, double dx) const noexcept
	{
		state_vector_t<Dim> wavePacket = {};
		// Compute the wave packet values and normalization factor
		double normSquared = 0.0;

		for (dimension_t n = 0; n < Dim; ++n)
		{
			// Position corresponding to index n
			const double x = (n + 1) * dx;

			// Gaussian envelope calculation 
			// exp(-((x - x0)^2) / (4 * sigma^2))
			const double exponent = -((x - x0) * (x - x0)) / (4.0 * sigma * sigma);
			const double envelope = exp_taylor<20>(exponent);

			// Plane wave component calculation: cos(k0 * x) + i * sin(k0 * x)
			const double realPart = ConstexprMath::cos(k0 * x);
			const double imagPart = ConstexprMath::sin(k0 * x);

			// Combine envelope and plane wave to form the complex amplitude
			wavePacket[n] = cplx_t(envelope * realPart, envelope * imagPart);
			normSquared += wavePacket[n].normSquared();
		}

		const double norm = sqrt_newton_raphson<10>(normSquared);
		for (cplx_t& c : wavePacket)
			c = c / norm;

		return wavePacket;
	}
};
