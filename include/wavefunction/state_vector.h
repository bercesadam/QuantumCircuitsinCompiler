#pragma once
#include "core_types.h"


/// @brief Represents a quantum state vector in a Hilbert space of given dimension.
/// @tparam HilbertDim  Dimension of the Hilbert space (number of basis states).
template <dimension_t HilbertDim>
struct StateVector
{
	/// Underlying state vector array
	state_vector_t<HilbertDim> m_StateVector;

public:
	/// @brief Indexing operator
	/// @return Reference to a complex number at the given state index
	constexpr cplx_t& operator[](dimension_t index) noexcept
	{
		return m_StateVector.at(index);
	}

	/// @brief Indexing operator (const)
	/// @return Const reference to a complex number at the given state index
	constexpr const cplx_t& operator[](dimension_t index) const noexcept
	{
		return m_StateVector.at(index);
	}

	/// @brief Get the probabilities of measuring the selected basis states.
	constexpr probability_vector_t<HilbertDim> getProbabilities() const noexcept
	{
		probability_vector_t<HilbertDim> Probabilities;
		for (int i = 0; i < HilbertDim; ++i)
		{
			Probabilities[i] = m_StateVector[i].normSquared();
		}
		return Probabilities;
	}

	/// @brief Normalizes the state vector (useful after initialization 
	/// with wavefunction functors to keep |ψ|² = 1).
	constexpr void normalize() noexcept
	{
		double normSquared = 0.0;

		for (cplx_t& c: m_StateVector)
		{
			normSquared += c.normSquared();
		}

		const double norm = ConstexprMath::sqrt(normSquared);
		for (cplx_t& c : m_StateVector)
		{
			c = c / norm;
		}
	}
};