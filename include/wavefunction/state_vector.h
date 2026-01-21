#pragma once
#include "core_types.h"


/// @brief Represents a quantum state vector in a Hilbert space of given dimension.
/// @tparam HilbertDim  Dimension of the Hilbert space (number of basis states).
template <dimension_t HilbertDim>
struct StateVector
{
	state_vector_t<HilbertDim> m_StateVector;

public:
	constexpr state_vector_t<HilbertDim> getVector() const noexcept
	{
		return m_StateVector;
	}

	constexpr cplx_t& operator[](dimension_t index) noexcept
	{
		return m_StateVector.at(index);
	}

	constexpr const cplx_t& operator[](dimension_t index) const noexcept
	{
		return m_StateVector.at(index);
	}

	/// @brief Get the probabilities of measuring selected basis states.
	/// @tparam SelectedStates  Indices of the basis states to get probabilities for.
	/// @return A probability vector with the probabilities of the selected states.
	template<dimension_t... SelectedStates>
	constexpr probability_vector_t<HilbertDim> getProbabilities() noexcept
	{
		constexpr dimension_t NumSelected = sizeof...(SelectedStates);
		constexpr subspace_indices_t<NumSelected> SelectedStateIndices = { SelectedStates... };
		probability_vector_t<NumSelected> Probabilities = {};
		for (dimension_t i : SelectedStateIndices)
		{
			const dimension_t StateIndex = SelectedStateIndices[i];
			// Compute probability as the squared norm of the amplitude
			Probabilities[i] = m_StateVector[StateIndex].normSquared();
		}
		return Probabilities;
	}


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