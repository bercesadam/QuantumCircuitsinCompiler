#pragma once
#include "types.h"
#include "hilbert.h"

/// @brief Represents a quantum state vector in a Hilbert space of given dimension.
/// @tparam HilbertDim  Dimension of the Hilbert space (number of basis states).
template <dimension_t HilbertDim>
struct StateVector
{
	state_vector_t<HilbertDim> m_StateVector;

public:
	/// @brief Get the probabilities of measuring selected basis states.
	/// @tparam SelectedStates  Indices of the basis states to get probabilities for.
	/// @return A probability vector with the probabilities of the selected states.
	template<dimension_t... SelectedStates>
	constexpr probability_vector_t getProbabilities() noexcept
	{
		constexpr dimension_t NumSelected = sizeof...(SelectedStates);
		constexpr subspace_indices_t<NumSelected> SelectedStateIndices = { SelectedStates... };
		probability_vector_t<NumSelected> Probabilities = {};
		for (dimension_t i : SelectedStateIndices)
		{
			const dimension_t StateIndex = SelectedStatesArray[i];
			// Compute probability as the squared norm of the amplitude
			Probabilities[i] = m_StateVector[StateIndex].normSquared();
		}
		return Probabilities;
	}

	void normalize() noexcept
	{
		double norm = 0.0;
		for (const auto& amplitude : m_StateVector)
		{
			norm += amplitude.normSquared();
		}
		if (norm > 0.0)
		{
			const double invSqrtNorm = 1.0 / std::sqrt(norm);
			for (auto& amplitude : m_StateVector)
			{
				amplitude = amplitude * cplx_t::fromReal(invSqrtNorm);
			}
		}
	}
};