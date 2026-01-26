#pragma once
#include "core_types.h"

namespace Ket
{
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
			float_t normSquared = 0.0;

			for (cplx_t& c : m_StateVector)
			{
				normSquared += c.normSquared();
			}

			const float_t norm = ConstexprMath::sqrt(normSquared);
			for (cplx_t& c : m_StateVector)
			{
				c = c / norm;
			}
		}
		/*
		/// @brief Compute the Kronecker product of this state vector with another.
		/// @tparam OtherDim  Dimension of the other state vector.
		/// @param other  The other state vector.
		/// @return The Kronecker product state vector.
		/// @note The resulting state vector has dimension HilbertDim * OtherDim
		template<dimension_t OtherDim>
		constexpr StateVector<Hilbertdim * OtherDim> kroneckerProduct(const StateVector<OtherDim>& other) const noexcept
		{
			StateVector<HilbertDim* OtherDim> Result;
			for (dimension_t i = 0; i < HilbertDim; ++i)
			{
				for (dimension_t j = 0; j < OtherDim; ++j)
				{
					Result[i * OtherDim + j] = m_StateVector[i] * other.m_StateVector[j];
				}
			}
			return Result;
		}
		*/
		/// @brief Multiply this state vector by a matrix.
		/// @param mat  The matrix to multiply with (HilbertDim x HilbertDim).
		/// @return The resulting state vector.
		constexpr StateVector<HilbertDim> matMul(const matrix_t<HilbertDim>& mat) const noexcept
		{
			StateVector<HilbertDim> Result;
			for (dimension_t i = 0; i < HilbertDim; ++i)
			{
				cplx_t Sum = cplx_t::zero();
				for (dimension_t j = 0; j < HilbertDim; ++j)
				{
					Sum = Sum + mat[i][j] * m_StateVector[j];
				}
				Result[i] = Sum;
			}
			return Result;
		}
	};
}