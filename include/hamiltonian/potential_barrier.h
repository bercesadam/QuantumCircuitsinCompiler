#pragma once

namespace Ket
{
	/// @brief Represents a potential function
	/// @details This class defines a potential function that can be used
	///          in the Hamiltonian of a quantum system. The potential is
	///          defined over a specific range and has a constant value within that range.
	///	Example usage:
	///	auto potentialWall = Potential(0.0, 5.0, 200.0);
	/// auto potentialWell = Potential(0.0, 5.0, -100.0);
	class PotentialBarrier
	{
		float_t m_Start;
		float_t m_End;
		float_t m_V0;

	public:
		constexpr PotentialBarrier() noexcept
			: m_Start(0.0), m_End(0.0), m_V0(0.0)
		{
		}

		/// @brief Constructs a potential barrier
		///
		/// @param start   Start position
		/// @param end     End position
		/// @param v0      The constant potential energy
		constexpr PotentialBarrier(float_t start, float_t end, float_t v0) noexcept
			: m_Start(start), m_End(end), m_V0(v0)
		{
		}

		/// @brief Evaluates the potential at a given spatial position.
		///
		/// @param x Spatial coordinate
		/// @return Potential energy V(x)
		constexpr float_t operator()(float_t position) const noexcept
		{
			if (position >= m_Start && position <= m_End)
			{
				return m_V0;
			}
			else
			{
				return 0.0;
			}
		}
	};

	/// @brief A zero potential barrier instance
	/// @details This constant represents a potential barrier with zero potential energy
	///          everywhere. It can be used in scenarios where no potential is desired.
	constexpr auto ZeroPotential = PotentialBarrier();
}