#pragma once


/// @brief Represents a potential function
/// @details This class defines a potential function that can be used
///          in the Hamiltonian of a quantum system. The potential is
///          defined over a specific range and has a constant value within that range.
///	Example usage:
///	auto potentialWall = Potential(0.0, 5.0, 200.0);
/// auto potentialWell = Potential(0.0, 5.0, -100.0);
class PotentialBarrier
{
	double m_Start;
	double m_End;
	double m_V0;

public:
	constexpr PotentialBarrier() noexcept
		: m_Start(0.0), m_End(0.0), m_V0(0.0)
	{
	}

	constexpr PotentialBarrier(double start, double end, double v0) noexcept
		: m_Start(start), m_End(end), m_V0(v0)
	{
	}

	constexpr double operator()(double position) const noexcept
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