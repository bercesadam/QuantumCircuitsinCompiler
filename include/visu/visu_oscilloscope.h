#pragma once
#include <chrono>
#include <thread>
#include <iostream>

#include "visu_iface.h"
#include "wavefunction/state_vector.h"

#ifdef _WIN32
#include <windows.h>
#undef max
#endif

namespace Ket::Visu
{
	inline const char* phaseToColor(float_t phase)
	{
		constexpr float_t PI = 3.14159265358979323846f;

		if (phase < -PI * 0.5f) return "\x1B[34m"; // sötétkék
		if (phase < -PI * 0.25f) return "\x1B[94m"; // világoskék
		if (phase < PI * 0.25f) return "\x1B[97m"; // fehér
		if (phase < PI * 0.5f) return "\x1B[91m"; // halványpiros
		return "\x1B[31m";                        // piros
	}

	/// @brief Simple terminal-based oscilloscope visualization implementation
	/// @tparam Dim Dimension of the state vector
	/// @detail Uses Unicode block characters to represent probability densities.
	template<dimension_t Dim>
	struct VisuOscilloscope
	{
		/// @brief Update the visualization with the current state vector
		/// @param s Current state vector
		void update(const StateVector<Dim>& s)
		{
			using namespace std::chrono_literals;

			std::cout << "\x1B[2J\x1B[H";

			static const char* bars[] = {
			"\xE2\x96\x81", "\xE2\x96\x82", "\xE2\x96\x83",
			"\xE2\x96\x84", "\xE2\x96\x85", "\xE2\x96\x86",
			"\xE2\x96\x87", "\xE2\x96\x88"
			};

			const auto Probabilities = s.getProbabilities();

			float_t MaxProba = 0.0;
			for (float_t p : Probabilities)
			{
				MaxProba = std::max(MaxProba, p);
			}

			std::cout << "Probability density: |";
			for (dimension_t i = 0; i < Probabilities.size(); ++i)
			{
				float_t p = Probabilities[i];
				p /= MaxProba;
				std::size_t idx = static_cast<std::size_t>(p * 7);
				//std::cout << bars[idx];
				
				float_t phase = std::atan2(s[i].im, s[i].re);
				const char* color = phaseToColor(phase);

				std::cout << color << bars[idx] << "\x1B[0m";
			}
			std::cout << "|\n";

			std::this_thread::sleep_for(100ms);
		}

		VisuOscilloscope()
		{
#ifdef _WIN32
			// Enable UTF-8 output on Windows console
			SetConsoleOutputCP(CP_UTF8);
#endif
		}
	};
}