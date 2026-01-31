#pragma once
#include <chrono>
#include <thread>
#include <iostream>
#include <array>
#include <tuple>
#include <algorithm>
#include <cmath>

#include "wavefunction/state_vector.h"
#include "constexprmath/constexpr_trigon.h"

#ifdef _WIN32
#include <windows.h>
#undef max
#endif

namespace KetCat::Visu
{
	/// @brief Enum to specify whether to use phase encoding in visualization
	enum class UsePhaseEncoding
	{
		YES,
		NO
	};

	/// @brief Enum to specify whether to clear the screen before updating visualization
	enum class ClearScreen
	{
		YES,
		NO
	};

	/// @brief Enum to enable visualization of real and imaginary parts
	enum class ShowComplexParts
	{
		YES,
		NO
	};

	/// @brief Check if an enum flag is enabled
	template<typename EnumType>
	constexpr inline bool enabled(EnumType e, EnumType yes = EnumType::YES) noexcept
	{
		return e == yes;
	}

	/// @brief Map phase angle arg(ψ) to ANSI color code
	inline const char* phaseToColor(float_t phase)
	{
		if (phase < -ConstexprMath::Pi * 0.5) return "\x1B[34m";  // dark blue
		if (phase < -ConstexprMath::Pi * 0.25) return "\x1B[94m"; // light blue
		if (phase < ConstexprMath::Pi * 0.25) return "\x1B[97m";  // white
		if (phase < ConstexprMath::Pi * 0.5) return "\x1B[91m";   // light red
		return "\x1B[31m";                                        // red
	}

	/// @brief Render a single oscilloscope line from (value, color) samples
	///
	/// @tparam Dim Dimension of the signal
	/// @param samples Array of tuples: (value, ANSI color)
	/// @param label Text label printed before the line
	///
	/// @details
	/// The values are visually normalized using:
	///   |vᵢ| / maxⱼ |vⱼ|
	///
	/// This normalization is purely for visualization and has
	/// no physical meaning.
	template<dimension_t Dim>
	inline void renderLine(
		const std::array<std::tuple<float_t, const char*>, Dim>& samples,
		const char* label)
	{
		static const char* bars[] = {
			"\xE2\x96\x81", "\xE2\x96\x82", "\xE2\x96\x83",
			"\xE2\x96\x84", "\xE2\x96\x85", "\xE2\x96\x86",
			"\xE2\x96\x87", "\xE2\x96\x88"
		};

		// Find maximum absolute value for normalization
		float_t maxVal = 0.0f;
		for (const auto& s : samples)
		{
			maxVal = std::max(maxVal, std::abs(std::get<0>(s)));
		}
		if (maxVal == 0.0f)
		{
			maxVal = 1.0f;
		}

		std::cout << label << " |";
		for (dimension_t i = 0; i < Dim; ++i)
		{
			const float_t value = std::get<0>(samples[i]);
			const char* color = std::get<1>(samples[i]);

			const float_t norm = std::abs(value) / maxVal;
			const std::size_t idx = static_cast<std::size_t>(norm * 7);

			std::cout << color << bars[idx] << "\x1B[0m";
		}
		std::cout << "|\n";
	}

	/// @brief Terminal-based oscilloscope visualization for 1D quantum states
	///
	/// @tparam Dim Dimension of the state vector
	///
	/// @details
	/// Displays:
	///  - Optional real part Re(ψ)  (yellow)
	///  - Optional imaginary part Im(ψ) (cyan)
	///  - Probability density |ψ|² (optionally phase-colored)
	template<dimension_t Dim>
	struct VisuOscilloscope
	{
		UsePhaseEncoding m_usePhaseEncoding;
		ClearScreen m_clearScreen;
		ShowComplexParts m_showComplex;

		/// @brief Construct a VisuOscilloscope with specified settings
		/// @param usePhaseEncoding Enable phase-based coloring of probability density
		/// @param clearScreen      Clear screen before rendering
		/// @param showComplex      Enable visualization of real and imaginary parts
		VisuOscilloscope(
			const UsePhaseEncoding usePhaseEncoding,
			const ClearScreen clearScreen,
			const ShowComplexParts showComplex)
			: m_usePhaseEncoding(usePhaseEncoding),
			  m_clearScreen(clearScreen),
			  m_showComplex(showComplex)
		{
#ifdef _WIN32
			// Enable UTF-8 output on Windows console
			SetConsoleOutputCP(CP_UTF8);
#endif
		}
		
		/// @brief Update the visualization with the current state vector
		/// @param s Current quantum state vector
		/// @param usePhaseEncoding Enable phase-based coloring of probability density
		/// @param cls Clear screen before rendering
		/// @param showComplex Enable visualization of real and imaginary parts
		void update(const StateVector<Dim>& s) const
		{
			using namespace std::chrono_literals;

			if (enabled(m_clearScreen))
			{
				std::cout << "\x1B[2J\x1B[H";
			}

			// --- Probability density |ψ|² ---
			std::array<std::tuple<float_t, const char*>, Dim> probLine{};

			for (dimension_t i = 0; i < Dim; ++i)
			{
				const float_t p = s[i].normSquared();

				const char* color = "\x1B[97m";
				if (enabled(m_usePhaseEncoding))
				{
					const float_t phase = std::atan2(s[i].im, s[i].re);
					color = phaseToColor(phase);
				}

				probLine[i] = { p, color };
			}

			renderLine<Dim>(probLine, "Proba: ");


			// --- Optional: Real and Imaginary parts ---
			if (enabled(m_showComplex))
			{
				std::array<std::tuple<float_t, const char*>, Dim> reLine{};
				std::array<std::tuple<float_t, const char*>, Dim> imLine{};

				for (dimension_t i = 0; i < Dim; ++i)
				{
					reLine[i] = { s[i].re, "\x1B[33m" }; // Yellow
					imLine[i] = { s[i].im, "\x1B[36m" }; // Cyan
				}

				renderLine<Dim>(reLine, "Real:  ");
				renderLine<Dim>(imLine, "Imag:  ");
			}

			// Small delay to allow visualization update
			std::this_thread::sleep_for(100ms);
		}
	};
}
