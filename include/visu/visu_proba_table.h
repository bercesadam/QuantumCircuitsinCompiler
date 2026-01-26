#pragma once

#include "visu_iface.h"

namespace Ket::Visu
{
	/// @brief Simple terminal-based probability table visualization implementation
	/// @tparam Dim Dimension of the state vector
	/// @detail Displays probability densities as a table of numerical values.
	template<dimension_t Dim>
	struct VisuProbaTable : public VisuIface<Dim>
	{
		/// @brief Update the visualization with the current state vector
		/// @param s Current state vector
		void update(const StateVector<Dim>& s) override
		{
			const auto Probabilities = s.getProbabilities();
			std::cout << "Probability densities:\n";
			for (std::size_t i = 0; i < Dim; ++i)
			{
				std::cout << "State |" << i << ">: " << Probabilities[i] << "\n";
			}
			std::cout << std::endl;
		}
	};
} // namespace Ket::Visu