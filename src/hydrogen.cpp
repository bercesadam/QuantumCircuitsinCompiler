#include <tuple>
#include <string>
#include <functional>

#include "visu/visu_oscilloscope.h"
#include "systems/particle_in_a_box.h"

using namespace KetCat;

int main()
{
	constexpr OneDimensionalParticleBoxConfig<96> cfg(1.0, 1E-4);

	// Hydrogen orbital constructor with fixed parameters
	constexpr auto hydrogenCtor = std::bind(HydrogenOrbital<cfg.M>(), std::placeholders::_1, 0.05, cfg.dx);

	// List of hydrogen orbitals to simulate: (StateVector, Name, l)
	constexpr std::array<std::tuple<StateVector<94>, std::string, unsigned int>, 6> hydrogenOrbitals =
	{
		 std::make_tuple(hydrogenCtor(QuantumNumber::_1s()), "1s", 0),
		 std::make_tuple(hydrogenCtor(QuantumNumber::_2s()), "2s", 0),
		 std::make_tuple(hydrogenCtor(QuantumNumber::_2p()), "2p", 1),
		 std::make_tuple(hydrogenCtor(QuantumNumber::_3s()), "3s", 0),
		 std::make_tuple(hydrogenCtor(QuantumNumber::_3p()), "3p", 1),
		 std::make_tuple(hydrogenCtor(QuantumNumber::_3d()), "3d", 2)
	};

	constexpr KetCat::float_t mass = 1.0;

	for (const auto& orbital : hydrogenOrbitals)
	{
		std::cout << "Hydrogen orbital: " << std::get<1>(orbital) << "\n";

		auto hamiltonian = Hamiltonian<cfg.M>(mass, cfg.dx, SoftCoulombRadialPotential(
			1.0,
			2e-2,
			std::get<2>(orbital),
			hBar,
			mass
		));

		auto box = OneDimensionalParticleBox<cfg.N>(
			cfg,
			hamiltonian,
			std::get<0>(orbital)
		);
		
		Visu::VisuOscilloscope<cfg.M>(
			Visu::UsePhaseEncoding::NO,
			Visu::ClearScreen::NO,
			Visu::ShowComplexParts::NO
		).update(box.evolve());
	}

}
