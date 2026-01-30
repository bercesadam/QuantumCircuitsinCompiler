#include <tuple>
#include <string>
#include <functional>

#include "visu/visu_oscilloscope.h"
#include "systems/particle_in_a_box.h"

using namespace KetCat;

int main()
{
	constexpr OneDimensionalParticleBoxConfig<96> cfg(1.0, 1E-4);

	constexpr auto hydrogenCtor = std::bind(HydrogenOrbital<cfg.M>(), std::placeholders::_1, 0.05, cfg.dx);
	constexpr std::array<std::tuple<StateVector<94>, std::string>, 6> hydrogenOrbitals =
	{
		 std::make_tuple(hydrogenCtor(QuantumNumber::_1s()), "1s"),
		 std::make_tuple(hydrogenCtor(QuantumNumber::_2s()), "2s"),
		 std::make_tuple(hydrogenCtor(QuantumNumber::_2p()), "2p"),
		 std::make_tuple(hydrogenCtor(QuantumNumber::_3s()), "3s"),
		 std::make_tuple(hydrogenCtor(QuantumNumber::_3p()), "3p"),
		 std::make_tuple(hydrogenCtor(QuantumNumber::_3d()), "3d")
	};

	constexpr KetCat::float_t mass = 1.0;
	constexpr auto hamiltonian = Hamiltonian<cfg.M>(mass, cfg.dx, SoftCoulombRadialPotential());

	for (const auto& orbital : hydrogenOrbitals)
	{
		std::cout << "Hydrogen orbital: " << std::get<1>(orbital) << "\n";

		auto box = OneDimensionalParticleBox<cfg.N>(
			cfg,
			hamiltonian,
			std::get<0>(orbital)
		);
		
		Visu::VisuOscilloscope<cfg.M>().update(box.evolve(), false);
	}

}
