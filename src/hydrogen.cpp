#include "visu/visu_oscilloscope.h"
#include "systems/particle_in_a_box.h"

using namespace Ket;

int main()
{
	constexpr OneDimensionalParticleBoxConfig<96> cfg(1.0, 1E-4);

	constexpr auto hydrogenOrbital = HydrogenOrbital<cfg.M>()(
		QuantumNumber::_2s(),
		0.05,
		cfg.dx
	);

	constexpr Ket::float_t mass = 1.0;
	constexpr auto hamiltonian = Hamiltonian<cfg.M>(mass, cfg.dx, SoftCoulombRadialPotential());

	OneDimensionalParticleBox<cfg.N> box(
		cfg,
		hamiltonian,
		hydrogenOrbital
	);
	
	Visu::VisuOscilloscope<cfg.M> visu;

	while (true)
	{
		auto p = box.evolve();
		visu.update(p);
	}
}
