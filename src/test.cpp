#include "visu/visu_oscilloscope.h"
#include "systems/particle_in_a_box.h"

using namespace Ket;

int main()
{
	constexpr OneDimensionalParticleBoxConfig<96> cfg(1.0, 1E-4);

    constexpr Ket::float_t x0 = 0.01;
	constexpr Ket::float_t sigma = 0.1;
	constexpr Ket::float_t k0 = ConstexprMath::Pi * 10;

	constexpr auto gaussianPacket = GaussianWavePacket<cfg.M>()(x0, k0, sigma, cfg.dx);

	constexpr PotentialBarrier potential{
		0.45, 0.55, // potential wall in the middle
		3000 // Joule
	};

	constexpr Ket::float_t mass = 1.0; // Electron mass in kg
	constexpr auto hamiltonian = Hamiltonian<cfg.M>(mass, cfg.dx, potential);

	OneDimensionalParticleBox<cfg.N> box(
		cfg,
		hamiltonian,
		gaussianPacket
	);
	
	
	Visu::VisuOscilloscope<cfg.M> visu;

	while (true)
	{
		const auto p = box.evolve();
		visu.update(p);
	}

}
