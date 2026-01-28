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

	constexpr PotentialBarrier potentialBarrier{
		0.45, 0.55, // potential wall in the middle
		3000
	};

	constexpr Ket::float_t mass = 1.0;
	constexpr auto hamiltonian = Hamiltonian<cfg.M>(mass, cfg.dx, potentialBarrier);

	OneDimensionalParticleBox<cfg.N> box(
		cfg,
		hamiltonian,
		gaussianPacket
	);
	
	Visu::VisuOscilloscope<cfg.M> visu;

	while (true)
	{
		auto p = box.evolve();
		p.normalize_with_dx(cfg.dx);
		visu.update(p);
	}

}
