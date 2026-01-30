#include "visu/visu_oscilloscope.h"
#include "systems/particle_in_a_box.h"

using namespace KetCat;

int main()
{
	constexpr OneDimensionalParticleBoxConfig<96> cfg(1.0, 1E-4);

    constexpr KetCat::float_t x0 = 0.01;
	constexpr KetCat::float_t sigma = 0.1;
	constexpr KetCat::float_t k0 = ConstexprMath::Pi * 10;

	constexpr auto gaussianPacKetCat = GaussianWavePacKetCat<cfg.M>()(x0, k0, sigma, cfg.dx);

	constexpr PotentialBarrier potentialBarrier{
		0.45, 0.55, // potential wall in the middle
		3000
	};

	constexpr KetCat::float_t mass = 1.0;
	constexpr auto hamiltonian = Hamiltonian<cfg.M>(mass, cfg.dx, potentialBarrier);

	OneDimensionalParticleBox<cfg.N> box(
		cfg,
		hamiltonian,
		gaussianPacKetCat
	);
	
	Visu::VisuOscilloscope<cfg.M> visu;

	while (true)
	{
		auto p = box.evolve();
		p.normalize_with_dx(cfg.dx);
		visu.update(p, Visu::UsePhaseEncoding::YES, Visu::ClearScreen::YES);
	}

}
