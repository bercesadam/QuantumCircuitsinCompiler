#include <chrono>
#include <iostream>
#include <thread>

#include "constexprmath/constexpr_trigon.h"

#include "hamiltonian/hamiltonian.h"
#include "hamiltonian/potential_barrier.h"
#include "wavefunction/gaussian_wave_packet.h"
#include "wavefunction/hydrogen.h"

#include "solvers/crank_nicolson_solver.h"

#ifdef _WIN32
#include <windows.h>
#undef max
#endif



template<dimension_t Dim>
void visualize(const state_vector_t<Dim>& s)
{
	static const char* bars[] = {
	"\xE2\x96\x81", "\xE2\x96\x82", "\xE2\x96\x83",
	"\xE2\x96\x84", "\xE2\x96\x85", "\xE2\x96\x86",
	"\xE2\x96\x87", "\xE2\x96\x88"
	};

	double maxProba = 0;
	for (const cplx_t& a : s)
	{
		maxProba = std::max(maxProba, a.normSquared());
	}

	std::cout << "Probability density: |";
	for (const cplx_t& a : s)
	{
		double p = a.normSquared() / maxProba;
		std::size_t idx = static_cast<std::size_t>(p * 7);
		std::cout << bars[idx];
	}
	std::cout << "|\n";

	std::cout << "Phase:               |";
	for (const cplx_t& a : s)
	{
		double phase = std::atan2(a.im, a.re); // -π..π
		double t = (phase + ConstexprMath::Pi) / (2 * ConstexprMath::Pi);
		std::size_t idx = static_cast<std::size_t>(t * 5);
		std::cout << bars[idx];
	}
	std::cout << "|\n";
}


int main()
{
#ifdef _WIN32
	SetConsoleOutputCP(CP_UTF8);
#endif

	// box length
	constexpr double L = 1.0;
	constexpr unsigned int N = 50;
	constexpr unsigned int M = N - 2; // Dirichlet peremfeltétel

    constexpr double x0 = 0.01;
	constexpr double sigma = 0.1;
	constexpr double k0 = 10 * ConstexprMath::Pi;
  
    constexpr double dx = L / (N - 1);
	constexpr double dt = 2E-4;

	constexpr auto gaussianPacket = GaussianWavePacket<M>()(x0, k0, sigma, dx); //(1, dx, L);//GaussianWavePacket<M>()(x0, k0, sigma, dx);
	constexpr auto hydrogen = HydrogenOrbitals<M>()(QuantumNumber::_2s(), 1.0, dx, 0.0);

	constexpr PotentialBarrier potential{
		0.45, 0.55, // potential wall in the middle
		1500 // Joule
	};
	constexpr HamiltonConstants constants{
		1.0,    // hBar
		1.0,    // m
		dx      // dx
	};
	constexpr auto hamiltonian = Hamiltonian<M>(constants);//, potential);

	constexpr CrankNicolsonSolver<M> evol(hamiltonian, dt);

	using namespace std::chrono_literals;
	auto res = gaussianPacket;
	while(true)
	{
		std::cout << "\x1B[2J\x1B[H";
		visualize(res.getVector());
		res = evol(res);
		
		std::this_thread::sleep_for(100ms);
	}

}
