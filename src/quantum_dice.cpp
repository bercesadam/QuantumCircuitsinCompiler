#include "systems/quantum_circuit.h"
#include "visu/visu_proba_table.h"

using namespace Ket::QCC;

int main()
{
	std::cout << "Fair dice from 3 qubits\n";

	// P(1) to 1/3, P(0) to 2/3, so sin?(?/2) = 1/3 -> ? = 2*arcsin(sqrt(1/3)) = 1.23095941734
	constexpr auto Ry13_23 = Gates::RotationY(1.23095941734);

	constexpr auto RyPi4 = Gates::RotationY(ConstexprMath::Pi / 4);
	constexpr auto RyMinusPi4 = Gates::RotationY(-ConstexprMath::Pi / 4);

    constexpr auto FairDice = QuantumCircuit<3>().withGates(
        QuantumGate<1, Ry13_23>().toBits(2),
		QuantumGate<1, Gates::X>().toBits(2),
		QuantumGate<1, RyPi4>().toBits(1),
        QuantumGate<2, Gates::CX>().toBits(2, 1),
		QuantumGate<1, RyMinusPi4>().toBits(1),
		QuantumGate<2, Gates::CX>().toBits(2, 1),
		QuantumGate<1, Gates::X>().toBits(1),
		QuantumGate<1, Gates::H>().toBits(0)
	);

    Ket::Visu::VisuProbaTable<8>().update<0, 1, 2>(FairDice.getStateVector());
}

