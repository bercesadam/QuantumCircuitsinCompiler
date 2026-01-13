#include <iostream>

#include "engine/circuit.h"
#include "gates/common_gates.h"
#include "gates/iqft_gate.h"


int main()
{
	// Sanity check: Create and run a GHZ state circuit on 3 qubits
	// Expected state: GHZ state (|000> + |111>)/sqrt(2)

	std::cout << "Greenberger-Horne-Zeilinger (GHZ) State Circuit demo (3 qubits)\n";

    constexpr auto GHZStateCircuitTest = QuantumCircuit<3>().withGates(
        QuantumGate<1, Gates::H>().toBits(0),
        QuantumGate<2, Gates::CX>().toBits(0, 1),
        QuantumGate<2, Gates::CX>().toBits(1, 2));

	GHZStateCircuitTest.printProbabilities<0, 1, 2>();
}

