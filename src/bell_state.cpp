#include <iostream>

#include "engine/circuit.h"
#include "gates/common_gates.h"
#include "gates/iqft_gate.h"


int main()
{
	// Sanity check: Create and run a GHZ state circuit on 3 qubits
	// Expected state: GHZ state (|000> + |111>)/sqrt(2)

	std::cout << "Bell State Circuit demo (2 qubits)\n";

    constexpr auto BellStateCircuit = QuantumCircuit<2>().withGates(
        QuantumGate<1, Gates::H>().toBits(0),
        QuantumGate<2, Gates::CX>().toBits(0, 1));

    BellStateCircuit.printProbabilities<0, 1>();
}

