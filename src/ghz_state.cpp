#include "systems/quantum_circuit.h"
#include "visu/visu_proba_table.h"

using namespace Ket::QCC;

int main()
{
	// Sanity check: Create and run a GHZ state circuit on 3 qubits
	// Expected state: GHZ state (|000> + |111>)/sqrt(2)

	std::cout << "Greenberger-Horne-Zeilinger (GHZ) State Circuit demo (3 qubits)\n";

    constexpr auto GHZStateCircuitTest = QuantumCircuit<3>().withGates(
        QuantumGate<1, Gates::H>().toBits(0),
        QuantumGate<2, Gates::CX>().toBits(0, 1),
        QuantumGate<2, Gates::CX>().toBits(1, 2));

	Ket::Visu::VisuProbaTable<8>().update<0, 1, 2>(GHZStateCircuitTest.getStateVector());
}

