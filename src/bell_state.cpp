#include "systems/quantum_circuit.h"
#include "visu/visu_proba_table.h"

using namespace KetCat::QCC;

int main()
{
	// Sanity check: Create and run a GHZ state circuit on 3 qubits
	// Expected state: GHZ state (|000> + |111>)/sqrt(2)

	std::cout << "Bell State Circuit demo (2 qubits)\n";

    constexpr auto BellStateCircuit = QuantumCircuit<2>().withGates(
        QuantumGate<1, Gates::H>().toBits(0),
        QuantumGate<2, Gates::CX>().toBits(0, 1));

    KetCat::Visu::VisuProbaTable<4>().update<0, 1>(BellStateCircuit.getStateVector());
}

