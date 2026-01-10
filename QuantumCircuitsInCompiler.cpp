// QuantumCircuitsInCompiler.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include "circuit.h"
#include "gates.h"
#include "iqft_gate.h"

int main()
{
	// Sanity check: Create and run a GHZ state circuit on 3 qubits
	// Expected state: GHZ state (|000> + |111>)/sqrt(2)
    constexpr auto GHZStateCircuitTest = QuantumCircuit<3>().withGates(
        QuantumGate<1>(Gates::H).toBits(0),
        QuantumGate<2>(Gates::CX).toBits(0, 1),
        QuantumGate<2>(Gates::CX).toBits(1, 2));
	GHZStateCircuitTest.printProbabilities<0, 1, 2>();

	// Construct Shor's algorithm circuit for N=21, a=2
	// 8 qubits: 3 for phase register, 5 for work register
	// The circuit performs controlled modular multiplications by
	// 2, 4, and 16 (mod 21) based on the phase qubits,
	// followed by an inverse QFT on the phase register
    constexpr auto shorCircuit_N21_a2 =
        QuantumCircuit<8>().withGates(

            // --------------------------------------------------
            // 1) Hadamard on phase register (qubits 0,1,2)
            // --------------------------------------------------
            QuantumGate<1>(Gates::H).toBits(0),
            QuantumGate<1>(Gates::H).toBits(1),
            QuantumGate<1>(Gates::H).toBits(2),

            // --------------------------------------------------
            // 2) Initialize work register to |1> (qubit 3)
            // --------------------------------------------------
            QuantumGate<1>(Gates::X).toBits(3),

            // --------------------------------------------------
            // 3) Controlled modular multiply by 2 (mod 21)
            // Control: phase qubit 0
            // Work: qubits 3..7
            // --------------------------------------------------
            QuantumGate<2>(Gates::CX).toBits(0, 3),
            QuantumGate<3>(Gates::TOFFOLI).toBits(0, 3, 4),
            QuantumGate<2>(Gates::CX).toBits(0, 4),
            QuantumGate<3>(Gates::TOFFOLI).toBits(0, 4, 5),
            QuantumGate<2>(Gates::CX).toBits(0, 5),
            QuantumGate<3>(Gates::TOFFOLI).toBits(0, 5, 6),
            QuantumGate<2>(Gates::CX).toBits(0, 6),
            QuantumGate<3>(Gates::TOFFOLI).toBits(0, 6, 7),

            // --------------------------------------------------
            // 4) Controlled modular multiply by 4 (mod 21)
            // Control: phase qubit 1
            // Work: qubits 3..7
            // --------------------------------------------------
            QuantumGate<2>(Gates::CX).toBits(1, 3),
            QuantumGate<3>(Gates::TOFFOLI).toBits(1, 3, 4),
            QuantumGate<2>(Gates::CX).toBits(1, 4),
            QuantumGate<3>(Gates::TOFFOLI).toBits(1, 4, 5),
            QuantumGate<2>(Gates::CX).toBits(1, 5),
            QuantumGate<3>(Gates::TOFFOLI).toBits(1, 5, 6),
            QuantumGate<2>(Gates::CX).toBits(1, 6),
            QuantumGate<3>(Gates::TOFFOLI).toBits(1, 6, 7),

            // --------------------------------------------------
            // 5) Controlled modular multiply by 16 (mod 21)
            // Control: phase qubit 2
            // Work: qubits 3..7
            // --------------------------------------------------
            QuantumGate<2>(Gates::CX).toBits(2, 3),
            QuantumGate<3>(Gates::TOFFOLI).toBits(2, 3, 4),
            QuantumGate<2>(Gates::CX).toBits(2, 4),
            QuantumGate<3>(Gates::TOFFOLI).toBits(2, 4, 5),
            QuantumGate<2>(Gates::CX).toBits(2, 5),
            QuantumGate<3>(Gates::TOFFOLI).toBits(2, 5, 6),
            QuantumGate<2>(Gates::CX).toBits(2, 6),
            QuantumGate<3>(Gates::TOFFOLI).toBits(2, 6, 7),

            // --------------------------------------------------
            // 6) Inverse QFT on phase register
            // --------------------------------------------------
            QuantumGate<3>(Gates::make_IQFT_matrix<3>()).toBits(0, 1, 2)
        );
    shorCircuit_N21_a2.printProbabilities<0, 1, 2>();
}

