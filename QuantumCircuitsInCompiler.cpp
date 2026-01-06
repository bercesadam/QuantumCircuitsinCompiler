// QuantumCircuitsInCompiler.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include "circuit.h"
#include "gates.h"
#include "cnot_gates.h"

int main()
{
    constexpr auto bellStateCircuitTest = QuantumCircuit<2>().withGates(
        QuantumGate<1>(Gates::H).toBits(0),
        QuantumGate<2>(Gates::CX).toBits(1, 0));
    bellStateCircuitTest.printProbabilities();
}

