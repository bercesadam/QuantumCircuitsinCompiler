#pragma once
#include <iostream>
#include <iomanip>

#include "wavefunction/state_vector.h"
#include "solvers/quantum_gate_solver.h"


/// @file
/// @brief Small executor and convenience API for composing and running quantum gates.
///

/// @brief Concept for types that behave like a gate: they are callable with a state vector.
/// @tparam GateType  Type to test.
///
/// The requirement checks one instantiation (with StateCount == 1) and uses that as a
/// proxy for gate-like behaviour. Implementers should ensure their operator() is
/// templated or overloaded appropriately for various `StateCount` values.
template<typename GateType>
concept QuantumGateLike =
    requires(GateType g)
{
    // Try invoking the gate with a sample 1-qubit state vector; the returned type
    // must be a state_vector_t<1>. This models the "callable that maps state->state".
    { g(state_vector_t<1>{}) } -> std::same_as<state_vector_t<1>>;
};

/// @brief Executor that schedules and runs a series of gate-like callables at compile-time (where possible).
/// @tparam QBitCount  Number of qubits in the circuit.
/// @tparam Gates      Variadic pack of gate-like callables accepted by the executor.
template<dimension_t QBitCount>
class QuantumCircuit {

    /// @brief Precompute 2^QBitCount for convenience
    constexpr dimension_t BasisStateCount = ConstexprMath::pow2(QBitCount);

    /// @brief The internal global state vector (amplitudes for 2^QBitCount basis states).
    StateVector<BasisStateCount> _stateVector;


    /// @brief Construct executor and immediately execute provided gates.
    /// @param gates  Variadic list of gate-like callables to apply in order.
    constexpr QuantumCircuitExecutor()
    {
        // Initialize to the |0...0> computational basis state
        StateVector = {};
        StateVector[0] = cplx_t::fromReal(1.0);

        // Apply the provided gates in sequence
        executeCircuit(gates...);
    }

    friend class QuantumCircuit<QBitCount>;

public:
    /// @brief Recursively apply gates: head then recurse on tail.
    /// @tparam Gate  First gate type.
    /// @tparam Rest  Remaining gate types.
    /// @param gate   Gate callable to apply.
    /// @param rest   Remaining gate callables.
    template<typename Gate, typename... Rest>
    constexpr void executeCircuit(const Gate& gate, const Rest&... rest)
    {
        static_assert(QuantumGateLike<Gate>);

        // Apply the gate to the current state vector (gate returns a new vector)
        StateVector = gate(StateVector);

        // Recurse for the remaining gates
        executeCircuit(rest...);
    }

    /// @brief Base case for recursion: no gates left to apply.
    constexpr void executeCircuit() {}

	/// @brief Print the probabilities of measuring selected qubits in each basis state.
	/// @tparam SelectedQBits  Indices of the qubits to consider (0-based).
    ///    
    template<dimension_t... SelectedQBits>
    void printProbabilities() const
    {
        constexpr dimension_t StateCount = ConstexprMath::pow2(QBitCount);
        constexpr dimension_t NumSelected = sizeof...(SelectedQBits);
        constexpr qbit_list_t<NumSelected> SelectedQubits =
            qbit_list_t<NumSelected>{ static_cast<dimension_t>(SelectedQBits)... };

        // Format numbers
        std::cout << std::fixed << std::setprecision(2);

        // Number of qubits we are considering
        constexpr dimension_t ReducedDim = ConstexprMath::pow2(NumSelected);

        // Array to accumulate probabilities for each reduced basis state
        std::array<double, ReducedDim> ReducedProbabilities = {};

        // Sum probabilities over all other qubits
        for (dimension_t i = 0; i < StateCount; ++i)
        {
			double p = StateVector[i].normSquared();

            dimension_t ReducedIdx = 0;
            for (dimension_t b = 0; b < NumSelected; ++b)
            {
                const dimension_t q = SelectedQubits[b];
                if (i & (1ULL << q))
                    ReducedIdx |= (1ULL << b);
            }

            ReducedProbabilities[ReducedIdx] += p;
        }

        // Print header
        std::cout << "| Binary | Decimal | Probability (%) |\n";
        std::cout << "|--------|---------|----------------|\n";

        // Print reduced probabilities
        for (dimension_t i = 0; i < ReducedDim; ++i)
        {
            // Binary string
            std::cout << "|";
            for (dimension_t b = NumSelected; b-- > 0;)
                std::cout << ((i & (1ULL << b)) ? '1' : '0');
            std::cout << "|";

            // Decimal
            std::cout << " " << i << " |";

            // Probability
            std::cout << " " << (ReducedProbabilities[i] * 100.0) << " % |\n";
        }
    }
};

/// @brief Facade class to create executors bound to a fixed qubit count.
template<index_t QBitCount>
class QuantumCircuit
{
public:
    /// @brief Create an executor with the provided gate sequence.
    /// @tparam Gates  Gate-like callables to include in the circuit.
    /// @param gates   Instances of the gate-like callables (passed by reference-to-const).
    /// @return        A `QuantumCircuitExecutor` that has already executed the gates on an initialized state.
    template<QuantumGateLike... Gates>
    constexpr QuantumCircuitExecutor<QBitCount, Gates...> withGates(const Gates& ... gates) const
    {
        return QuantumCircuitExecutor<QBitCount, Gates...>(gates...);
    }
};
