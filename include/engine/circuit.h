#pragma once
#include "quantum_gate.h"

#include <iostream>
#include <iomanip>

/// @file
/// @brief Small executor and convenience API for composing and running quantum gates.
///
/**
 * @details
 * `QuantumCircuit` is a lightweight compile-time helper to execute a sequence of
 * gate-like callables on an initialized state vector. The executor is templated
 * on the number of qubits and accepts any callable satisfying the `QuantumGateLike`
 * concept (i.e. callables that accept and return a `state_vector_t<N>`).
 */

template<index_t QBitCount>
class QuantumCircuit;

/// @brief Concept for types that behave like a gate: they are callable with a state vector.
/// @tparam G  Type to test.
///
/// The requirement checks one instantiation (with StateCount == 1) and uses that as a
/// proxy for gate-like behaviour. Implementers should ensure their operator() is
/// templated or overloaded appropriately for various `StateCount` values.
template<typename G>
concept QuantumGateLike =
    requires(G g)
{
    // Try invoking the gate with a sample 1-qubit state vector; the returned type
    // must be a state_vector_t<1>. This models the "callable that maps state->state".
    { g(state_vector_t<1>{}) } -> std::same_as<state_vector_t<1>>;
};

/// @brief Executor that schedules and runs a series of gate-like callables at compile-time (where possible).
/// @tparam QBitCount  Number of qubits in the circuit.
/// @tparam Gates      Variadic pack of gate-like callables accepted by the executor.
template<dimension_t QBitCount, QuantumGateLike... Gates>
class QuantumCircuitExecutor {

    /// @brief The internal global state vector (amplitudes for 2^QBitCount basis states).
    state_vector_t<ConstexprMath::pow2(QBitCount)> StateVector;

    /// @brief Construct executor and immediately execute provided gates.
    /// @param gates  Variadic list of gate-like callables to apply in order.
    constexpr QuantumCircuitExecutor(const Gates& ... gates)
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
	/// @param os              Output stream to write to (defaults to std::cout).
    ///    
    template<dimension_t... SelectedQBits>
    void printProbabilities(std::ostream& os = std::cout) const
    {
        constexpr dimension_t StateCount = ConstexprMath::pow2(QBitCount);
        constexpr dimension_t numQBits = sizeof...(SelectedQBits);
        constexpr qbit_list_t<numQBits> selectedQubits =
            qbit_list_t<numQBits>{ static_cast<dimension_t>(SelectedQBits)... };

        // Format numbers
        os << std::fixed << std::setprecision(10);

        // Number of qubits we are considering
        constexpr dimension_t numSelected = numQBits;
        constexpr dimension_t reducedDim = ConstexprMath::pow2(numSelected);

        // Array to accumulate probabilities for each reduced basis state
        std::array<double, reducedDim> reducedProbabilities = {};

        // Sum probabilities over all other qubits
        for (dimension_t i = 0; i < StateCount; ++i)
        {
            double p = static_cast<double>(
                StateVector[i].re * StateVector[i].re +
                StateVector[i].im * StateVector[i].im
                );

            dimension_t reducedIdx = 0;
            for (dimension_t b = 0; b < numSelected; ++b)
            {
                const dimension_t q = selectedQubits[b];
                if (i & (1ULL << q))
                    reducedIdx |= (1ULL << b);
            }

            reducedProbabilities[reducedIdx] += p;
        }

        // Print header
        os << "| Binary | Decimal | Probability (%) |\n";
        os << "|--------|---------|----------------|\n";

        // Print reduced probabilities
        for (dimension_t i = 0; i < reducedDim; ++i)
        {
            // Binary string
            os << "|";
            for (dimension_t b = numSelected; b-- > 0;)
                os << ((i & (1ULL << b)) ? '1' : '0');
            os << "|";

            // Decimal
            os << " " << i << " |";

            // Probability
            os << " " << (reducedProbabilities[i] * 100.0) << " % |\n";
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
