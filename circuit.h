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

    // ------------------------------------------------------------
    // Prints measurement probabilities of all computational states
    //
    // Example output:
    // |00> : 50.00 %
    // |11> : 50.00 %
    // ------------------------------------------------------------
    /// @brief Print the probability (in percent) of each computational basis state.
    /// @param os  Output stream (defaults to std::cout).
    void printProbabilities(std::ostream& os = std::cout) const
    {
        constexpr dimension_t StateCount = ConstexprMath::pow2(QBitCount);

        // Format numbers with two decimal places
        os << std::fixed << std::setprecision(2);

        for (dimension_t i = 0; i < StateCount; ++i)
        {
            const cplx_t& amp = StateVector[i];

            // Probability = |amplitude|^2 = re^2 + im^2
            const double probability =              
                static_cast<double>(amp.re * amp.re +
                    amp.im * amp.im);

            // Print basis state in binary ket notation: |0101>
            os << "|";

            for (dimension_t bit = 0; bit < QBitCount; ++bit)
            {
                // Mask selects bit from most-significant to least for human-readable order
                const dimension_t mask = dimension_t(1) << (QBitCount - 1 - bit);
                os << ((i & mask) ? '1' : '0');
            }

            os << "> : " << (probability * 100.0) << " %\n";
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
