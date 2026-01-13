#pragma once
#include "types.h"
#include "quantum_gate_helpers.h"

template<dimension_t QBitCount, auto GateMatrix>
	requires (is_gate_matrix_v<std::remove_cvref_t<decltype(GateMatrix)>> && is_unitary<GateMatrix>())
struct QuantumGate;

/**
 * @brief     Represents an application of a quantum gate to a specific set of qubits.
 *
 * @tparam QBitCount  Number of qubits the gate matrix acts on (compile-time).
 *
 * This class is a lightweight, non-copyable, non-movable handle that stores
 * a compile-time sized unitary matrix for the gate and a fixed list of
 * target qubit indices. Calling the object with a global state vector
 * applies the gate to the specified qubits and returns the transformed
 * global state vector (functional style).
 */
template<dimension_t QBitCount>
class QuantumGateOp
{
	/**
	 * @brief  The unitary matrix that defines the gate.
	 *
	 * Size: 2^QBitCount × 2^QBitCount. Stored as a compile-time sized matrix type.
	 */
	const matrix_t<ConstexprMath::pow2(QBitCount), ConstexprMath::pow2(QBitCount)> gateMatrix;

	/**
	 * @brief  Fixed-size list of qubit indices affected by this gate.
	 *
	 * The order of indices in this list determines the mapping between the
	 * local basis bits ([0 .. 2^QBitCount)) and the global qubit positions.
	 */
	const qbit_list_t<QBitCount> affectedBits;

	// Delete copy constructor and copy assignment
	QuantumGateOp(const QuantumGateOp&) = delete;
	QuantumGateOp& operator=(const QuantumGateOp&) = delete;

	// Delete move constructor and move assignment
	QuantumGateOp(QuantumGateOp&&) = delete;
	QuantumGateOp& operator=(QuantumGateOp&&) = delete;

	/**
	 * @brief     Construct an operation from a gate matrix and affected qubits.
	 * @param U   The unitary gate matrix (must be square and unitary).
	 * @param affectedBits  The list of target qubit indices.
	 *
	 * Both `U` and `affectedBits` are stored by value (immutable after construction).
	 * Static assertions verify matrix shape and unitarity at compile-time where possible.
	 */
	constexpr QuantumGateOp(
		const matrix_t<ConstexprMath::pow2(QBitCount), ConstexprMath::pow2(QBitCount)>& U,
		const qbit_list_t<QBitCount>& affectedBits)
		: gateMatrix(U), affectedBits(affectedBits)
	{

	}

	// Allow the factory QuantumGate to access the private ctor
	template<dimension_t _QBitCount, auto _GateMatrix>
		requires (is_gate_matrix_v<std::remove_cvref_t<decltype(_GateMatrix)>>&& is_unitary<_GateMatrix>())
	friend struct QuantumGate;

public:
	/**
	 * @brief     Apply the stored gate to a global state vector and return the result.
	 *
	 * @tparam StateCount  Dimension of the global state vector (2^number_of_global_qubits).
	 * @param state        The input global state vector (passed by value — functional style).
	 * @return             A new global state vector with the gate applied to `affectedBits`.
	 *
	 * The algorithm partitions the global state into independent blocks of size 2^QBitCount.
	 * Each block corresponds to a fixed assignment of the unaffected qubits; the affected
	 * qubits span the local 2^QBitCount basis inside each block. For each block we:
	 *  1) gather the local amplitudes into `localIn`,
	 *  2) compute `localOut = gatematrix_t * localIn`,
	 *  3) write `localOut` back into the corresponding positions of the global state.
	 *
	 * Only indices where all targeted qubits are zero are treated as block "bases" to avoid
	 * redundant processing.
	 */
	template<dimension_t StateCount>
	constexpr state_vector_t<StateCount>
		operator()(state_vector_t<StateCount> state) const
	{
		constexpr dimension_t k = QBitCount;
		constexpr dimension_t dim = ConstexprMath::pow2(k);

		state_vector_t<StateCount> result = state;

		// ------------------------------------------------------------
		// Precompute masks for affected qubits
		// ------------------------------------------------------------

		dimension_t affectedMask = 0;
		std::array<dimension_t, k> bitPos{};

		for (dimension_t i = 0; i < k; ++i) {
			bitPos[i] = affectedBits[i];
			affectedMask |= (dimension_t(1) << bitPos[i]);
		}

		// ------------------------------------------------------------
		// Iterate over all block bases
		// Each base corresponds to one fixed configuration of
		// the non-affected qubits
		// ------------------------------------------------------------

		for (dimension_t base = 0; base < StateCount; ++base) {

			// Only process each block once:
			// base must have all affected qubits = 0
			if (base & affectedMask)
				continue;

			// --------------------------------------------------------
			// Gather local (2^k) amplitudes
			// --------------------------------------------------------

			state_vector_t<dim> localIn{};

			for (dimension_t localIdx = 0; localIdx < dim; ++localIdx) {

				dimension_t globalIdx = base;

				// Map local basis |b0 b1 ... bk-1>
				// to physical qubits affectedBits[]
				for (dimension_t q = 0; q < k; ++q) {
					if (localIdx & (dimension_t(1) << q)) {
						globalIdx |= (dimension_t(1) << bitPos[q]);
					}
				}

				localIn[localIdx] = result[globalIdx];
			}

			// --------------------------------------------------------
			// Apply k-qubit unitary in the local subspace
			// --------------------------------------------------------

			const state_vector_t<dim> localOut =
				apply_unitary(gateMatrix, localIn);

			// --------------------------------------------------------
			// Scatter results back to the full statevector
			// --------------------------------------------------------

			for (dimension_t localIdx = 0; localIdx < dim; ++localIdx) {

				dimension_t globalIdx = base;

				for (dimension_t q = 0; q < k; ++q) {
					if (localIdx & (dimension_t(1) << q)) {
						globalIdx |= (dimension_t(1) << bitPos[q]);
					}
				}

				result[globalIdx] = localOut[localIdx];
			}
		}

		return result;
	}
};

/**
 * @brief     Factory for creating `QuantumGateOp` objects from a gate matrix.
 *
 * @tparam QBitCount  Number of qubits the gate matrix acts on.
 *
 * `QuantumGate` owns the gate matrix and exposes `toBits(...)` to bind the
 * matrix to a concrete set of qubit indices producing a `QuantumGateOp`.
 */
template<dimension_t QBitCount, auto GateMatrix>
	requires (is_gate_matrix_v<std::remove_cvref_t<decltype(GateMatrix)>> && is_unitary<GateMatrix>())
struct QuantumGate
{
	/**
	 * @brief     Construct a gate from its unitary matrix.
	 * @param U   The gate matrix (must be 2^QBitCount × 2^QBitCount and unitary).
	 *//*
	constexpr explicit QuantumGate(
		const matrix_t<ConstexprMath::pow2(QBitCount), ConstexprMath::pow2(QBitCount)>& U)
		: gateMatrix(U) {
	}*/

	/**
	 * @brief     Bind this gate to a list of qubit indices and return an operation.
	 *
	 * @tparam QBits  Variadic list of indices convertible to `dimension_t`. The
	 *                number of provided indices must equal `QBitCount`.
	 * @param qbits   The qubit indices to which the gate will be applied.
	 * @return         A `QuantumGateOp<QBitCount>` object ready to apply the gate.
	 *
	 * Example:
	 *   auto h2 = QuantumGate<1>(HADAMARD).toBits(1);
	 */
	template<std::convertible_to<dimension_t>... QBits>
	constexpr QuantumGateOp<QBitCount> toBits(QBits... qbits) const
	{
		static_assert(sizeof...(qbits) == QBitCount);
        //static_assert(is_valid_square_matrix(GateMatrix), "The provided matrix is not a valid square matrix.");
       // static_assert(is_unitary(GateMatrix), "The provided matrix is not unitary.");
		return QuantumGateOp<QBitCount>(GateMatrix, qbit_list_t<QBitCount>{ static_cast<dimension_t>(qbits)... });
	}
};


/// @brief     Factory for creating `QuantumGateOp` objects from a 1-qubit gate matrix.
/// @tparam QBitCount  Number of qubits the gate matrix acts on.
template<dimension_t QBitCount, auto GateMatrix>
	requires (is_1_qbit_gate_matrix_v<std::remove_cvref_t<decltype(GateMatrix)>>&& is_unitary<GateMatrix>())
struct ParallelSingleQubitGate
{
	//TODO
	template<std::convertible_to<dimension_t>... QBits>
	constexpr auto toBits(QBits... qbits) const
	{
		static_assert(sizeof...(qbits) == QBitCount);

		return;
	}
};

