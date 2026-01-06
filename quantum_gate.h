#pragma once
#include "types.h"
#include "quantum_gate_helpers.h"

template<index_t QBitCount>
class QuantumGate;

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
	friend class QuantumGate<QBitCount>;

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
        // Total number of amplitudes in the global state vector
        constexpr dimension_t N = StateCount;

        // Number of qubits the gate operates on
        constexpr dimension_t k = QBitCount;

        // Local block size: 2^k amplitudes for the affected qubits
        constexpr dimension_t dim = ConstexprMath::pow2(k);

        // Work on a copy so the operation is functional (no side effects on input)
        state_vector_t<StateCount> result = state;

        // ------------------------------------------------------------------
        // 1) Build a bitmask identifying all qubits affected by this gate.
        //
        // The `targetMask` has bits set at positions of `affectedBits`. It
        // is used to detect indices where ALL affected qubits are zero
        // (these indices serve as base indices for the block iteration).
        // ------------------------------------------------------------------
        dimension_t targetMask = 0;
        for (dimension_t q : affectedBits)
        {
            // Set the bit corresponding to qubit `q`.
            targetMask |= (dimension_t(1) << q);
        }

        // ------------------------------------------------------------------
        // 2) Temporary local state vectors for one block
        //
        // `localIn` holds the amplitudes for the 2^k local basis states
        // (for the current configuration of the unaffected qubits).
        // `localOut` receives the gate-transformed amplitudes.
        // ------------------------------------------------------------------
        state_vector_t<dim> localIn{ cplx_t::zero() };
        state_vector_t<dim> localOut{ cplx_t::zero() };

        // ------------------------------------------------------------------
        // 3) Iterate over all global basis indices
        //
        // Only process indices `base` where none of the targeted qubit bits
        // are set. Each such `base` identifies one independent block of
        // size `dim` that will be transformed.
        // ------------------------------------------------------------------
        for (dimension_t base = 0; base < N; ++base)
        {
            // If any targeted qubit is 1 in `base` skip — we only treat the
            // canonical representative (where targeted qubits are all 0).
            if (base & targetMask)
                continue;

            // --------------------------------------------------------------
            // 4) Gather amplitudes for the current block into `localIn`.
            //
            // For each local index i ∈ [0, 2^k) interpret `i` as a bitmask
            // over the k affected qubits. Map those local bits into their
            // global positions (given by `affectedBits`) and read the
            // corresponding amplitude from `result`.
            // --------------------------------------------------------------
            for (dimension_t i = 0; i < dim; ++i)
            {
                // Start from the base index (all targeted qubits zero)
                dimension_t idx = base;

                // Map local bit b of `i` to global qubit position `affectedBits[b]`.
                for (dimension_t b = 0; b < k; ++b)
                {
                    if (i & (dimension_t(1) << b))
                    {
                        // Set the bit in the global index corresponding to the
                        // b-th local qubit being 1.
                        idx |= (dimension_t(1) << affectedBits[b]);
                    }
                }

                // Copy the amplitude into the local input vector.
                localIn[i] = result[idx];
            }

            // --------------------------------------------------------------
            // 5) Apply the quantum gate to the local block:
            //
            //    localOut = gatematrix_t × localIn
            //
            // This uses the project's constexpr matrix-vector multiplication helper.
            // --------------------------------------------------------------
            localOut = apply_unitary(gateMatrix, localIn);

            // --------------------------------------------------------------
            // 6) Scatter the transformed amplitudes back into the global state.
            //
            // Reverse the mapping used in step 4 to place each `localOut[i]`
            // into the correct global index.
            // --------------------------------------------------------------
            for (dimension_t i = 0; i < dim; ++i)
            {
                dimension_t idx = base;

                // Map local index bits back to global qubit positions.
                for (dimension_t b = 0; b < k; ++b)
                {
                    if (i & (dimension_t(1) << b))
                    {
                        idx |= (dimension_t(1) << affectedBits[b]);
                    }
                }

                // Store the transformed amplitude.
                result[idx] = localOut[i];
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
template<dimension_t QBitCount>
class QuantumGate
{
	/**
	 * @brief  The unitary matrix for this gate.
	 *
	 * The matrix is owned by the factory and used to construct `QuantumGateOp`
	 * instances. Kept immutable after construction.
	 */
	const matrix_t<ConstexprMath::pow2(QBitCount), ConstexprMath::pow2(QBitCount)> gateMatrix;

public:

	/**
	 * @brief     Construct a gate from its unitary matrix.
	 * @param U   The gate matrix (must be 2^QBitCount × 2^QBitCount and unitary).
	 */
	constexpr explicit QuantumGate(
		const matrix_t<ConstexprMath::pow2(QBitCount), ConstexprMath::pow2(QBitCount)>& U)
		: gateMatrix(U) {
	}

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
        static_assert(is_valid_square_matrix(gateMatrix), "The provided matrix is not a valid square matrix.");
        static_assert(is_unitary(gateMatrix), "The provided matrix is not unitary.");
		return QuantumGateOp<QBitCount>(gateMatrix, qbit_list_t<QBitCount>{ static_cast<dimension_t>(qbits)... });
	}
};
