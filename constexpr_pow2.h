#pragma once
#include <concepts>

/// @file
/// @brief Small constexpr integer utilities used for dimensions and bit manipulations.
///
/**
 * @details
 * This header provides:
 *  - `pow2(n)` : compute 2^n at compile time via bit-shift (returns an unsigned integral).
 *  - `is_power_of_two(x)` : test whether an integer is a (positive) power of two.
 *
 * Both functions are constexpr and noexcept and intended for use in compile-time
 * dimension computations (matrix sizes, state vector lengths) and static assertions.
 */

namespace ConstexprMath
{
    /**
     * @brief Compute 2 raised to the power of `n` at compile time.
     *
     * @tparam UIntType  Unsigned integral type used for the exponent and result.
     * @param n          Exponent (non-negative).
     * @return 2^n as `UIntType`.
     *
     * @note The implementation uses a left bit-shift on `1` which is defined for
     *       unsigned integral types and is constexpr-friendly.
     *
     * Example:
     * @code
     * constexpr auto four = ConstexprMath::pow2<std::size_t>(2); // four == 4
     * @endcode
     */
    template <std::unsigned_integral UIntType>
    constexpr UIntType pow2(UIntType n) noexcept
    {
        // Left-shift 1 by n positions: 1 << n == 2^n for unsigned types.
        return UIntType{ 1 } << n;
    }

    /**
     * @brief Determine whether a value is an exact power of two.
     *
     * @tparam UIntType  Unsigned integral type for the argument.
     * @param x          Value to test.
     * @return true if x is a power of two (1, 2, 4, 8, ...); false otherwise.
     *
     * @note This uses the classic bit trick: x > 0 && (x & (x - 1)) == 0.
     *       The operation is constexpr and does not allocate memory.
     */
    template <std::unsigned_integral UIntType>
    constexpr bool is_power_of_two(UIntType x) noexcept
    {
        // Fast test: powers of two have exactly one bit set.
        return x > 0 && (x & (x - 1)) == 0;
    }

}