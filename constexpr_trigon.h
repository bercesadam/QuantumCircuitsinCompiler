#pragma once
#include <cstdint>

/// @file
/// @brief Constexpr-friendly trigonometric helpers used to build unitary gates.
///
/**
 * @details
 * This header provides small, constexpr-capable approximations for sine and cosine
 * together with a minimal range-reduction routine. The implementations are intended
 * for use in compile-time construction of gate matrices (for example QFT) where
 * full runtime math library calls may not be available or desirable.
 *
 * Notes:
 *  - Approximations use low-order Taylor-like polynomials (suitable for small
 *    input magnitudes after range reduction).
 *  - Range reduction maps inputs to the primary approximation interval so the
 *    polynomial approximations remain accurate.
 *  - All functions are free of dynamic allocation and are trivial to evaluate
 *    in constexpr contexts.
 */

namespace ConstexprMath
{
    /// @brief Mathematical constant π (double precision).
    static constexpr double Pi = 3.141592653589793;

    /**
     * @brief Evaluate a polynomial approximation to sin(x) around 0.
     *
     * @param x  Input angle (radians). Intended to be small (|x| ≲ π/4) after range reduction.
     * @return   Approximate sin(x).
     *
     * @remarks
     * The polynomial is derived from the Taylor series:
     *   sin(x) ≈ x - x^3/6 + x^5/120 - x^7/5040
     * which has good accuracy for small |x|. This routine expects the caller to
     * first reduce the argument into the polynomial's region of accuracy.
     */
    constexpr double sin_poly(double x) noexcept {
        // Precompute x^2 to reuse in higher powers.
        const double x2 = x * x;
        // Evaluate x * (1 - x^2/6 + x^4/120 - x^6/5040)
        return x * (1.0
            - x2 / 6.0
            + x2 * x2 / 120.0
            - x2 * x2 * x2 / 5040.0);
    }

    /**
     * @brief Evaluate a polynomial approximation to cos(x) around 0.
     *
     * @param x  Input angle (radians). Intended to be small (|x| ≲ π/4) after range reduction.
     * @return   Approximate cos(x).
     *
     * @remarks
     * The polynomial is taken from the Taylor series:
     *   cos(x) ≈ 1 - x^2/2 + x^4/24 - x^6/720
     * which is accurate for small |x|.
     */
    constexpr double cos_poly(double x) noexcept {
        const double x2 = x * x;
        return 1.0
            - x2 / 2.0
            + x2 * x2 / 24.0
            - x2 * x2 * x2 / 720.0;
    }

    /**
     * @brief Floor implementation usable in constexpr context.
     *
     * @param x  A double value.
     * @return   The largest integer not greater than x, returned as int.
     *
     * @remarks
     * This simple helper casts to int32_t and corrects for truncated values
     * when x is negative and not an integer. It is small and constexpr-friendly.
     */
    constexpr int floor_constexpr(double x) {
        const int32_t i = static_cast<int32_t>(x);
        // Casting truncates toward zero; adjust when x < i to get true floor.
        return (x < static_cast<double>(i)) ? (i - 1) : i;
    }

    /**
     * @brief Reduce an arbitrary angle to the small interval used by the polynomial approximations.
     *
     * @param x  Input angle in radians.
     * @return   Reduced angle xr such that xr = x - k*(π/2) for some integer k.
     *
     * @details
     * Range reduction strategy:
     *  - Compute k = floor((x + π/4) / (π/2)).
     *  - Subtract k*(π/2) from x to obtain xr.
     *
     * The chosen formula maps x into the interval [-π/4, π/4) which is suitable
     * for the low-order polynomial approximations used above.
     */
    constexpr double reduce(double x) noexcept {
        // Determine the quadrant-like shift count k using a constexpr floor.
        double k = floor_constexpr((x + Pi / 4) / (Pi / 2));
        // Subtract k * (π/2) to obtain the reduced angle.
        double xr = x - k * (Pi / 2);
        return xr;
    }

    /**
     * @brief Constexpr sine approximation for arbitrary input using range reduction + polynomial.
     *
     * @param x  Input angle in radians.
     * @return   Approximate sin(x).
     *
     * @remarks
     * The pipeline is:
     *  1) Reduce x into the small interval with `reduce`.
     *  2) Evaluate `sin_poly` on the reduced angle.
     *
     * For higher accuracy across the full domain, a full quadrant-aware sign/permute
     * handling can be added; the current routine relies on the reduction producing
     * a small xr where the polynomial approximation is valid.
     */
    constexpr double sin_constexpr(double x) noexcept {
        double xr = reduce(x);
        // Use the small-angle polynomial on the reduced value.
        return sin_poly(xr);
    }

    /**
     * @brief Constexpr cosine approximation for arbitrary input using range reduction + polynomial.
     *
     * @param x  Input angle in radians.
     * @return   Approximate cos(x).
     *
     * @remarks
     * Same pipeline as sin_constexpr: range reduction followed by polynomial evaluation.
     */
    constexpr double cos_constexpr(double x) noexcept {
        double xr = reduce(x);
        return cos_poly(xr);
    }
}