#pragma once
#include <concepts>


/// @file
/// @brief Constexpr-capable complex number type used throughout the project.
///
/**
 * @details
 * This header defines `ConstexprMath::Complex<T>`, a small constexpr-friendly
 * complex number type with a minimal set of operations required by the
 * simulator (construction, addition, subtraction, multiplication, in-place
 * accumulation and conjugation). All members and functions are `constexpr`
 * and `noexcept` where appropriate so they can be used in compile-time
 * contexts (constexpr initialization of gate matrices, etc.).
 */

namespace ConstexprMath
{
    /**
     * @brief Simple constexpr complex number for numeric computations.
     *
     * @tparam T  Floating-point type to store real and imaginary parts (e.g. double).
     *
     * The type exposes public data members `re` and `im` and lightweight
     * arithmetic operators implemented in a constexpr-friendly way. It is
     * intentionally minimal (no exceptions, no heap allocations).
     */
    template <std::floating_point T>
    struct Complex
    {
        /// @brief The real part of the complex number.
        T re = 0.0;
        /// @brief The imaginary part of the complex number.
        T im = 0.0;

        /// @brief Default constructor producing 0 + 0i.
        constexpr Complex() noexcept = default;

        /// @brief Construct a complex number from a real value (imaginary part = 0).
        /// @param r  Real part.
        constexpr Complex(T r) noexcept : re(r), im(0.0) {}

        /// @brief Construct a complex number from real and imaginary parts.
        /// @param r  Real part.
        /// @param i  Imaginary part.
        constexpr Complex(T r, T i) noexcept : re(r), im(i) {}

        /// @brief Return the additive identity 0 + 0i.
        /// @return Complex zero.
        static constexpr Complex zero() noexcept { return {}; }

        /// @brief Create a complex number from a real value (helper alias).
        /// @param r  Real part.
        /// @return Complex number r + 0i.
        static constexpr Complex fromReal(T r) noexcept { return { r, 0.0 }; }

        /// @brief Return the imaginary unit +i (0 + 1i).
        /// @return Complex representing i.
        static constexpr Complex plus_i() noexcept { return { 0.0, 1.0 }; }

        /// @brief Return the negative imaginary unit -i (0 - 1i).
        /// @return Complex representing -i.
        static constexpr Complex minus_i() noexcept { return { 0.0, -1.0 }; }

        /// @brief Component-wise addition.
        /// @param other  Addend.
        /// @return The sum (this + other).
        ///
        /// Example:
        /// @code
        /// Complex<double> a{1.0, 2.0};
        /// Complex<double> b{0.5, -1.0};
        /// auto c = a + b; // c == {1.5, 1.0}
        /// @endcode
        constexpr Complex operator+ (Complex other) const noexcept
        {
            // Add real and imaginary parts separately.
            return { re + other.re, im + other.im };
        }

        /// @brief Component-wise subtraction.
        /// @param other  Subtrahend.
        /// @return The difference (this - other).
        constexpr Complex operator- (Complex other) const noexcept
        {
            // Subtract real and imaginary parts separately.
            return { re - other.re, im - other.im };
        }

        /// @brief Complex multiplication following (a+bi)*(c+di) = (ac - bd) + (ad + bc)i.
        /// @param other  The multiplier.
        /// @return The product.
        constexpr Complex operator* (Complex other) const noexcept
        {
            // Compute using the standard formula to avoid temporaries where possible.
            return { re * other.re - im * other.im, re * other.im + im * other.re };
        }

        /// @brief In-place addition (accumulate another complex into this).
        /// @param other  Value to add.
        /// @return A reference to *this after addition (returned by value here for constexpr convenience).
        constexpr Complex operator+= (Complex other) noexcept
        {
            // Mutate both components; returning *this by value is acceptable for constexpr use.
            re += other.re;
            im += other.im;
            return *this;
        }

		/// @brief Compute the complex conjugate.
		/// @return The conjugate of this complex number.
        constexpr Complex conj() const noexcept
        {
            // Conjugate: invert the sign of the imaginary part.
            return { re, -im };
		}

		/// @brief Compute the squared magnitude (norm) of the complex number.
		/// @return The squared norm (re^2 + im^2).
        constexpr T normSquared() const noexcept
        {
            // Return the squared magnitude (norm) of the complex number.
            return re * re + im * im;
		}
    };
}