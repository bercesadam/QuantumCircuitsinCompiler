#pragma once

#include "core_types.h"


/// @brief Represents a Hilbert space of a given dimension.
/// @details This class encapsulates the concept of a Hilbert space
/// with a specified dimension to provide a type-safe way to handle
/// state vector dimensions.

template <typename Derived>
class Visualization
{
public:
    template <dimension_t Dim>
    void update(const state_vector_t<N>& m) {
        static_cast<Derived*>(this)->update_impl(m);
    }
};