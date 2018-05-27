// Copyright (C) 2017-2018  Waterloo Quantitative Consulting Group, Inc.
// See COPYING and LICENSE files at project root for more details.

/**
 * @file Types.hpp
 * @brief Defines the basic types used throughout this project.
 */
#pragma once
#include <blitz/array.h>

#ifndef M_PI
#define M_PI (3.141592653589793238)
#endif

namespace blitzdg {
    template <typename T>
    using vector_type = blitz::Array<T, 1>;

    template <typename T>
    using matrix_type = blitz::Array<T, 2>;

    using index_type = int;
    using real_type = double;
    using real_vector_type = vector_type<double>;
    using real_matrix_type = matrix_type<double>;
    using index_vector_type = vector_type<int>;
    using index_matrix_type = matrix_type<int>;
}