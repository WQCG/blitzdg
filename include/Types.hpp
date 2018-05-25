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
    using index_type = int;
    using real_type = double;
    using vector_type = blitz::Array<double, 1>;
    using matrix_type = blitz::Array<double, 2>;
    using index_vector_type = blitz::Array<int, 1>;
    using index_matrix_type = blitz::Array<int, 2>;
}