// Copyright (C) 2017-2020  Waterloo Quantitative Consulting Group, Inc.
// See COPYING and LICENSE files at project root for more details.

/**
 * @file Types.hpp
 * @brief Defines the basic types used throughout this project.
 */
#pragma once
#include <blitz/array.h>
#include <unordered_map>
#include <vector>
#include <memory>

namespace blitzdg {
    template <typename T>
    using vector_type = blitz::Array<T, 1>;

    template <typename T>
    using matrix_type = blitz::Array<T, 2>;

    template <typename T>
    using tensor3_type = blitz::Array<T, 3>;

    using index_type = int;
    using real_type = double;
    using real_vector_type = vector_type<double>;
    using real_matrix_type = matrix_type<double>;
    using index_vector_type = vector_type<int>;
    using index_matrix_type = matrix_type<int>;
    using index_tensor3_type = tensor3_type<index_type>;
    using real_tensor3_type = tensor3_type<real_type>;

    using real_mat_smart_ptr = std::unique_ptr<real_matrix_type>;
    using real_vec_smart_ptr = std::unique_ptr<real_vector_type>;

    using index_hashmap = std::unordered_map<index_type, std::vector<index_type>>;
}