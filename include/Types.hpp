#pragma once
#include <blitz/array.h>

// Defines the basic types used throughout this project.
namespace blitzdg {
    using index_type = int;
    using real_type = double;
    using vector_type = blitz::Array<double, 1>;
    using matrix_type = blitz::Array<double, 2>;
    using index_vector_type = blitz::Array<int, 1>;
    using index_matrix_type = blitz::Array<int, 2>;
}