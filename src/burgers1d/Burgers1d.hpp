// Copyright (C) 2017-2018  Waterloo Quantitative Consulting Group, Inc.
// See COPYING and LICENSE files at project root for more details.

/**
 * @file Burgers1d.hpp
 * @brief Short header file for 1D Burger's solver's main.cpp file.
 */
#pragma once
#include "Nodes1DProvisioner.hpp"
#include "Types.hpp"
#include <string>

namespace blitzdg {
    namespace burgers1d {
        real_type Burgers2(real_type x, real_type t, real_type alpha, real_type nu, real_type c);
        void Burgers2(real_matrix_type & u, const real_matrix_type & x, real_type t, real_type alpha, real_type nu, real_type c);
        void computeRHS(const real_matrix_type & u, const real_matrix_type & x, real_type t, real_type c, real_type alpha, real_type nu, Nodes1DProvisioner & nodes1D, real_matrix_type & RHS);
    }
}
