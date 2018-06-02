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
		void computeRHS(const real_matrix_type & u, real_type c, real_type alph, real_type nu, Nodes1DProvisioner & nodes1D, real_matrix_type & RHS);
	}
}
