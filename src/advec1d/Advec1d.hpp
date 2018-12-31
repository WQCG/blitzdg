// Copyright (C) 2017-2019  Waterloo Quantitative Consulting Group, Inc.
// See COPYING and LICENSE files at project root for more details.

/**
 * @file Advec1d.hpp
 * @brief Short header file for 1D advection solver's main.cpp file.
 */
#pragma once
#include "Nodes1DProvisioner.hpp"
#include "Types.hpp"
#include <string>

namespace blitzdg {
	namespace advec1d {
		void computeRHS(const real_matrix_type & u, real_type c, Nodes1DProvisioner & nodes1D, real_matrix_type & RHS);
	}
}