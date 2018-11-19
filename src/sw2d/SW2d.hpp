// Copyright (C) 2017-2018  Waterloo Quantitative Consulting Group, Inc.
// See COPYING and LICENSE files at project root for more details.

/**
 * @file SW2d.hpp
 * @brief Short header file for 2D shallow water solver's main.cpp file.
 */
#pragma once
#include "Nodes1DProvisioner.hpp"
#include "Types.hpp"
#include <string>

namespace blitzdg {
	namespace sw2d {
		void computeRHS(const real_matrix_type & u, real_type g, Nodes1DProvisioner & nodes1D, real_matrix_type & RHS);
	}
}