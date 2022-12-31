// Copyright (C) 2017-2022  Waterloo Quantitative Consulting Group, Inc.
// See COPYING and LICENSE files at project root for more details.

/**
 * @file ins2d.hpp
 * @brief Short header file for 2D incompressible Navier-Stokes solver's main.cpp file.
 */
#pragma once
#include "QuadNodesProvisioner.hpp"
#include "Types.hpp"
#include <string>

namespace blitzdg {
	namespace ins2d {
		void computeRHS(real_matrix_type h, real_matrix_type hu, real_matrix_type hv, real_type g, QuadNodesProvisioner& quadNodesProvisioner, real_matrix_type& RHS1, real_matrix_type& RHS2, real_matrix_type& RHS3);
	}
}