// Copyright (C) 2017-2019  Waterloo Quantitative Consulting Group, Inc.
// See COPYING and LICENSE files at project root for more details.

/**
 * @file SW2d.hpp
 * @brief Short header file for 2D shallow water solver's main.cpp file.
 */
#pragma once
#include "TriangleNodesProvisioner.hpp"
#include "Types.hpp"
#include <string>

namespace blitzdg {
	namespace sw2d {
		void computeRHS(real_matrix_type& h, real_matrix_type& hu, real_matrix_type& hv, real_type g, real_matrix_type& H, real_matrix_type& Hx, real_matrix_type& Hy, real_type CD, real_type f, TriangleNodesProvisioner& triangleNodesProvisioner, real_matrix_type& RHS1, real_matrix_type& RHS2, real_matrix_type& RHS3, real_type t);
	}
}