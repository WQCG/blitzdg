// Copyright (C) 2017-2018  Waterloo Quantitative Consulting Group, Inc.
// See COPYING and LICENSE files at project root for more details.

/**
 * @file Poisson1d.hpp
 * @brief Short header file for 1D Poisson solver's main.cpp file.
 */
#pragma once
#include "Nodes1DProvisioner.hpp"
#include "Types.hpp"
#include <string>

namespace blitzdg {
	namespace poisson1d {
		void applyPoissonOperator(const real_matrix_type & u, const Nodes1DProvisioner & nodes1D, real_matrix_type & result);
	}
}