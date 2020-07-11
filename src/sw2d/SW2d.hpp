// Copyright (C) 2017-2020  Waterloo Quantitative Consulting Group, Inc.
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
		struct physParams {
			const real_type g = 9.81;
			real_type CD;
			real_type f;
			real_type initTime;
			real_type finalTime;
		};

			// Numerical parameters 
		struct numParams {
			index_type N ; // N = Order of polynomials
			real_type CFL;
			index_type outputInterval;
			real_type filterPercent = 0.95;
			index_type filterOrder = 4;
		};

		// Physical fields
		struct fields {
			real_matrix_type h;
			real_matrix_type hu;
			real_matrix_type hv;
			real_matrix_type H;
			real_matrix_type Hx;
			real_matrix_type Hy;
			
			// 'non-conservative' variables.
			real_matrix_type eta;
			real_matrix_type u;
			real_matrix_type v;

			// RHS and RK Residual storage
			real_matrix_type RHS1, RHS2, RHS3;
			real_matrix_type resRK1, resRK2, resRK3;
		};

		// helper declarations.
		double computeTimeStep(fields& fds, const physParams& phys, const numParams& num, const DGContext2D& dg);
		void computeRHS(fields fds, const numParams& num, const physParams& phys, const DGContext2D& dg, real_type t);
		void readDepthData(const std::string& depthFile, real_matrix_type& H);
		void buildSpongeCoeff(const DGContext2D& dg, real_type spongeStrength, real_type radInfl, real_matrix_type& spongeCoeff);
	}
}