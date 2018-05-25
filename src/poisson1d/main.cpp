// Copyright (C) 2017-2018  Waterloo Quantitative Consulting Group, Inc.
// See COPYING and LICENSE files at project root for more details.

/**
 *  @file poisson1d/main.cpp ...
 */

#include "Nodes1DProvisioner.hpp"
#include "CsvOutputter.hpp"
#include "Types.hpp"
#include "Warning.hpp"
#include "Poisson1d.hpp"
#include "LSERK4.hpp"
#include "LinAlgHelpers.hpp"
#include <blitz/array.h>
#include <cmath>
#include <string>
#include <stdexcept>

using blitz::ColumnMajorArray;
using blitz::firstIndex;
using blitz::secondIndex;
using blitz::thirdIndex;
using blitz::sum;
using std::abs;
using std::string;
using std::cout;

int main(int argc, char **argv) {
	using namespace blitzdg;
	// Physical parameters
	const real_type xmin =-1.0;
	const real_type xmax = 1.0;

	// Numerical parameters:
	const index_type N = 4;
	const index_type K = 15;

	const real_type tau = 1.0;

	// Build dependencies.
	Nodes1DProvisioner nodes1DProvisioner(N, K, xmin, xmax);
	CsvOutputter outputter;

	// Pre-processing. Build grid, and get data we need for initialization.
	nodes1DProvisioner.buildNodes();
	nodes1DProvisioner.computeJacobian();

	index_type Np = nodes1DProvisioner.get_NumLocalPoints();

	const matrix_type & x = nodes1DProvisioner.get_xGrid();

	real_type min_dx = x(1,0) - x(0,0);

	matrix_type u(Np, K);
	matrix_type uexact(Np, K);
	matrix_type RHS(Np, K);

	firstIndex ii;
	secondIndex jj;

	printDisclaimer();

	// Intialize fields.
	u = 0*jj;
	uexact = sin(M_PI*x(ii,jj));
	RHS = -M_PI*M_PI*sin(M_PI*x(ii,jj));
	
	const char delim = ' ';
	outputter.writeFieldToFile("x.dat", x, delim);

	// Calculate Right-hand side.
	poisson1d::applyPoissonOperator(u, nodes1DProvisioner, tau, RHS);

	return 0;
} // end main

namespace blitzdg {
	namespace poisson1d {
		void applyPoissonOperator(const matrix_type & u, const Nodes1DProvisioner & nodes1D, const real_type tau, matrix_type & result) {
			// Blitz indices
			firstIndex ii;
			secondIndex jj;
			thirdIndex kk;
			const matrix_type& Dr = nodes1D.get_Dr();
			const matrix_type& rx = nodes1D.get_rx();
			const matrix_type& Lift = nodes1D.get_Lift();
			const matrix_type& Fscale = nodes1D.get_Fscale();
			const matrix_type& nx = nodes1D.get_nx();

			// Get volume to surface maps.
			const index_vector_type& vmapM = nodes1D.get_vmapM();
			const index_vector_type& vmapP = nodes1D.get_vmapP();

			// boundary indices.
			index_type mapO = nodes1D.get_mapO();
			index_type mapI = nodes1D.get_mapI();

			index_type numFaces = nodes1D.NumFaces;
			index_type Nfp = nodes1D.NumFacePoints;
			index_type Np = nodes1D.get_NumLocalPoints();
			index_type K = nodes1D.get_NumElements();

			vector_type du(numFaces*Nfp*K);
			vector_type uM(numFaces*Nfp*K);
			vector_type uP(numFaces*Nfp*K);
			vector_type nxVec(numFaces*Nfp*K);

			du = 0.*ii;
			uM = 0.*ii;
			uP = 0.*ii;
			nxVec = 0.*ii;

			matrix_type nxCol(numFaces*Nfp,K, ColumnMajorArray<2>());
			nxCol = nx;		

			matrix_type uCol(Np, K, ColumnMajorArray<2>());
			uCol = u;

			matrix_type uxCol(Np, K, ColumnMajorArray<2>());
			uxCol = rx*sum(Dr(ii,kk)*u(kk,jj), kk);

			index_type count = 0;
			for( index_type k=0; k < K; k++) {
				for ( index_type f=0; f < Nfp*numFaces; f++) {
					nxVec(count) = nxCol(f,k);
					uM(count) = uCol(vmapM(count));
					uP(count) = uCol(vmapP(count));
					count++;
				}
			}

			// BC's
			uP(mapO) = uM(mapO);
			uP(mapI) = 0;
					
			// Compute jump in flux:
			du = (uM - uP); 

			matrix_type duMat(Nfp*numFaces, K);

			count = 0;
			for( index_type k=0; k < K; k++) {
				for ( index_type f=0; f < Nfp*numFaces; f++) {
					duMat(f,k) = du(count);
					count++;
				}
			}

			matrix_type surfaceRHS(Np, K);
			surfaceRHS = Fscale*duMat;
			result = sum(Lift(ii,kk)*surfaceRHS(kk,jj), kk);
		} // computeRHS
	} // namespace poisson1d
} // namespace blitzdg
