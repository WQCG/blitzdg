// Copyright (C) 2017-2020  Waterloo Quantitative Consulting Group, Inc.
// See COPYING and LICENSE files at project root for more details.

/**
 *  @file poisson1d/main.cpp ...
 */

#include "Nodes1DProvisioner.hpp"
#include "CsvOutputter.hpp"
#include "Types.hpp"
#include "Constants.hpp"
#include "Warning.hpp"
#include "Poisson1d.hpp"
#include "LSERK4.hpp"
#include "GMRESSolver.hpp"
#include <blitz/array.h>
#include <cmath>
#include <string>
#include <stdexcept>

using blitz::firstIndex;
using blitz::secondIndex;
using blitz::thirdIndex;
using blitz::sum;
using blitzdg::normMax;
using blitzdg::poisson1d::PoissonOperator;
using blitzdg::poisson1d::Precon;
using blitzdg::fullToVector;
using blitzdg::constants::pi;
using std::abs;
using std::string;
using std::cout;
using std::endl;

int main(int argc, char **argv) {
	using namespace blitzdg;
	// Physical parameters
	const real_type xmin =-1.0;
	const real_type xmax = 1.0;

	// Numerical parameters:
	const index_type N = 4;
	const index_type K = 15;

	// Build dependencies.
	Nodes1DProvisioner nodes1DProvisioner(N, K, xmin, xmax);
	CsvOutputter outputter;

	// Pre-processing. Build grid, and get data we need for initialization.
	nodes1DProvisioner.buildNodes();
	nodes1DProvisioner.computeJacobian();

	index_type Np = nodes1DProvisioner.get_NumLocalPoints();

	const real_matrix_type & x = nodes1DProvisioner.get_xGrid();

	real_matrix_type u(Np, K);
	real_matrix_type uexact(Np, K);
	real_matrix_type RHS(Np, K);
	real_matrix_type RHSexact(Np, K);

	const index_type numGlobalNodes = Np*K;

	real_vector_type MMRHSVec(numGlobalNodes);
	real_vector_type outVec(numGlobalNodes);

	firstIndex ii;
	secondIndex jj;
	thirdIndex kk;

	printDisclaimer();

	// Intialize fields.
	uexact = sin(pi*x(ii,jj));
	RHSexact = -pi*pi*sin(pi*x(ii,jj));

	const bool byRowsOpt = false;

	real_matrix_type MassMatrix(Np,Np);
	real_matrix_type MMRHS(Np, K);

	const real_matrix_type & J = nodes1DProvisioner.get_J();
	const real_matrix_type & Vinv = nodes1DProvisioner.get_Vinv();

	MassMatrix = sum(Vinv(kk,ii)*Vinv(kk,jj), kk);
	MMRHS = J*(sum(MassMatrix(ii,kk)*RHSexact(kk,jj),kk));

	fullToVector(MMRHS, MMRHSVec, byRowsOpt);
	
	const char delim = ' ';
	outputter.writeFieldToFile("x.dat", x, delim);

	// Calculate Right-hand side for funs.
	GMRESSolver gmres;
	GMRESParams params;
	params.verbose = true;

	// Need to initialize to zero to tell gmres we aren't making an initial guess.
	outVec = 0*ii;
	GMRESOut result = gmres.solve(PoissonOperator(nodes1DProvisioner), Precon(), MMRHSVec, outVec, params);

	vectorToFull(outVec, u, byRowsOpt);

	outputter.writeFieldToFile("u.dat", u, delim);
	outputter.writeFieldToFile("uexact.dat", uexact, delim);
	outputter.writeFieldToFile("RHS.dat", RHSexact, delim);

	u -= uexact;
	double err = normMax(u);
	cout << "Error: " << err << endl;

	return 0;
} // end main