// Copyright (C) 2017-2018  Derek Steinmoeller. 
// See COPYING and LICENSE files at project root for more details. 

/** \mainpage blitzdg documentation
  *
  * \section intro_sec Introduction
  *
  * This is the auto-generated API docs page for the blitzdg project.
  *
  * Click on 'Classes' to start familarizing yourself with the API
  */

#include <iostream>
#include <fstream>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>

#include <blitz/array.h>
#include <MeshManager.hpp>
#include <Nodes1DProvisioner.hpp>
#include <SparseMatrixConverter.hpp>
#include <EigenSolver.hpp>
#include <DirectSolver.hpp>

using namespace std;
using namespace blitz;
using namespace blitzdg;
using namespace boost;

// Blitz indices
firstIndex ii;
secondIndex jj;
thirdIndex kk;

void computeRHS(const matrix_type & u, const double c, Nodes1DProvisioner & nodes1D, matrix_type & RHS) {
	matrix_type & Dr = nodes1D.get_Dr();
	matrix_type & rx = nodes1D.get_rx();
	matrix_type & Lift = nodes1D.get_Lift();
	const matrix_type & Fscale = nodes1D.get_Fscale();
	const matrix_type & nx = nodes1D.get_nx();

	// Get volume to surface maps.
	const Array<int,1> & vmapM = nodes1D.get_vmapM();
	const Array<int,1> & vmapP = nodes1D.get_vmapP();

	// boundary indices;
	const index_type vmapI = nodes1D.get_vmapI();
	const index_type mapO = nodes1D.get_mapO();
	const index_type mapI = nodes1D.get_mapI();


	int numFaces = nodes1D.NumFaces;
	int Nfp = nodes1D.NumFacePoints;
	int Np = nodes1D.get_NumLocalPoints();
	int K = nodes1D.get_NumElements();

	double alpha = 0;   // 1 == central flux, 0 == upwind flux.

	matrix_type du(numFaces*Nfp, K);
	du = 0.*jj;
	matrix_type uM(numFaces*Nfp, K);

	matrix_type uCol(Np, K, ColumnMajorArray<2>());
	uCol = u; // is this gross?


	// maybe this loop can be blitz-ified...
	index_type count = 0;
	for (index_type k1=0; k1 < K; k1++) {
		for (index_type f1=0; f1 < numFaces; f1++) {
			index_type vM = vmapM(count);
			index_type vP = vmapP(count);

			double uM = uCol(vM);
			double uP = uCol(vP);

			// Inflow BC:
			if (count == mapI)
				uP = uCol(vmapI);

			// Outflow BC;
			if (count == mapO)
				uP = uM; // exit stage left.

			// Compute jump in flux:
			du(f1,k1) = (uM - uP)*0.5*(c*nx(f1,k1) - (1-alpha)*abs(c*nx(f1,k1))); 
			count++;
		}
	}

	// Assumes PDE has been left-multiplied by local inverse mass matrix, so all we have left
	// is the differentiation matrix contribution, and the surface integral
	RHS =-c*rx*(sum(Dr(ii,kk)*u(kk,jj), kk)) + Lift*(Fscale*(du));;
}

int main(int argc, char **argv) {
	
	// Physical parameters
	double xmin =-1.0;
	double xmax = 1.0;
	double c = 0.1;

	const double finalTime = 10.0;
	double t = 0.0;

	// Numerical parameters:
	int N = 4;
	int K = 20;
	double CFL = 0.05;

	// Build dependencies.
	SparseMatrixConverter matrixConverter;
	EigenSolver eigenSolver(matrixConverter);
	DirectSolver directSolver(matrixConverter);
	Nodes1DProvisioner nodes1DProvisioner(N, K, xmin, xmax, matrixConverter, eigenSolver, directSolver);
	
	// Pre-processing. Build grid, and get data we need for initialization.
	nodes1DProvisioner.buildNodes();
	nodes1DProvisioner.computeJacobian();

	int Np = nodes1DProvisioner.get_NumLocalPoints();

	matrix_type & x = nodes1DProvisioner.get_xGrid();

	double min_dx = x(1,0) - x(0,0);

	double dt = CFL*min_dx/c;

	matrix_type u(Np, K);
	matrix_type RHS(Np, K);

	u = exp(-10*(x(ii,jj)*x(ii,jj)));

	index_type count = 0;
	while (t < finalTime) {
		// Toy outputting for now.
		if ((count % 10) == 0) {
			ofstream outFile;
			outFile.open("advec1d." + lexical_cast<string>(count));
			outFile << u;
			outFile.close();
		}

		// Calculate Righ-hand side at current time-step.
		computeRHS(u, c, nodes1DProvisioner, RHS);

		// Forward Euler time-step for now, to be replaced.
		u = u + dt*RHS;

		t += dt;
		count++;
	}

	return 0;
}