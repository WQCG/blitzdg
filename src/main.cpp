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
#include <blitz/array.h>
#include <MeshManager.hpp>
#include <Nodes1DProvisioner.hpp>
#include <SparseMatrixConverter.hpp>
#include <EigenSolver.hpp>
#include <DirectSolver.hpp>

using namespace std;
BZ_USING_NAMESPACE(blitz)
using namespace blitzdg;

// Blitz indices
firstIndex ii;
secondIndex jj;
thirdIndex kk;

void computeRHS(const Array<double,2> & u, const double c, Nodes1DProvisioner & nodes1D, Array<double,2> & RHS) {
	//something like: RHS = -c*rx*(Dr*u) + c*LIFT*Fscale*numFlux
	Array<double,2> & Dr = nodes1D.get_Dr();
	Array<double,2> & rx = nodes1D.get_rx();
	Array<double,2> & Lift = nodes1D.get_Lift();
	const Array<double,2> & Fscale = nodes1D.get_Fscale();
	const Array<double,2> & nx = nodes1D.get_nx();

	const Array<int,1> & vmapM = nodes1D.get_vmapM();
	const Array<int,1> & vmapP = nodes1D.get_vmapP();

	int numFaces = nodes1D.NumFaces;
	int Nfp = nodes1D.NumFacePoints;
	int Np = nodes1D.get_NumLocalPoints();
	int K = nodes1D.get_NumElements();

	double alpha = 0;   //1 == central flux, 0 == upwind flux.

	Array<double,2> du(numFaces*Nfp, K);
	du = 0.*jj;
	Array<double,2> uM(numFaces*Nfp, K);

	Array<double,2> uCol(Np, K, ColumnMajorArray<2>());
	uCol = u; // is this gross?

	index_type count = 0;
	for (index_type k1=0; k1 < K; k1++) {
		for (index_type f1=0; f1 < numFaces; f1++) {
			index_type vM = vmapM(count);
			index_type vP = vmapP(count);

			// Compute jump in flux:
			du(f1,k1) = (uCol(vM) - uCol(vP))*0.5*(c*nx(f1,k1) - (1-alpha)*abs(c*nx(f1,k1))); 
			count++;
		}
	}

	RHS =-c*rx*(sum(Dr(ii,jj)*u(jj,kk), jj)) + Lift*(Fscale*(du));
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
	double CFL = 0.5;

	// Build dependencies.
	SparseMatrixConverter matrixConverter;
	EigenSolver eigenSolver(matrixConverter);
	DirectSolver directSolver(matrixConverter);
	Nodes1DProvisioner nodes1DProvisioner(N, K, xmin, xmax, matrixConverter, eigenSolver, directSolver);
	
	// Pre-processing. Build grid, and get data we need for initialization.
	nodes1DProvisioner.buildNodes();
	nodes1DProvisioner.computeJacobian();

	int Np = nodes1DProvisioner.get_NumLocalPoints();

	Array<double,2> & x = nodes1DProvisioner.get_xGrid();

	double min_dx = x(1,0) - x(0,0);

	double dt = CFL*min_dx/c;

	Array<double, 2> u(Np, K);
	Array<double, 2> RHS(Np, K);

	u = exp(-5*(x(ii)*x(ii)));
	while (t < finalTime) {
		// Calculate Righ-hand side at current time-step.
		computeRHS(u, c, nodes1DProvisioner, RHS);

		// Forward Euler time-step for now, to be replaced.
		u = u + dt*RHS;
		t += dt;
	}

	return 0;
}