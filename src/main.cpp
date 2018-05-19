// Copyright (C) 2017-2018  Waterloo Quantitative Consulting Group, Inc.
// See COPYING and LICENSE files at project root for more details.

/** \mainpage blitzdg documentation
  *
  * \section intro_sec Introduction
  *
  * This is the auto-generated API docs page for the blitzdg project.
  *
  * Click on 'Classes' to start familarizing yourself with the API
  */

#include "MeshManager.hpp"
#include "Nodes1DProvisioner.hpp"
#include "SparseMatrixConverter.hpp"
#include "EigenSolver.hpp"
#include "DirectSolver.hpp"
#include "Types.hpp"
#include "Warning.hpp"
#include <blitz/array.h>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>

using blitz::ColumnMajorArray;
using blitz::firstIndex;
using blitz::secondIndex;
using blitz::thirdIndex;
using blitz::sum;
using std::endl;
using std::ofstream;
using std::setfill;
using std::setw;
using std::string;
using std::stringstream;

namespace blitzdg {
	void computeRHS(const matrix_type & u, real_type c, Nodes1DProvisioner & nodes1D, matrix_type & RHS) {
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

		// boundary indices;
		index_type vmapI = nodes1D.get_vmapI();
		index_type mapO = nodes1D.get_mapO();
		index_type mapI = nodes1D.get_mapI();


		index_type numFaces = nodes1D.NumFaces;
		index_type Nfp = nodes1D.NumFacePoints;
		index_type Np = nodes1D.get_NumLocalPoints();
		index_type K = nodes1D.get_NumElements();

		real_type alpha = 0;   // 1 == central flux, 0 == upwind flux.

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

				real_type uM = uCol(vM);
				real_type uP = uCol(vP);

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
		RHS =-c*rx*sum(Dr(ii,kk)*u(kk,jj), kk) + Lift*Fscale*du;
	}

	void writeFieldToFile(const string fileName, const matrix_type & field ) {
		ofstream outFile;
		outFile.open(fileName);
		for(index_type i=0; i < field.rows(); i++) {
			for(index_type k=0; k < field.cols(); k++) {
				real_type num = field(i, k);
				outFile << num << " ";
			}
			outFile << endl;
		}
		outFile.close();
	}
} // namespace blitzdg

int main(int argc, char **argv) {
	using namespace blitzdg;
	// Physical parameters
	const real_type xmin =-1.0;
	const real_type xmax = 1.0;
	const real_type c = 0.1;

	const real_type finalTime = 10.0;
	real_type t = 0.0;

	// Numerical parameters:
	const index_type N = 4;
	const index_type K = 100;
	const real_type CFL = 0.05;

	// Build dependencies.
	Nodes1DProvisioner nodes1DProvisioner(N, K, xmin, xmax);
	
	// Pre-processing. Build grid, and get data we need for initialization.
	nodes1DProvisioner.buildNodes();
	nodes1DProvisioner.computeJacobian();

	index_type Np = nodes1DProvisioner.get_NumLocalPoints();

	const matrix_type & x = nodes1DProvisioner.get_xGrid();

	real_type min_dx = x(1,0) - x(0,0);

	real_type dt = CFL*min_dx/abs(c);

	matrix_type u(Np, K);
	matrix_type RHS(Np, K);

	firstIndex ii;
	secondIndex jj;
	u = exp(-10*(x(ii,jj)*x(ii,jj)));

	printDisclaimer();

	index_type count = 0;

	writeFieldToFile("x.dat", x);

	while (t < finalTime) {
		// Toy outputting for now.
		if ((count % 10) == 0) {
			stringstream fileNameStrm;
			fileNameStrm << "u" << setfill('0') << setw(7) << count << ".dat";

			writeFieldToFile(fileNameStrm.str(), u);
		}

		// Calculate Righ-hand side at current time-step.
		computeRHS(u, c, nodes1DProvisioner, RHS);

		// Forward Euler time-step for now, to be replaced.
		u += dt*RHS;

		t += dt;
		count++;
	}

	return 0;
} // end main