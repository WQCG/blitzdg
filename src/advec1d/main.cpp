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

#include "Nodes1DProvisioner.hpp"
#include "CsvOutputter.hpp"
#include "Types.hpp"
#include "Warning.hpp"
#include "Advec1d.hpp"
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
	const real_type c = 0.1;

	const real_type finalTime = 20.0;
	real_type t = 0.0;

	// Numerical parameters:
	const index_type N = 4;
	const index_type K = 15;
	const real_type CFL = 0.8;

	// Build dependencies.
	Nodes1DProvisioner nodes1DProvisioner(N, K, xmin, xmax);
	CsvOutputter outputter;

	// Pre-processing. Build grid, and get data we need for initialization.
	nodes1DProvisioner.buildNodes();
	nodes1DProvisioner.computeJacobian();

	index_type Np = nodes1DProvisioner.get_NumLocalPoints();

	const matrix_type & x = nodes1DProvisioner.get_xGrid();

	real_type min_dx = x(1,0) - x(0,0);

	real_type dt = CFL*min_dx/abs(c);

	matrix_type u(Np, K);
	matrix_type RHS(Np, K);
	matrix_type resRK(Np, K);

	firstIndex ii;
	secondIndex jj;

	printDisclaimer();

	// Intialize fields.
	u = exp(-10*(x(ii,jj)*x(ii,jj)));
	//u = sin(M_PI*x(ii,jj));
	RHS = 0*jj;
	resRK = 0*jj;

	index_type count = 0;

	const char delim = ' ';
	outputter.writeFieldToFile("x.dat", x, delim);

	while (t < finalTime) {
		if ((count % 10) == 0) {
			string fileName = outputter.generateFileName("u", count);
			outputter.writeFieldToFile(fileName, u, delim);
		}	

		for (index_type i=0; i < LSERK4::numStages;  i++ ) {

			// Calculate Right-hand side.
			advec1d::computeRHS(u, c, nodes1DProvisioner, RHS);

			// Compute Runge-Kutta Residual.
			resRK = LSERK4::rk4a[i]*resRK + dt*RHS;

			// Update solution.
			u += LSERK4::rk4b[i]*resRK;
		}

		real_type u_max = normMax(u);
		if ( u_max > 1e8  || std::isnan(u_max) ) {
			throw std::runtime_error("A numerical instability has occurred!");
		}

		t += dt;
		count++;
	}

	return 0;
} // end main

namespace blitzdg {
	namespace advec1d {
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

			// boundary indices.
			index_type mapO = nodes1D.get_mapO();
			index_type mapI = nodes1D.get_mapI();

			index_type numFaces = nodes1D.NumFaces;
			index_type Nfp = nodes1D.NumFacePoints;
			index_type Np = nodes1D.get_NumLocalPoints();
			index_type K = nodes1D.get_NumElements();

			real_type alpha = 0;   // 1 == central flux, 0 == upwind flux.

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
			uP(mapO) = uM(mapO); // outflow - exit stage left.
			uP(mapI) = 0;        // inflow - assumed 0 (or use exact solution at x=0).
					
			// Compute jump in flux:
			du = (uM - uP)*0.5*(c*nxVec - (1-alpha)*fabs(c*nxVec)); 

			matrix_type duMat(Nfp*numFaces, K);

			count = 0;
			for( index_type k=0; k < K; k++) {
				for ( index_type f=0; f < Nfp*numFaces; f++) {
					duMat(f,k) = du(count);
					count++;
				}
			}

			// Assumes PDE has been left-multiplied by local inverse mass matrix, so all we have left
			// is the differentiation matrix contribution, and the surface integral
			RHS =-c*rx*sum(Dr(ii,kk)*u(kk,jj), kk);

			matrix_type surfaceRHS(Np, K);
			surfaceRHS = Fscale*duMat;
			RHS += sum(Lift(ii,kk)*surfaceRHS(kk,jj), kk);
		} // computeRHS
	} // namespace advec1d
} // namespace blitzdg
