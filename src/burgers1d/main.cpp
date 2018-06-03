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
#include "DenseMatrixHelpers.hpp"
#include "Warning.hpp"
#include "Burgers1d.hpp"
#include "LSERK4.hpp"
#include "LinAlgHelpers.hpp"
#include <blitz/array.h>
#include <cmath>
#include <string>
#include <stdexcept>

using blitz::firstIndex;
using blitz::secondIndex;
using blitz::thirdIndex;
using blitz::sum;
using std::abs;
using std::string;
using std::cout;
using std::endl;

int main(int argc, char **argv) {
	using namespace blitzdg;
	// Physical parameters
	const real_type xmin =-5.0;
	const real_type xmax = 5.0;

	const real_type alpha = 1.0;
	const real_type nu = 0.1;
	const real_type c = 0.5;

	const real_type finalTime = 0.5;
	real_type t = 0.0;

	// Numerical parameters:
	const index_type N = 4;
	const index_type K = 40;
	const real_type CFL = 0.01;

	// Build dependencies.
	Nodes1DProvisioner nodes1DProvisioner(N, K, xmin, xmax);
	CsvOutputter outputter;

	// Pre-processing. Build grid, and get data we need for initialization.
	nodes1DProvisioner.buildNodes();
	nodes1DProvisioner.computeJacobian();

	index_type Np = nodes1DProvisioner.get_NumLocalPoints();

	const real_matrix_type & x = nodes1DProvisioner.get_xGrid();

	real_type min_dx = x(1,0) - x(0,0);

	real_type dt = CFL*min_dx/abs(c);

	real_matrix_type u(Np, K);
	real_matrix_type RHS(Np, K);
	real_matrix_type resRK(Np, K);
	real_matrix_type ue(u.shape());

	firstIndex ii;
	secondIndex jj;

	printDisclaimer();

	// Intialize u field at t=0
    burgers1d::Burgers2(u, x, t, alpha, nu, c);

	RHS = 0*jj;
	resRK = 0*jj;

	index_type count = 0;

	const char delim = ' ';
	outputter.writeFieldToFile("x.dat", x, delim);

	while (t < finalTime) {
		if ((count % 10) == 0) {
			string fileName = outputter.generateFileName("u", count);
			outputter.writeFieldToFile(fileName, u, delim);

            // Why doesn't this work?
			//string fileName2 = outputter.generateFileName("ue", count);
            //burgers1d::Burgers2(ue, x, t, alpha, nu, c);
			//outputter.writeFieldToFile(fileName2, ue, delim);
		}	

		for (index_type i=0; i < LSERK4::numStages;  i++ ) {

			// Calculate Right-hand side.
			burgers1d::computeRHS(u, x, t, c, alpha, nu, nodes1DProvisioner, RHS);

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


    burgers1d::Burgers2(ue, x, t, alpha, nu, c);
	u -= ue;
	double err = normMax(u);
	cout << "Error: " << err << endl;

	return 0;
} // end main

namespace blitzdg {
	namespace burgers1d {
        real_type Burgers2(real_type x, real_type t, real_type alpha, real_type nu, real_type c) {
            // Eq (2) from https://www.hindawi.com/journals/mpe/2015/414808/
            return (c/alpha) - (c/alpha)*tanh( (c/(2.0*nu)) * (x - c*t) );
        }
        void Burgers2(real_matrix_type & u, const real_matrix_type & x, real_type t, real_type alpha, real_type nu, real_type c) {
             for (index_type i=0; i < u.extent(0); i++) {
               for (index_type j=0; j <u.extent(1); j++) {
                 u(i,j) = Burgers2(x(i,j), t, alpha, nu, c);
               }
             }
        }

		void computeRHS(const real_matrix_type & u, const real_matrix_type & x, real_type t, real_type c, real_type alpha, real_type nu, Nodes1DProvisioner & nodes1D, real_matrix_type & RHS) {
			// Blitz indices
			firstIndex ii;
			secondIndex jj;
			thirdIndex kk;
			const real_matrix_type& Dr = nodes1D.get_Dr();
			const real_matrix_type& rx = nodes1D.get_rx();
			const real_matrix_type& Lift = nodes1D.get_Lift();
			const real_matrix_type& Fscale = nodes1D.get_Fscale();
			const real_matrix_type& nx = nodes1D.get_nx();

			// Get volume to surface maps.
			const index_vector_type& vmapM = nodes1D.get_vmapM();
			const index_vector_type& vmapP = nodes1D.get_vmapP();

			// boundary indices.
			index_type mapO = nodes1D.get_mapO();
			index_type mapI = nodes1D.get_mapI();
			index_type vmapO = nodes1D.get_vmapO();
			index_type vmapI = nodes1D.get_vmapI();

			index_type numFaces = nodes1D.NumFaces;
			index_type Nfp = nodes1D.NumFacePoints;
			index_type Np = nodes1D.get_NumLocalPoints();
			index_type K = nodes1D.get_NumElements();

			real_vector_type du(numFaces*Nfp*K);
			real_vector_type uM(numFaces*Nfp*K);
			real_vector_type uP(numFaces*Nfp*K);
			real_vector_type dq(numFaces*Nfp*K);
			real_vector_type qM(numFaces*Nfp*K);
			real_vector_type qP(numFaces*Nfp*K);
			real_vector_type nxVec(numFaces*Nfp*K);
			real_vector_type xVec(Np*K);
			real_vector_type uVec(Np*K);
			real_vector_type qVec(Np*K);

            real_type uL, uR;
			du = 0.*ii;
			uM = 0.*ii;
			uP = 0.*ii;
			nxVec = 0.*ii;
			xVec = 0.*ii;

			// We want to apply maps to column-wise ordering of the nodes.
			const bool byRowsOpt = false;

			fullToVector(nx, nxVec, byRowsOpt);
			fullToVector(u, uVec, byRowsOpt);
			fullToVector(x, xVec, byRowsOpt);
			applyIndexMap(uVec, vmapM, uM);
			applyIndexMap(uVec, vmapP, uP);

            // Following implementation of
            // https://github.com/dsteinmo/euler-nullproj/blob/master/NUDG/Codes1D/BurgersRHS1D.m

 			// BCs
            uL = Burgers2(xVec(vmapI), t, alpha, nu, c);
			uR = Burgers2(xVec(vmapO), t, alpha, nu, c);

			// Compute jump in flux
			du = uM - uP;

            // Why factor of 2 here?
            du(mapI) = 2*(uVec(vmapI)-uL);
            du(mapO) = 2*(uVec(vmapO)-uR);

			real_matrix_type duMat(Nfp*numFaces, K);
			vectorToFull(du, duMat, byRowsOpt);

            // Derivatives
            real_matrix_type ux(u.shape());
            ux  = rx*sum(Dr(ii,kk)*u(kk,jj), kk);

            // compute q and q stuff
			real_matrix_type tmp(Nfp*numFaces, K);
			real_matrix_type q(u.shape());

            tmp = 0.5*Fscale*nx*duMat;
            q = nu*(ux - sum(Lift(ii,kk)*tmp(kk,jj), kk));

			fullToVector(q, qVec, byRowsOpt);
			applyIndexMap(qVec, vmapM, qM);
			applyIndexMap(qVec, vmapP, qP);

            dq = (qM-qP)*0.5;
            dq(mapI) = 0.0;
            dq(mapO) = 0.0;

            // Evaluate nonlinear flux, boundary conditions
            real_vector_type du2(numFaces*Nfp*K);
            du2 = 0.5*(uM*uM - uP*uP);
            du2(mapI) = uVec(vmapI)*uVec(vmapI) - uL*uL;
            du2(mapO) = uVec(vmapO)*uVec(vmapO) - uR*uR;

            real_type maxvel = max(abs(u)); 

            real_vector_type flux(numFaces*Nfp*K); 
            flux = nxVec*(du2/2.0 - nu*dq) - (maxvel/2.0)*du;

			real_matrix_type fluxMat(Nfp*numFaces, K);
			vectorToFull(flux, fluxMat, byRowsOpt);

            real_matrix_type dfdr(u.shape());
            real_matrix_type tmp2(u.shape());
            tmp2 = (0.5*u*u - nu*q);
            dfdr = sum(Dr(ii,kk)*tmp2(kk,jj), kk);

            // Form RHS
    		RHS = -rx*dfdr;

            real_matrix_type tmp3(fluxMat.shape());
            tmp3 = Fscale*fluxMat;

			RHS += sum(Lift(ii,kk)*tmp3(kk,jj), kk);
		} // computeRHS
	} // namespace burgers1d
} // namespace blitzdg
