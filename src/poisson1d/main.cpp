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
#include "DenseMatrixHelpers.hpp"
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

	firstIndex ii;
	secondIndex jj;

	printDisclaimer();

	// Intialize fields.
	u = 0*jj;
	uexact = sin(M_PI*x(ii,jj));
	RHS = -M_PI*M_PI*sin(M_PI*x(ii,jj));
	
	const char delim = ' ';
	outputter.writeFieldToFile("x.dat", x, delim);

	// Calculate Right-hand side for funs.
	poisson1d::applyPoissonOperator(uexact, nodes1DProvisioner, RHS);

	return 0;
} // end main

namespace blitzdg {
	namespace poisson1d {
		void applyPoissonOperator(const real_matrix_type & u, const Nodes1DProvisioner & nodes1D, real_matrix_type & result) {
			// Blitz indices
			firstIndex ii;
			secondIndex jj;
			thirdIndex kk;
			const real_matrix_type& Dr = nodes1D.get_Dr();
			const real_matrix_type& rx = nodes1D.get_rx();
			const real_matrix_type& J = nodes1D.get_J();
			const real_matrix_type& Lift = nodes1D.get_Lift();
			const real_matrix_type& Fscale = nodes1D.get_Fscale();
			const real_matrix_type& nx = nodes1D.get_nx();
			const real_matrix_type& Vinv = nodes1D.get_Vinv();

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

			real_matrix_type VinvTrans(Np, Np);
			VinvTrans = Vinv(jj,ii);

			real_vector_type du(numFaces*Nfp*K);
			real_vector_type uM(numFaces*Nfp*K);
			real_vector_type uP(numFaces*Nfp*K);
			real_vector_type nxVec(numFaces*Nfp*K);

			du = 0.*ii;
			uM = 0.*ii;
			uP = 0.*ii;

			const bool byRowsOpt = false;
			fullToVector(nx, nxVec, byRowsOpt);

			real_matrix_type ux(Np, K);
			real_vector_type uVec(Np*K);
			real_vector_type uxVec(Np*K);

			fullToVector(u, uVec, byRowsOpt);

			ux = rx*sum(Dr(ii,kk)*u(kk,jj), kk);
			fullToVector(ux, uxVec, byRowsOpt);

			applyIndexMap(uVec, vmapM, uM);
			applyIndexMap(uVec, vmapP, uP);

			// impose boundary condition -- Dirichlet BC's
			uP(mapI) = -uM(mapI);
			uP(mapO) = -uM(mapO);

			// Compute jump in flux:
			du = (uM - uP);

			real_matrix_type duMat(Nfp*numFaces, K);
			vectorToFull(du, duMat, byRowsOpt);

			real_matrix_type surfaceRHS(Nfp*numFaces, K);
			real_matrix_type q(Np, K);
			real_vector_type qVec(Np, K);

			surfaceRHS = Fscale*(nx*duMat/2.0);

			q  = rx*(Dr(ii,kk)*u(kk,jj));
			q -= sum(Lift(ii,kk)*(surfaceRHS(kk,jj)), kk);

			fullToVector(q, qVec, byRowsOpt);
			real_vector_type dq(numFaces*Nfp*K);
			real_vector_type qM(numFaces*Nfp*K);
			real_vector_type qP(numFaces*Nfp*K);
			real_vector_type uxM(numFaces*Nfp*K);
			real_vector_type uxP(numFaces*Nfp*K);

			applyIndexMap(qVec, vmapM, qM);
			applyIndexMap(uxVec, vmapM, uxM);
			applyIndexMap(uxVec, vmapP, uxP);

			// Neumann BC's on q
			uxP(mapI) = uxM(mapI);
			uxP(mapO) = uxM(mapO);

			qP = 0.5*(uxM + uxP);

			dq = qM -qP;

			real_matrix_type dqMat;
			vectorToFull(dq, dqMat, byRowsOpt);

			const double hmin = 2.0/normMax(rx);
			const double tau = (Np*Np)/hmin;

			surfaceRHS = nx*(dqMat + tau*nx*duMat);

			real_matrix_type MassMatrix(Np, Np);
			MassMatrix = sum(Vinv(kk,ii)*Vinv(kk,jj), kk);

			real_matrix_type tmp(Np,K);

			tmp  = rx*(Dr(ii,kk)*q(kk,jj));
			tmp -= sum(Lift(ii,kk)*(surfaceRHS(kk,jj)), kk);
			result = J*(MassMatrix(ii,kk)*tmp(jj,kk));
		} // computeRHS
	} // namespace poisson1d
} // namespace blitzdg
