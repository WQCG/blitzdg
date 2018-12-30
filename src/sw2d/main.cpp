// Copyright (C) 2017-2018  Waterloo Quantitative Consulting Group, Inc.
// See COPYING and LICENSE files at project root for more details.

#include "Nodes1DProvisioner.hpp"
#include "CsvOutputter.hpp"
#include "Types.hpp"
#include "BlitzHelpers.hpp"
#include "Warning.hpp"
#include "SW2d.hpp"
#include "LSERK4.hpp"
#include "LinAlgHelpers.hpp"
#include "TriangleNodesProvisioner.hpp"
#include <blitz/array.h>
#include <math.h>
#include <string>
#include <stdexcept>

using blitz::firstIndex;
using blitz::secondIndex;
using blitz::thirdIndex;
using blitz::sum;
using blitz::sqrt;
using std::abs;
using std::string;
using std::cout;
using std::endl;
using std::sqrt;
using std::max;

int main(int argc, char **argv) {
	using namespace blitzdg;

	printDisclaimer();

	// Physical parameters
	const real_type g = 9.81;

	const real_type finalTime = 2.0;
	real_type t = 0.0;

	// Numerical parameters (N = Order of polynomials)
	const index_type N = 4;
	const real_type CFL = 0.8;

	// Build dependencies.
	MeshManager meshManager;
	meshManager.readMesh("input/coarse_box.msh");
	const index_type K = meshManager.get_NumElements();

	// Dependency-inject mesh manager to nodes provisioner.
	TriangleNodesProvisioner triangleNodesProvisioner(N, meshManager);

	// Pre-processing steps.
	triangleNodesProvisioner.buildNodes();
	triangleNodesProvisioner.buildLift();
	triangleNodesProvisioner.buildPhysicalGrid();
	triangleNodesProvisioner.buildMaps();

	CsvOutputter outputter;

	const real_matrix_type& x = triangleNodesProvisioner.get_xGrid();
	const real_matrix_type& y = triangleNodesProvisioner.get_yGrid();
	index_type Np = triangleNodesProvisioner.get_NumLocalPoints();

	real_matrix_type H(Np,K);
	real_matrix_type eta(Np, K);
	real_matrix_type h(Np, K);
	real_matrix_type u(Np, K);
	real_matrix_type v(Np, K);
	real_matrix_type hu(Np, K);
	real_matrix_type hv(Np, K);
	real_matrix_type RHS1(Np, K), RHS2(Np, K), RHS3(Np, K);
	real_matrix_type resRK1(Np, K), resRK2(Np, K), resRK3(Np, K);

	firstIndex ii;
	secondIndex jj;

	// Intialize fields.
	H = 10.0 + 0*jj;
	eta = exp(-10.0*(x*x));
	u = 0*jj;
	v = 0*jj;
	h = H + eta;
	hu = h*u;
	hv = h*v;

	real_matrix_type Fscale = triangleNodesProvisioner.get_Fscale();
	
	const real_type dt = CFL/max((N+1)*(N+1)*0.5*abs(Fscale)*sqrt(g*h));

	cout << "dt=" << dt << endl;

	RHS1 = 0*jj;
	RHS2 = 0*jj;
	RHS3 = 0*jj;

	resRK1= 0*jj;
	resRK2= 0*jj;
	resRK3= 0*jj;

	index_type count = 0;

	const char delim = ' ';
	outputter.writeFieldToFile("x.dat", x, delim);
	outputter.writeFieldToFile("y.dat", y, delim);

	while (t < finalTime) {
		if ((count % 10) == 0) {
			eta = h-H;
			string fileName = outputter.generateFileName("eta", count);
			outputter.writeFieldToFile(fileName, eta, delim);
		}	

		for (index_type i=0; i < LSERK4::numStages;  i++ ) {

			// Calculate Right-hand side.
			sw2d::computeRHS(h, hu, hv, g, triangleNodesProvisioner, RHS1, RHS2, RHS3);

			// Compute Runge-Kutta Residual.
			resRK1 = LSERK4::rk4a[i]*resRK1 + dt*RHS1;
			resRK2 = LSERK4::rk4a[i]*resRK2 + dt*RHS2;
			resRK3 = LSERK4::rk4a[i]*resRK3 + dt*RHS3;

			// Update solution.
			h  += LSERK4::rk4b[i]*resRK1;
			hu += LSERK4::rk4b[i]*resRK2;
			hv += LSERK4::rk4b[i]*resRK3;
		}

		real_type u_max = normMax(u);
		if ( u_max > 1e8  || std::isnan(u_max) ) {
			throw std::runtime_error("A numerical instability has occurred!");
		}

		t += dt;
		count++;
	}
	real_matrix_type etafinal(Np, K);
	etafinal = h - H;

	return 0;
} // end main

namespace blitzdg {
	namespace sw2d {
		void computeRHS(real_matrix_type h, real_matrix_type hu, real_matrix_type hv, real_type g, TriangleNodesProvisioner& triangleNodesProvisioner, real_matrix_type& RHS1, real_matrix_type& RHS2, real_matrix_type& RHS3) {
			// Blitz indices
			firstIndex ii;
			secondIndex jj;
			thirdIndex kk;

			// Differentiation matrices and scaling factors.
			const real_matrix_type& Dr = triangleNodesProvisioner.get_Dr();
			const real_matrix_type& Ds = triangleNodesProvisioner.get_Ds();

			const real_matrix_type& rx = triangleNodesProvisioner.get_rx();
			const real_matrix_type& ry = triangleNodesProvisioner.get_ry();
			const real_matrix_type& sx = triangleNodesProvisioner.get_sx();
			const real_matrix_type& sy = triangleNodesProvisioner.get_sy();

			const real_matrix_type &nx = triangleNodesProvisioner.get_nx();
			const real_matrix_type &ny = triangleNodesProvisioner.get_ny();

			// Get volume to surface maps.
			const index_vector_type& vmapM = triangleNodesProvisioner.get_vmapM();
			const index_vector_type& vmapP = triangleNodesProvisioner.get_vmapP();

			// boundary indices.
			const index_hashmap bcHash = triangleNodesProvisioner.get_bcMap();
			const std::vector<index_type> mapW = bcHash.at(3);

			index_type numFaces = triangleNodesProvisioner.NumFaces;
			index_type Nfp = triangleNodesProvisioner.get_NumFacePoints();
			index_type Np = triangleNodesProvisioner.get_NumLocalPoints();
			index_type K = triangleNodesProvisioner.get_NumElements();

			index_type numFaceNodes = numFaces*Nfp*K;
			index_type numTotalNodes = Np*K;

			real_vector_type dh(numFaceNodes), dhu(numFaceNodes), dhv(numFaceNodes);

			real_vector_type hM(numFaceNodes), hP(numFaceNodes);
			real_vector_type huM(numFaceNodes), huP(numFaceNodes);
			real_vector_type hvM(numFaceNodes), hvP(numFaceNodes);

			real_vector_type nxVec(numFaceNodes), nyVec(numFaceNodes);

			real_vector_type hVec(numTotalNodes), huVec(numTotalNodes), hvVec(numTotalNodes);

			dh = 0.*ii;
			dhu = 0.*ii;
			dhv = 0.*ii;

			huM = 0.*ii;
			huP = 0.*ii;
			nxVec = 0.*ii;
			nyVec = 0.*ii;

			// We want to apply maps to column-wise ordering of the nodes.
			const bool byRowsOpt = false;

			fullToVector(nx, nxVec, byRowsOpt);
			fullToVector(ny, nyVec, byRowsOpt);

			fullToVector(h, hVec, byRowsOpt);
			fullToVector(hu, huVec, byRowsOpt);
			fullToVector(hv, hvVec, byRowsOpt);

			applyIndexMap(hVec, vmapM, hM);
			applyIndexMap(hVec, vmapP, hP);

			applyIndexMap(huVec, vmapM, huM);
			applyIndexMap(huVec, vmapP, huP);

			applyIndexMap(hvVec, vmapM, hvM);
			applyIndexMap(hvVec, vmapP, hvP);

			// BC's - no flow through walls.
			for (index_type i=0; i < static_cast<index_type>(mapW.size()); ++i) {
				index_type w = mapW[i];
				huP(w) = huM(w) - 2*(huM(w)*nx(w) + hvM(w)*ny(w));
				hvP(w) = hvM(w) - 2*(huM(w)*nx(w) + hvM(w)*ny(w));
			}

			// Compute jump in state vector
			dh = hM - hP;
			dhu = huM - huP;
			dhv = hvM - hvP;

			// calculate flux tensor: '-' traces.
			real_vector_type& F1M = huM, G1M = hvM;
			real_vector_type F2M(numFaceNodes), G2M(numFaceNodes), G3M(numFaceNodes);

			F2M = (huM*huM)/hM + 0.5*g*hM*hM; G2M = (huM*hvM)/hM;
			real_vector_type& F3M = G2M;  G3M = (hvM*hvM)/hM + 0.5*g*hM*hM;

			// '+' traces:
			real_vector_type& F1P = huP, G1P = hvP;
			real_vector_type F2P(numFaceNodes), G2P(numFaceNodes), G3P(numFaceNodes);

			F2P = (huP*huP)/hP + 0.5*g*hP*hP; G2P = (huP*hvP)/hP;
			real_vector_type& F3P = G2P;  G3P = (hvP*hvP)/hP + 0.5*g*hM*hM;

			// Full fields
			real_matrix_type& F1 = hu, G1 = hv;
			real_matrix_type F2(numFaceNodes), G2(numFaceNodes), G3(numFaceNodes);

			F2 = (hu*hu)/h + 0.5*g*h*h; G2 = (hu*hv)/h;
			real_matrix_type& F3 = G2;  G3 = (hv*hv)/h + 0.5*g*h*h;

			// compute maximum linearized wave-speed - for numerical flux stabilization.
			real_vector_type uM(numFaceNodes), vM(numFaceNodes);
			real_vector_type uP(numFaceNodes), vP(numFaceNodes);

			uM = huM/hM; vM= hvM/hM;

			real_vector_type spdM(numFaceNodes), spdM2(numFaceNodes);
			real_vector_type spdP(numFaceNodes), spdP2(numFaceNodes);
			real_vector_type spdMax(numFaceNodes);

			spdM = sqrt(uM*uM + vM*vM) + sqrt(g*hM);
			spdP = sqrt(uP*uP + vP*vP) + sqrt(g*hP);

			// Compute 'trace maximum' over '-' and '+'.
			for (index_type i=0; i < numFaceNodes; ++i) {
				spdMax(i) = max(spdM(i), spdP(i));
			}

			real_vector_type dFlux1(numFaceNodes), dFlux2(numFaceNodes), dFlux3(numFaceNodes);

			// strong form: Compute flux jump vector. (fluxM - numericalFlux ) dot n
			dFlux1 = 0.5*((F1M - F1P)*nxVec + (G1M-G1P)*nyVec - spdMax*dh);
			dFlux2 = 0.5*((F2M - F2P)*nxVec + (G2M-G2P)*nyVec - spdMax*dhu);
			dFlux3 = 0.5*((F3M - F3P)*nxVec + (G2M-G2P)*nyVec - spdMax*dhv);

			real_matrix_type dFlux1Mat(Nfp*numFaces, K), dFlux2Mat(Nfp*numFaces, K), dFlux3Mat(Nfp*numFaces, K);

			vectorToFull(dFlux1, dFlux1Mat, byRowsOpt);
			vectorToFull(dFlux2, dFlux2Mat, byRowsOpt);
			vectorToFull(dFlux3, dFlux3Mat, byRowsOpt);

			// Assumes PDE has been left-multiplied by local inverse mass matrix, so all we have left
			// is the differentiation matrix contribution, and the surface integral.
			// RHS == -Flux divergence + surface integral contributions.

			// Flux divergence:
			RHS1 = -(rx*sum(Dr(ii,kk)*F1(kk,jj), kk) + sx*sum(Ds(ii,kk)*F1(kk,jj), kk));
			RHS1+= -(ry*sum(Dr(ii,kk)*F1(kk,jj), kk) + sy*sum(Ds(ii,kk)*F1(kk,jj), kk));

			RHS2 = -(rx*sum(Dr(ii,kk)*F2(kk,jj), kk) + sx*sum(Ds(ii,kk)*F2(kk,jj), kk));
			RHS2+= -(ry*sum(Dr(ii,kk)*F2(kk,jj), kk) + sy*sum(Ds(ii,kk)*F2(kk,jj), kk));

			RHS3 = -(rx*sum(Dr(ii,kk)*F3(kk,jj), kk) + sx*sum(Ds(ii,kk)*F3(kk,jj), kk));
			RHS3+= -(ry*sum(Dr(ii,kk)*F3(kk,jj), kk) + sy*sum(Ds(ii,kk)*F3(kk,jj), kk));

			// Surface integral contributions
			const real_matrix_type& Fscale = triangleNodesProvisioner.get_Fscale();
			const real_matrix_type& Lift   = triangleNodesProvisioner.get_Lift();

			real_matrix_type surfaceRHS1(Nfp*numFaces, K), surfaceRHS2(Nfp*numFaces, K), surfaceRHS3(Nfp*numFaces, K);

			// Apply Jacobian scaling.
			surfaceRHS1 = Fscale*dFlux1Mat;
			surfaceRHS2 = Fscale*dFlux2Mat;
			surfaceRHS3 = Fscale*dFlux3Mat;

			// Integrate using Lifting operator, add to RHS
			RHS1+= sum(Lift(ii,kk)*surfaceRHS1(kk,jj), kk);
			RHS2+= sum(Lift(ii,kk)*surfaceRHS2(kk,jj), kk);
			RHS3+= sum(Lift(ii,kk)*surfaceRHS3(kk,jj), kk);
		} // computeRHS
	} // namespace advec1d
} // namespace blitzdg
