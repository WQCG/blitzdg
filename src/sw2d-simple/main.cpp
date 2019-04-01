
// Copyright (C) 2017-2019  Waterloo Quantitative Consulting Group, Inc.
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
#include "VtkOutputter.hpp"
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
	const real_type finalTime = 1000.0;
	real_type t = 0.0;

	// Numerical parameters (N = Order of polynomials)
	const index_type N = 4;
	const real_type CFL = 0.65;

	// Build dependencies.
	MeshManager meshManager;
	meshManager.readMesh("input/figure8_mesh.msh");
	const index_type K = meshManager.get_NumElements();

	// Dependency-inject mesh manager to nodes provisioner.
	TriangleNodesProvisioner triangleNodesProvisioner(N, meshManager);

	// Pre-processing step - build polynomial dealiasing filter.
	triangleNodesProvisioner.buildFilter(0.9*N, N);

#ifndef __MINGW32__
	VtkOutputter outputter(triangleNodesProvisioner);
#else
	CsvOutputter outputter;
#endif

	const real_matrix_type& x = triangleNodesProvisioner.get_xGrid();
	const real_matrix_type& y = triangleNodesProvisioner.get_yGrid();
	index_type Np = triangleNodesProvisioner.get_NumLocalPoints();

	const real_matrix_type& Filt = triangleNodesProvisioner.get_Filter();

	real_matrix_type H(Np,K);
	real_matrix_type eta(Np, K);
	real_matrix_type h(Np, K);
	real_matrix_type u(Np, K);
	real_matrix_type v(Np, K);
	real_matrix_type hu(Np, K);
	real_matrix_type hv(Np, K);
	real_matrix_type spd(Np, K);
	real_matrix_type RHS1(Np, K), RHS2(Np, K), RHS3(Np, K);
	real_matrix_type resRK1(Np, K), resRK2(Np, K), resRK3(Np, K);

	firstIndex ii;
	secondIndex jj;
	thirdIndex kk;

	// Intialize fields.
	H = 10.0 + 0*jj;
	eta = -1.*(x/1500.0);
	//eta = exp(-(x/200)*(x/200) - (y/200)*(y/200));

	u = 0*jj;
	v = 0*jj;
	h = H + eta;
	hu = h*u;
	hv = h*v;

	real_matrix_type Fscale = triangleNodesProvisioner.get_Fscale();

	spd = blitz::sqrt(u*u + v*v) + blitz::sqrt(g*h);
	const index_vector_type& vmapM = triangleNodesProvisioner.get_vmapM();

	real_vector_type spdM(triangleNodesProvisioner.NumFaces*triangleNodesProvisioner.get_NumFacePoints()*K);
	real_vector_type spdVec(Np*K), fsVec(triangleNodesProvisioner.NumFaces*triangleNodesProvisioner.get_NumFacePoints()*K);

	fullToVector(spd, spdVec, false);
	fullToVector(Fscale, fsVec, false);
	applyIndexMap(spdVec, vmapM, spdM);

	real_type Fsc_max = blitz::max(abs(fsVec)*spdM);
	real_type dt = CFL/((N+1)*(N+1)*0.5*Fsc_max);

	RHS1 = 0*jj;
	RHS2 = 0*jj;
	RHS3 = 0*jj;

	resRK1= 0*jj;
	resRK2= 0*jj;
	resRK3= 0*jj;

	index_type count = 0;

	while (t < finalTime) {
		if ((count % 10) == 0) {
			eta = h-H;
			cout << "t=" << t << ", eta_max=" << max(eta) << ", dt=" << dt << "\n";
			string fileName = outputter.generateFileName("eta", count);

			std::map<std::string, real_matrix_type> fields;
			fields.insert({"eta", eta});
			outputter.writeFieldsToFiles(fields, count);
		}

		// SSP RK2
		sw2d::computeRHS(h, hu, hv, g, triangleNodesProvisioner, RHS1, RHS2, RHS3);
		//RHS1 = sum(Filt(ii,kk)*RHS1(kk,jj), kk);
		//RHS2 = sum(Filt(ii,kk)*RHS2(kk,jj), kk);
		//RHS3 = sum(Filt(ii,kk)*RHS3(kk,jj), kk);

		real_matrix_type h1(Np,K), hu1(Np,K), hv1(Np,K), u1(Np,K), v1(Np,K);

		h1 =  h + 0.5*dt*RHS1;
		hu1 = hu + 0.5*dt*RHS2;
		hv1 = hv + 0.5*dt*RHS3;

		sw2d::computeRHS(h1, hu1, hv1, g, triangleNodesProvisioner, RHS1, RHS2, RHS3);
		//RHS1 = sum(Filt(ii,kk)*RHS1(kk,jj), kk);
		//RHS2 = sum(Filt(ii,kk)*RHS2(kk,jj), kk);
		//RHS3 = sum(Filt(ii,kk)*RHS3(kk,jj), kk);

		h  += dt*RHS1;
		hu += dt*RHS2;
		hv += dt*RHS3;

		eta = h-H;


		real_type eta_max = normMax(eta);
		if ( std::abs(eta_max) > 1e8  || std::isnan(eta_max) )
			throw std::runtime_error("A numerical instability has occurred!");

		u = hu / h;
		v = hv / h;

		spd = blitz::sqrt(u*u + v*v) + blitz::sqrt(g*h);
		fullToVector(spd, spdVec, false);
		applyIndexMap(spdVec, vmapM, spdM);

		real_type Fsc_max = max(abs(fsVec)*spdM);
		dt = CFL/((N+1)*(N+1)*0.5*Fsc_max);

		

		t += dt;
		++count;

	}
	real_matrix_type etafinal(Np, K);
	etafinal = eta;

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
				huP(w) = huM(w) - 2*nxVec(w)*(huM(w)*nxVec(w) + hvM(w)*nyVec(w));
				hvP(w) = hvM(w) - 2*nyVec(w)*(huM(w)*nxVec(w) + hvM(w)*nyVec(w));
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
			real_vector_type& F3P = G2P;  G3P = (hvP*hvP)/hP + 0.5*g*hP*hP;

			// Full fields
			real_matrix_type& F1 = hu, G1 = hv;
			real_matrix_type F2(Np, K), G2(Np, K), G3(Np, K);

			F2 = (hu*hu)/h + 0.5*g*h*h; G2 = (hu*hv)/h;
			real_matrix_type& F3 = G2;  G3 = (hv*hv)/h + 0.5*g*h*h;

			// compute maximum linearized wave-speed - for numerical flux stabilization.
			real_vector_type uM(numFaceNodes), vM(numFaceNodes);
			real_vector_type uP(numFaceNodes), vP(numFaceNodes);

			uM = huM/hM; vM = hvM/hM;
			uP = huP/hP; vP = hvP/hP;

			real_vector_type spdM(numFaceNodes), spdP(numFaceNodes);
			real_vector_type spdMax(numFaceNodes);

			spdM = sqrt(uM*uM + vM*vM) + sqrt(g*hM);
			spdP = sqrt(uP*uP + vP*vP) + sqrt(g*hP);

			// Compute 'trace max' over '-' and '+'.
			for (index_type i=0; i < numFaceNodes; ++i)
				spdMax(i) = max(spdM(i), spdP(i));

			// Compute max over each face to get maximum linearized eigenvalue.
			real_matrix_type lambda(Nfp, numFaces*K);
			vectorToFull(spdMax, lambda, false);
		
			real_matrix_type lambdaMaxMat(Nfp, numFaces*K);
			lambdaMaxMat = (1 + 0*ii)*(blitz::max(lambda(kk, jj), kk));

			fullToVector(lambdaMaxMat, spdMax, byRowsOpt);

			real_vector_type dFlux1(numFaceNodes), dFlux2(numFaceNodes), dFlux3(numFaceNodes);

			// strong form: Compute flux jump vector. (fluxM - numericalFlux ) dot n
			dFlux1 = 0.5*((F1M - F1P)*nxVec + (G1M-G1P)*nyVec - spdMax*dh);
			dFlux2 = 0.5*((F2M - F2P)*nxVec + (G2M-G2P)*nyVec - spdMax*dhu);
			dFlux3 = 0.5*((F3M - F3P)*nxVec + (G3M-G3P)*nyVec - spdMax*dhv);

			real_matrix_type dFlux1Mat(Nfp*numFaces, K), dFlux2Mat(Nfp*numFaces, K), dFlux3Mat(Nfp*numFaces, K);

			vectorToFull(dFlux1, dFlux1Mat, byRowsOpt);
			vectorToFull(dFlux2, dFlux2Mat, byRowsOpt);
			vectorToFull(dFlux3, dFlux3Mat, byRowsOpt);

			// Assumes PDE has been left-multiplied by local inverse mass matrix, so all we have left
			// is the differentiation matrix contribution, and the surface integral.
			// RHS == -Flux divergence + surface integral contributions.

			// Flux divergence:
			RHS1 = -(rx*sum(Dr(ii,kk)*F1(kk,jj), kk) + sx*sum(Ds(ii,kk)*F1(kk,jj), kk));
			RHS1+= -(ry*sum(Dr(ii,kk)*G1(kk,jj), kk) + sy*sum(Ds(ii,kk)*G1(kk,jj), kk));

			RHS2 = -(rx*sum(Dr(ii,kk)*F2(kk,jj), kk) + sx*sum(Ds(ii,kk)*F2(kk,jj), kk));
			RHS2+= -(ry*sum(Dr(ii,kk)*G2(kk,jj), kk) + sy*sum(Ds(ii,kk)*G2(kk,jj), kk));

			RHS3 = -(rx*sum(Dr(ii,kk)*F3(kk,jj), kk) + sx*sum(Ds(ii,kk)*F3(kk,jj), kk));
			RHS3+= -(ry*sum(Dr(ii,kk)*G3(kk,jj), kk) + sy*sum(Ds(ii,kk)*G3(kk,jj), kk));

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
