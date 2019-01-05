
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
#include "Constants.hpp"
#include <blitz/array.h>
#include <math.h>
#include <string>
#include <stdexcept>
#include <unordered_set>

using blitz::firstIndex;
using blitz::secondIndex;
using blitz::thirdIndex;
using blitz::sum;
using blitz::sqrt;
using blitzdg::constants::pi;
using std::abs;
using std::string;
using std::cout;
using std::endl;
using std::sqrt;
using std::max;
using std::unordered_set;
using std::cos;

int main(int argc, char **argv) {
	using namespace blitzdg;

	printDisclaimer();

	// Physical parameters
	const real_type g = 9.81;
	const real_type CD = 2.5e-3;
	const real_type f = 1.0070e-4;

	const real_type finalTime = 24.0*3600;
	real_type t = 0.0;

	// Numerical parameters (N = Order of polynomials)
	const index_type N = 1;
	const real_type CFL = 0.85;

	// Build dependencies.
	MeshManager meshManager;
	meshManager.readVertices("input/vh_verts_z.dat");
	meshManager.readElements("input/vh_els_0.oct");
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

	real_matrix_type H(Np,K), eta(Np, K), h(Np, K), u(Np, K), v(Np, K), hu(Np, K), hv(Np, K);
	real_matrix_type RHS1(Np, K), RHS2(Np, K), RHS3(Np, K);
	real_matrix_type resRK1(Np, K), resRK2(Np, K), resRK3(Np, K);

	unordered_set<index_type> obcNodes = {0,1,2,3,5,6,8,11,12,15,19,20,25,28,31,37,39,44,51,52,60,65,69,78,81,88,98,99,110,117,122,137,157,178,199,220,242,264,287,310,334,358,383,409,435,461,488,516,544,573,602,631,660,690,721,753,784,816,849,933,1018,1106,1193,1281,1369,1455,1538,1619,1699,1781,1813,1846,1879,1912,1946,1980,2014,2049,2084,2120,2156,2192,2228,2265,2302,2340,2374,2412,2448,2488,2526,2565,2603,2642,2682,2723,2765,2804,2845,2886,2929,2972,3014,3054,3097,3144,3190,3236,3284,3329,3376,3492,3607,3724,3835,3953,4074,4192,4321,4461,4625,4691,4759,4829,4900,4971,5039,5107,5175,5253,5319,5383,5440,5494,5562,5620,5681,5746,5904,6064,6221,6374,6541,6704,6872,7028,7186,7374,7580,7788,8009,8233,8482,8783,9138,9587,9762,9920,10068,10202,10338,10461,10582,10699,10817,10918,11015,11118,11217,11310,11564,11829,12073,12326,12584,12866,13191,13559,14017,14554,15146,15717,16037,16309,16533,16750,16878,17009,17147,17280,17554,17870,18093,18276,18471,18712,18940,19142,19294,19424,19540,19652,19740,19836,19929,20013,20111,20216,20341,20465,20605,20733,20852,20970,21066,21157,21246,21345,21469,21591,21726,21847,21952,22043,22127,22195,22257,22308,22351,22391,22429,22467,22500,22529,22557,22586,22617};

	firstIndex ii;
	secondIndex jj;

	// Intialize fields.
	H = 10.0 + 0*jj;
	eta = 0*jj;
	u = 0*jj;
	v = 0*jj;
	h = H + eta;
	hu = h*u;
	hv = h*v;

	const real_matrix_type& Fscale = triangleNodesProvisioner.get_Fscale();
	real_type Fsc_max = max(abs(Fscale));
		
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

	const index_vector_type& EToV = meshManager.get_Elements();

	// Make copy of bcMap, for hacking it.
    index_vector_type bcType = meshManager.get_BCType();

	// deal with open boundary conditions.
	for (index_type k=0; k < K; ++k) {
			index_type v1 = EToV(3*k);
			index_type v2 = EToV(3*k+1);
			index_type v3 = EToV(3*k+2);

			if (obcNodes.count(v1) > 0 && obcNodes.count(v2) > 0)
				bcType(3*k) = 2;
			if (obcNodes.count(v2) > 0 && obcNodes.count(v3) > 0)
				bcType(3*k+1) = 2;
			if (obcNodes.count(v3) > 0 && obcNodes.count(v1) > 0)
				bcType(3*k+2) = 2;
	}

	triangleNodesProvisioner.buildBCHash(bcType);

	real_type dt;
	while (t < finalTime) {
		real_matrix_type u(Np,K), v(Np,K);
		u = hu/h; v = hv/h;
		real_type spdmax = blitz::max(blitz::sqrt(u*u+v*v) + sqrt(g*h));
		dt = CFL/((N+1)*(N+1)*0.5*Fsc_max*spdmax);

		if ((count % 50) == 0) {
			cout << "dt=" << dt << endl;
			eta = h-H;
			u = hu/h;
			v = hv/h;
			cout << "t=" << t << ", eta_max=" << max(eta) << endl;
			string fileName = outputter.generateFileName("eta", count);
			outputter.writeFieldToFile(fileName, eta, delim);
		}	

		for (index_type i=0; i < LSERK4::numStages;  i++ ) {

			// Calculate Right-hand side.
			sw2d::computeRHS(h, hu, hv, g, H, CD, f, triangleNodesProvisioner, RHS1, RHS2, RHS3, t);

			// Compute Runge-Kutta Residual.
			resRK1 = LSERK4::rk4a[i]*resRK1 + dt*RHS1;
			resRK2 = LSERK4::rk4a[i]*resRK2 + dt*RHS2;
			resRK3 = LSERK4::rk4a[i]*resRK3 + dt*RHS3;

			// Update solution.
			h  += LSERK4::rk4b[i]*resRK1;
			hu += LSERK4::rk4b[i]*resRK2;
			hv += LSERK4::rk4b[i]*resRK3;
		}

		real_type eta_max = normMax(eta);
		if ( std::abs(eta_max) > 1e8  || std::isnan(eta_max) )
			throw std::runtime_error("A numerical instability has occurred!");

		t += dt;
		count++;
	}
	real_matrix_type etafinal(Np, K);
	etafinal = h - H;

	return 0;
} // end main

namespace blitzdg {
	namespace sw2d {
		void computeRHS(real_matrix_type& h, real_matrix_type& hu, real_matrix_type& hv, real_type g, real_matrix_type& H, real_type CD, real_type f, TriangleNodesProvisioner& triangleNodesProvisioner, real_matrix_type& RHS1, real_matrix_type& RHS2, real_matrix_type& RHS3, real_type t) {

			real_type T_tide = 3600*8; // or something?
			real_type om_tide = 2.0*pi/T_tide;
			real_type amp_tide = 1.0;

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
			const index_hashmap& bcHash = triangleNodesProvisioner.get_bcMap();
			const std::vector<index_type>& mapW = bcHash.at(3);
			const std::vector<index_type>& mapO = bcHash.at(2);

			index_type numFaces = triangleNodesProvisioner.NumFaces;
			index_type Nfp = triangleNodesProvisioner.get_NumFacePoints();
			index_type Np = triangleNodesProvisioner.get_NumLocalPoints();
			index_type K = triangleNodesProvisioner.get_NumElements();

			index_type numFaceNodes = numFaces*Nfp*K;
			index_type numTotalNodes = Np*K;

			real_vector_type dh(numFaceNodes), dhu(numFaceNodes), dhv(numFaceNodes);

			real_vector_type hM(numFaceNodes), hP(numFaceNodes), HM(numFaceNodes), HP(numFaceNodes);
			real_vector_type huM(numFaceNodes), huP(numFaceNodes);
			real_vector_type hvM(numFaceNodes), hvP(numFaceNodes);

			real_vector_type nxVec(numFaceNodes), nyVec(numFaceNodes);

			real_vector_type hVec(numTotalNodes), huVec(numTotalNodes), hvVec(numTotalNodes), Hvec(numTotalNodes);

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

			fullToVector(H, Hvec, byRowsOpt);

			applyIndexMap(hVec, vmapM, hM);
			applyIndexMap(hVec, vmapP, hP);

			applyIndexMap(huVec, vmapM, huM);
			applyIndexMap(huVec, vmapP, huP);

			applyIndexMap(hvVec, vmapM, hvM);
			applyIndexMap(hvVec, vmapP, hvP);

			applyIndexMap(Hvec, vmapM, HM);
			applyIndexMap(Hvec, vmapP, HP);

			// BC's - no flow through walls.
			for (index_type i=0; i < static_cast<index_type>(mapW.size()); ++i) {
				index_type w = mapW[i];
				huP(w) = huM(w) - 2*nxVec(w)*(huM(w)*nxVec(w) + hvM(w)*nyVec(w));
				hvP(w) = hvM(w) - 2*nyVec(w)*(huM(w)*nxVec(w) + hvM(w)*nyVec(w));
			}

			// OBC's - free surface moves up and down according to the tidal forcing.
			for (index_type i=0; i < static_cast<index_type>(mapO.size()); ++i) {
				index_type o = mapO[i];
				huP(o) = huM(o);
				hvP(o) = hvM(o);
				hP(o) = HM(o) + amp_tide*std::cos(om_tide*t);
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

			// Add source terms
			real_matrix_type u(Np, K), v(Np, K); 
			u = hu/h;
			v = hv/h;

			// bottom drag
			RHS2+= -CD*u*blitz::abs(u);
			RHS3+= -CD*v*blitz::abs(v);

			// Coriolis force
			RHS2+=  f*hv;
			RHS3+= -f*hu;
		} // computeRHS
	} // namespace advec1d
} // namespace blitzdg
