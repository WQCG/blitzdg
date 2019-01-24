
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
#include "CSVFileReader.hpp"
#include "BCtypes.hpp"
#ifndef __MINGW32__
	#include "VtkOutputter.hpp"
#endif
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
using std::hypot;

int main(int argc, char **argv) {
	using namespace blitzdg;

	printDisclaimer();

	const index_type rkStages = 2;

	sw2d::physParams p; 
	p.CD = 2.5e-3;
	p.f = 1.0070e-4;
	p.initTime = 0.0;
	p.finalTime = 24.0*3600;

	sw2d::numParams n;
	n.N = 1;
	n.CFL = 0.6;
	n.outputInterval = 20;
	n.filterPercent = 0.95;
	n.filterOrder = 4;

	sw2d::fields fds[2];

	real_type t = p.initTime;

	// Build dependencies.
	MeshManager meshManager;
	meshManager.readVertices("input/vh_verts_z.dat");
	meshManager.readElements("input/vh_els_0.oct");
	// Dependency-inject mesh manager to nodes provisioner.
	TriangleNodesProvisioner triangleNodesProvisioner(n.N, meshManager);

	// Pre-processing step - build polynomial dealiasing filter.
	triangleNodesProvisioner.buildFilter(n.filterPercent*n.N, n.filterOrder);

	DGContext2D dg = triangleNodesProvisioner.get_DGContext();


#ifndef __MINGW32__
	VtkOutputter vtkOutputter(triangleNodesProvisioner);
#endif
	
	CsvOutputter outputter;

	// Allocate memory for fields.
	for (index_type i=0; i < rkStages; ++i) {
		fds[i].h   = real_matrix_type(dg.Np, dg.K);
		fds[i].hu  = real_matrix_type(dg.Np, dg.K);
		fds[i].hv  = real_matrix_type(dg.Np, dg.K);
		fds[i].eta = real_matrix_type(dg.Np, dg.K);
		fds[i].Hx  = real_matrix_type(dg.Np, dg.K);
		fds[i].Hy  = real_matrix_type(dg.Np, dg.K);
		fds[i].u   = real_matrix_type(dg.Np, dg.K);
		fds[i].v   = real_matrix_type(dg.Np, dg.K);
		fds[i].RHS1 = real_matrix_type(dg.Np, dg.K);
		fds[i].RHS2 = real_matrix_type(dg.Np, dg.K);
		fds[i].RHS3 = real_matrix_type(dg.Np, dg.K);
		fds[i].resRK1 = real_matrix_type(dg.Np, dg.K);
		fds[i].resRK2 = real_matrix_type(dg.Np, dg.K);
	}

	// Make amazing references to the fields at different time-levels.
	sw2d::fields& fields_n = fds[0];
	sw2d::fields& fields_np1 = fds[1];

	unordered_set<index_type> obcNodes = {0,1,2,3,5,6,8,11,12,15,19,20,25,28,31,37,39,44,51,52,60,65,69,78,81,88,98,99,110,117,122,137,157,178,199,220,242,264,287,310,334,358,383,409,435,461,488,516,544,573,602,631,660,690,721,753,784,816,849,933,1018,1106,1193,1281,1369,1455,1538,1619,1699,1781,1813,1846,1879,1912,1946,1980,2014,2049,2084,2120,2156,2192,2228,2265,2302,2340,2374,2412,2448,2488,2526,2565,2603,2642,2682,2723,2765,2804,2845,2886,2929,2972,3014,3054,3097,3144,3190,3236,3284,3329,3376,3492,3607,3724,3835,3953,4074,4192,4321,4461,4625,4691,4759,4829,4900,4971,5039,5107,5175,5253,5319,5383,5440,5494,5562,5620,5681,5746,5904,6064,6221,6374,6541,6704,6872,7028,7186,7374,7580,7788,8009,8233,8482,8783,9138,9587,9762,9920,10068,10202,10338,10461,10582,10699,10817,10918,11015,11118,11217,11310,11564,11829,12073,12326,12584,12866,13191,13559,14017,14554,15146,15717,16037,16309,16533,16750,16878,17009,17147,17280,17554,17870,18093,18276,18471,18712,18940,19142,19294,19424,19540,19652,19740,19836,19929,20013,20111,20216,20341,20465,20605,20733,20852,20970,21066,21157,21246,21345,21469,21591,21726,21847,21952,22043,22127,22195,22257,22308,22351,22391,22429,22467,22500,22529,22557,22586,22617};

	firstIndex ii;
	secondIndex jj;
	thirdIndex kk;

	const string depthFile = "input/H0_try2.oct";
	sw2d::readDepthData(depthFile, fields_n.H);

#ifndef __MINGW32__
	string vtkFileName = "H.vtu";
	vtkOutputter.writeFieldToFile(vtkFileName, fields_n.H, "H");
#endif

	// Set initial values for fields
	fields_n.eta = 0*jj;
	fields_n.h = fields_n.H + fields_n.eta;
	fields_n.u = 0*jj;
	fields_n.v = 0*jj;
	fields_n.hu = fields_n.h*fields_n.u;
	fields_n.hv = fields_n.h*fields_n.v;
		
	fields_n.RHS1 = 0*jj;
	fields_n.RHS2 = 0*jj;
	fields_n.RHS3 = 0*jj;

	fields_n.resRK1= 0*jj;
	fields_n.resRK2= 0*jj;
	fields_n.resRK3= 0*jj;

	index_type count = 0;

	const char delim = ' ';
	outputter.writeFieldToFile("x.dat", dg.x, delim);
	outputter.writeFieldToFile("y.dat", dg.y, delim);

	outputter.writeFieldToFile("H.dat", fields_n.H, delim);

	const index_vector_type& EToV = meshManager.get_Elements();

	// Make copy of bcMap, for hacking it.
    index_vector_type bcType = meshManager.get_BCType();

	// Get bed slopes
	fields_n.Hx = (dg.rx*sum(dg.Dr(ii,kk)*fields_n.H(kk,jj), kk) + dg.sx*sum(dg.Ds(ii,kk)*fields_n.H(kk,jj), kk));
	fields_n.Hy = (dg.ry*sum(dg.Dr(ii,kk)*fields_n.H(kk,jj), kk) + dg.sy*sum(dg.Ds(ii,kk)*fields_n.H(kk,jj), kk));

	fields_n.Hx = sum(dg.Filt(ii,kk)*fields_n.Hx(kk,jj), kk);
	fields_n.Hy = sum(dg.Filt(ii,kk)*fields_n.Hy(kk,jj), kk);

	// deal with open boundary conditions.
	for (index_type k=0; k < dg.K; ++k) {
			index_type v1 = EToV(3*k);
			index_type v2 = EToV(3*k+1);
			index_type v3 = EToV(3*k+2);

			if (obcNodes.count(v1) > 0 && obcNodes.count(v2) > 0)
				bcType(3*k) = BCTag::Out;
			if (obcNodes.count(v2) > 0 && obcNodes.count(v3) > 0)
				bcType(3*k+1) = BCTag::Out;
			if (obcNodes.count(v3) > 0 && obcNodes.count(v1) > 0)
				bcType(3*k+2) = BCTag::Out;
	}

	triangleNodesProvisioner.buildBCHash(bcType);



	// sponge layer -- build that wall!
	real_matrix_type spongeCoeff(dg.Np, dg.K);
	spongeCoeff = 0.0*jj;
	sw2d::buildSpongeCoeff(dg, 10000.0, 1000.0, spongeCoeff);

	outputter.writeFieldToFile("sponge.dat", spongeCoeff, delim);

	index_type numFaceNodes = dg.NumFaces*dg.Nfp*dg.K;
	real_type dt;
	while (t < p.finalTime) {
		real_vector_type uVec(dg.Np*dg.K), vVec(dg.Np*dg.K), hVec(dg.Np*dg.K);
		real_vector_type uM(numFaceNodes), vM(numFaceNodes), hM(numFaceNodes), fsVec(numFaceNodes);

		fields_n.u = fields_n.hu/fields_n.h; 
		fields_n.v = fields_n.hv/fields_n.h;

		fullToVector(fields_n.u, uVec, false);
		fullToVector(fields_n.v, vVec, false);
		fullToVector(fields_n.h, hVec, false);
		fullToVector(dg.Fscale , fsVec, false);
		applyIndexMap(uVec, dg.vmapM, uM);
		applyIndexMap(vVec, dg.vmapM, vM);
		applyIndexMap(hVec, dg.vmapM, hM);
		real_vector_type spd(numFaceNodes);
		spd = blitz::sqrt(uM*uM + vM*vM) + blitz::sqrt(p.g*hM);
		real_type spdFscaleMax = blitz::max(fsVec*spd);
		dt = n.CFL/((n.N+1)*(n.N+1)*0.5*spdFscaleMax);

		if ((count % n.outputInterval) == 0) {
			cout << "dt=" << dt << endl;
			cout << "t=" << t << ", h_min=" << blitz::min(fields_n.h) << ", h_max=" << normMax(fields_n.h) << ", hu_max=" << normMax(fields_n.hu) << ", hv_max=" << normMax(fields_n.hv) << endl;
#ifndef __MINGW32__
			string fileName = vtkOutputter.generateFileName("eta", count);
			vtkOutputter.writeFieldToFile(fileName, fields_n.eta, "eta");
#endif
		}	

		// 2nd order SSP Runge-Kutta
		sw2d::computeRHS(fields_n, n, p, dg, t);
		// Update solution.
		fields_np1.h  = fields_n.h  + dt*fields_n.RHS1;
		fields_np1.hu = fields_n.hu + dt*fields_n.RHS2;
		fields_np1.hv = fields_n.hv + dt*fields_n.RHS3;

		// sponge layer relaxation -- to control insane velocities near the open boundary.		
		fields_np1.hu /= (1.0 + spongeCoeff*fields_np1.hu*fields_np1.hu);
		fields_np1.hv /= (1.0 + spongeCoeff*fields_np1.hv*fields_np1.hv);

		sw2d::computeRHS(fields_np1, n, p, dg, t);

		fields_n.h  = 0.5*(fields_n.h  + fields_np1.h  + dt*fields_np1.RHS1);
		fields_n.hu = 0.5*(fields_n.hu + fields_np1.hu + dt*fields_np1.RHS2);
		fields_n.hv = 0.5*(fields_n.hv + fields_np1.hv + dt*fields_np1.RHS3);

		fields_n.hu /= (1.0 + spongeCoeff*fields_n.hu*fields_n.hu);
		fields_n.hv /= (1.0 + spongeCoeff*fields_n.hv*fields_n.hv);

		fields_n.eta = fields_n.h-fields_n.H;
		real_type eta_max = normMax(fields_n.eta);
		if ( std::abs(eta_max) > 1e8  || std::isnan(sum(fields_n.eta)))
			throw std::runtime_error("A numerical instability has occurred!");

		t += dt;
		count++;
	}
	real_matrix_type etafinal(dg.Np, dg.K);
	etafinal = fields_n.h - fields_n.H;

	return 0;
} // end main

namespace blitzdg {
	namespace sw2d {
		void computeRHS(fields fds, const numParams& num, const physParams& phys, const DGContext2D& dg, real_type t) {

			real_type T_tide = 3600*12.42; // or something?
			real_type om_tide = 2.0*pi/T_tide;
			real_type amp_tide = 0*3.e-4; // 3.
			real_type g = phys.g;

			// Blitz indices
			firstIndex ii;
			secondIndex jj;
			thirdIndex kk;

			// boundary indices.
			const index_hashmap& bcHash = dg.bcHash;
			const std::vector<index_type>& mapW = bcHash.at(BCTag::Wall);
			const std::vector<index_type>& mapO = bcHash.at(BCTag::Out);

			index_type numFaces = dg.NumFaces;
			index_type Nfp = dg.Nfp;
			index_type Np = dg.Np;
			index_type K = dg.K;

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

			fullToVector(dg.nx, nxVec, byRowsOpt);
			fullToVector(dg.ny, nyVec, byRowsOpt);

			fullToVector(fds.h, hVec, byRowsOpt);
			fullToVector(fds.hu, huVec, byRowsOpt);
			fullToVector(fds.hv, hvVec, byRowsOpt);

			fullToVector(fds.H, Hvec, byRowsOpt);

			applyIndexMap(hVec, dg.vmapM, hM);
			applyIndexMap(hVec, dg.vmapP, hP);

			applyIndexMap(huVec, dg.vmapM, huM);
			applyIndexMap(huVec, dg.vmapP, huP);

			applyIndexMap(hvVec, dg.vmapM, hvM);
			applyIndexMap(hvVec, dg.vmapP, hvP);

			applyIndexMap(Hvec, dg.vmapM, HM);
			applyIndexMap(Hvec, dg.vmapP, HP);

			// BC's - no flow through walls.
			for (index_type i=0; i < static_cast<index_type>(mapW.size()); ++i) {
				index_type w = mapW[i];
				huP(w) = huM(w) - 2*nxVec(w)*(huM(w)*nxVec(w) + hvM(w)*nyVec(w));
				hvP(w) = hvM(w) - 2*nyVec(w)*(huM(w)*nxVec(w) + hvM(w)*nyVec(w));
			}

			// OBC's - free surface moves up and down according to the tidal forcing.
			for (index_type i=0; i < static_cast<index_type>(mapO.size()); ++i) {
				index_type o = mapO[i];
				huP(o) = huM(o) - 2*nxVec(o)*(huM(o)*nxVec(o) + hvM(o)*nyVec(o));
				hvP(o) = hvM(o) - 2*nyVec(o)*(huM(o)*nxVec(o) + hvM(o)*nyVec(o));
				//huP(o) = huM(o);
				//hvP(o) = hvM(o);
				//hP(o) = HM(o) + amp_tide*std::cos(om_tide*t)*0.5*(std::tanh(0.15/3600*(t-T_tide))+1);
			}

			// well-balancing scheme (star variables).
			real_vector_type hMstar(numFaceNodes), hPstar(numFaceNodes),
				bM(numFaceNodes), bP(numFaceNodes);

			bM = -HM; bP = -HP;

			for (index_type i=0; i < numFaceNodes; ++i) {
				hMstar(i) = std::max(0.0, hM(i) + bM(i) - std::max(bP(i), bM(i)));
				hPstar(i) = std::max(0.0, hP(i) + bP(i) - std::max(bP(i), bM(i)));
			}

			hM = hMstar; hP = hPstar;
			huM = hMstar*(huM/hM); huP = hPstar*(huP/hP);
			hvM = hMstar*(hvM/hM); hvP = hPstar*(hvP/hP);

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
			real_matrix_type& F1 = fds.hu, G1 = fds.hv;
			real_matrix_type F2(Np, K), G2(Np, K), G3(Np, K);

			F2 = (fds.hu*fds.hu)/fds.h + 0.5*g*fds.h*fds.h; G2 = (fds.hu*fds.hv)/fds.h;
			real_matrix_type& F3 = G2;  G3 = (fds.hv*fds.hv)/fds.h + 0.5*g*fds.h*fds.h;

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
			dFlux2 = 0.5*((F2M - F2P)*nxVec + (G2M-G2P)*nyVec - spdMax*dhu - (0.5*g*hM*hM - 0.5*g*hMstar*hMstar)*nxVec);
			dFlux3 = 0.5*((F3M - F3P)*nxVec + (G2M-G2P)*nyVec - spdMax*dhv - (0.5*g*hM*hM - 0.5*g*hMstar*hMstar)*nyVec);

			real_matrix_type dFlux1Mat(Nfp*numFaces, K), dFlux2Mat(Nfp*numFaces, K), dFlux3Mat(Nfp*numFaces, K);

			vectorToFull(dFlux1, dFlux1Mat, byRowsOpt);
			vectorToFull(dFlux2, dFlux2Mat, byRowsOpt);
			vectorToFull(dFlux3, dFlux3Mat, byRowsOpt);

			// Assumes PDE has been left-multiplied by local inverse mass matrix, so all we have left
			// is the differentiation matrix contribution, and the surface integral.
			// RHS == -Flux divergence + surface integral contributions.

			// Flux divergence:
			fds.RHS1 = -(dg.rx*sum(dg.Dr(ii,kk)*F1(kk,jj), kk) + dg.sx*sum(dg.Ds(ii,kk)*F1(kk,jj), kk));
			fds.RHS1+= -(dg.ry*sum(dg.Dr(ii,kk)*G1(kk,jj), kk) + dg.sy*sum(dg.Ds(ii,kk)*G1(kk,jj), kk));

			fds.RHS2 = -(dg.rx*sum(dg.Dr(ii,kk)*F2(kk,jj), kk) + dg.sx*sum(dg.Ds(ii,kk)*F2(kk,jj), kk));
			fds.RHS2+= -(dg.ry*sum(dg.Dr(ii,kk)*G2(kk,jj), kk) + dg.sy*sum(dg.Ds(ii,kk)*G2(kk,jj), kk));

			fds.RHS3 = -(dg.rx*sum(dg.Dr(ii,kk)*F3(kk,jj), kk) + dg.sx*sum(dg.Ds(ii,kk)*F3(kk,jj), kk));
			fds.RHS3+= -(dg.ry*sum(dg.Dr(ii,kk)*G3(kk,jj), kk) + dg.sy*sum(dg.Ds(ii,kk)*G3(kk,jj), kk));

			// Surface integral contributions
			const real_matrix_type& Fscale = dg.Fscale;
			const real_matrix_type& Lift   = dg.Lift;

			real_matrix_type surfaceRHS1(Nfp*numFaces, K), surfaceRHS2(Nfp*numFaces, K), surfaceRHS3(Nfp*numFaces, K);

			// Apply Jacobian scaling.
			surfaceRHS1 = Fscale*dFlux1Mat;
			surfaceRHS2 = Fscale*dFlux2Mat;
			surfaceRHS3 = Fscale*dFlux3Mat;

			// Integrate using Lifting operator, add to RHS
			fds.RHS1+= sum(Lift(ii,kk)*surfaceRHS1(kk,jj), kk);
			fds.RHS2+= sum(Lift(ii,kk)*surfaceRHS2(kk,jj), kk);
			fds.RHS3+= sum(Lift(ii,kk)*surfaceRHS3(kk,jj), kk);

			// Add source terms
			real_matrix_type u(Np, K), v(Np, K); 
			u = fds.hu/fds.h;
			v = fds.hv/fds.h;

			// bottom topography
			real_matrix_type sourcex(Np, K);
			real_matrix_type sourcey(Np, K);

			sourcex = g*fds.h*fds.Hx;
			sourcey = g*fds.h*fds.Hy;

			fds.RHS2+= sourcex;
			fds.RHS3+= sourcey;

			// bottom drag
			real_matrix_type norm_u(Np,K);
			norm_u = blitz::sqrt(u*u + v*v);
			fds.RHS2+= -phys.CD*u*norm_u;
			fds.RHS3+= -phys.CD*v*norm_u;

			// Coriolis force
			fds.RHS2+=  phys.f*fds.hv;
			fds.RHS3+= -phys.f*fds.hu;
		} // computeRHS

		void readDepthData(const std::string& depthFile, real_matrix_type& H) {
			firstIndex ii;
			CSVFileReader reader(depthFile);

			string line;
			index_type count = 0;
			index_type K = H.cols();
			index_type Np = H.rows();
			real_type val;
			real_vector_type depthData(Np*K);
			depthData = 0*ii;

			while (reader.parseRowValues(val)) {
				depthData(count) = val;
				++count;
			}

			count = 0;
			for (index_type k=0; k < K; ++k) {
				for (index_type n=0; n < Np; ++n) {
					real_type val = depthData(count);

					if (val < 300.0)
						val = 300.0;

					H(n,k) = val;
					++count;
				}
			}
		}

		void buildSpongeCoeff(const DGContext2D& dg, real_type spongeStrength, real_type radInfl, real_matrix_type& spongeCoeff) {

			const index_hashmap& bcHash = dg.bcHash;
			const std::vector<index_type>& mapO = bcHash.at(BCTag::Out);
			index_type outLength = static_cast<index_type>(mapO.size());

			index_vector_type vmapO(outLength);
			const index_vector_type& vmapM = dg.vmapM;

			index_type numFaceNodes = dg.NumFaces*dg.Nfp*dg.K;

			real_vector_type xVec(numFaceNodes), yVec(numFaceNodes);

			fullToVector(dg.x, xVec, false);
			fullToVector(dg.y, yVec, false);

			for(index_type k=0; k < dg.K; ++k) {
				for(index_type n=0; n < dg.Np; ++n) {
					real_type x0 = dg.x(n,k), y0 = dg.y(n,k);
					real_type closestDist = 1.0e12;

					for (index_type i=0; i < outLength; ++i) {
						index_type o = vmapM(mapO[i]);
						real_type dist = std::hypot(x0-xVec(o), y0-yVec(o));
						if ( dist < radInfl && dist < closestDist )
							closestDist = dist;
					}

					if (closestDist < 1.0e12) {
						real_type c_sponge = spongeStrength*(1.0-(closestDist/radInfl));
						spongeCoeff(n, k) = c_sponge;
					}
				}
			}
		}
	} // namespace sw2d
} // namespace blitzdg
