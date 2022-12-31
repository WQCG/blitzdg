
// Copyright (C) 2017-2022  Waterloo Quantitative Consulting Group, Inc.
// See COPYING and LICENSE files at project root for more details.

#include "Nodes1DProvisioner.hpp"
#include "CSCMatrix.hpp"
#include "CsvOutputter.hpp"
#include "Types.hpp"
#include "BlitzHelpers.hpp"
#include "Warning.hpp"
#include "ins2d.hpp"
#include "LSERK4.hpp"
#include "LinAlgHelpers.hpp"
#include "QuadNodesProvisioner.hpp"
#include "Poisson2DSparseMatrix.hpp"
#include "VtkOutputter.hpp"
#include "DGContext2D.hpp"
#include "LUSolver.hpp"
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
using std::unordered_map;

string trim(string& str)
{
    size_t first = str.find_first_not_of(' ');
    if (first == std::string::npos)
        return "";
    size_t last = str.find_last_not_of(' ');
    return str.substr(first, (last-first+1));
}

std::string str_toupper(std::string s) {
    std::transform(s.begin(), s.end(), s.begin(),
                    [](unsigned char c){ return std::toupper(c); }
                  );
    return s;
}

int main(int argc, char **argv) {
	using namespace blitzdg;

	printDisclaimer();

    if (argc != 2) {
        cout << "\nPlease provide a namelist file as input.\n";
        cout << "\nUsage: " << argv[0] << " <NameList-File>\n\n";
        return -1;
    }
    std::ifstream rfile;
    std::string line, segment;
    unordered_map<string, string> config;
    rfile.open(argv[1]);
    if (rfile.is_open()) {
        while (std::getline(rfile, line)) {
            if (line.length() == 0) {
                continue;
            }
            if (line.rfind("#", 0) == 0) {
                continue;
            }
            std::stringstream line_strm(line);
            std::vector<std::string> seglist;
            while(std::getline(line_strm, segment, '=')) {
                seglist.push_back(trim(segment).c_str());
            }
            if ((seglist.size() % 2 ) == 1) {
                std::cout << "Error parsing namelist file!\n";
                throw std::runtime_error("Error parsing namelist file!");
            }
            config.insert({str_toupper(seglist[0]), seglist[1]});
        }
        rfile.close();
    }

	// Physical parameters
	const real_type g = std::stof(config["GRAVITATIONALACCELERATION"]);
	const real_type finalTime = std::stof(config["FINALTIME"]);
	real_type t = std::stof(config["INITIALTIME"]);

	// Numerical parameters (N = Order of polynomials)
	const index_type N = std::stoi(config["POLYNOMIALORDER"]);
	const real_type CFL = std::stof(config["CFL"]);
    const real_type filtCutOff = std::stof(config["FILTERCUTOFF"]);
    const real_type filtOrder = std::stof(config["FILTERORDER"]);

	// Build dependencies.
	MeshManager meshManager;
	meshManager.readMesh(config["MESHFILE"].c_str());
	const index_type K = meshManager.get_NumElements();

	// Dependency-inject mesh manager to nodes provisioner.
	QuadNodesProvisioner quadNodesProvisioner(N, meshManager);

	// Pre-processing step - build polynomial dealiasing filter.
	quadNodesProvisioner.buildFilter(filtCutOff*N, filtOrder);

#ifndef __MINGW32__
	VtkOutputter outputter(quadNodesProvisioner);
#else
	CsvOutputter outputter;
#endif

	const real_matrix_type& x = quadNodesProvisioner.get_xGrid();
	const real_matrix_type& y = quadNodesProvisioner.get_yGrid();
	index_type Np = quadNodesProvisioner.get_NumLocalPoints();

	const real_matrix_type& Filt = quadNodesProvisioner.get_Filter();

	real_matrix_type rho(Np, K);
	real_matrix_type u(Np, K);
	real_matrix_type v(Np, K);
	real_matrix_type spd(Np, K);
	real_matrix_type RHS1(Np, K), RHS2(Np, K), RHS3(Np, K);
	real_matrix_type resRK1(Np, K), resRK2(Np, K), resRK3(Np, K);

	firstIndex ii;
	secondIndex jj;
	thirdIndex kk;

	// Intialize fields.
	//rho = -1.*(x/1500.0);
    real_type drho = 0.003;
    rho = 1 + drho*(-tanh(5*(y+.005*x+.075))+1)/2;

	u = 0*jj;
	v = 0*jj;

	real_matrix_type Fscale = quadNodesProvisioner.get_Fscale();

	spd = blitz::sqrt(u*u + v*v);
	const index_vector_type& vmapM = quadNodesProvisioner.get_vmapM();

	real_vector_type spdM(quadNodesProvisioner.NumFaces*quadNodesProvisioner.get_NumFacePoints()*K);
	real_vector_type spdVec(Np*K), fsVec(quadNodesProvisioner.NumFaces*quadNodesProvisioner.get_NumFacePoints()*K);
    real_vector_type yVec(Np*K), rhoVec(Np*K);
    rhoVec = 0.0*ii;
    yVec = 0.0*ii;
    spdM = 0.0*ii;

	fullToVector(spd, spdVec, false);
	fullToVector(Fscale, fsVec, false);
    fullToVector(y, yVec, false);
    fullToVector(rho, rhoVec, false);
	applyIndexMap(spdVec, vmapM, spdM);

    real_type H = abs(blitz::max(yVec) - blitz::min(yVec));
    real_type delta_rho = blitz::max(rhoVec)-blitz::min(rhoVec);


    real_type buoyancy_freq = sqrt(g*delta_rho/H);
    real_type c0_bound = abs(buoyancy_freq*H);


	blitzdg::DGContext2D ctx = quadNodesProvisioner.get_DGContext();
	const index_type bordered = 1;
	Poisson2DSparseMatrix borderedPoisson(ctx, meshManager, bordered);

	const blitzdg::CSCMat& OP = borderedPoisson.getOP(), MM = borderedPoisson.getMM();



	blitzdg::LUSolver solver;


	cout << "LUSolver" << endl;
	
	cout << "Done building LUSolver" << endl;
	real_vector_type soln(Np*K);

	// Compute LU factors.
	cout << "Starting factorization" << endl;
	solver.factorize(OP);
	cout << "Factorization successful" << endl;
	cout << "Solve successful" << endl;



	// Time-step calculation
	real_vector_type Fsc(quadNodesProvisioner.NumFaces*quadNodesProvisioner.get_NumFacePoints()*K);
    Fsc = blitz::abs(fsVec)*spdM;
	real_type dt = std::min(1/blitz::max(CFL/((N+1)*(N+1)*0.5*(Fsc+c0_bound))), CFL/buoyancy_freq);

	RHS1 = 0*jj;
	RHS2 = 0*jj;
	RHS3 = 0*jj;

	resRK1= 0*jj;
	resRK2= 0*jj;
	resRK3= 0*jj;

	index_type count = 0;

	while (t < finalTime) {
		if ((count % 10) == 0) {
			cout << "t=" << t << ", rho_max=" << max(rho) << ", dt=" << dt << "\n";
			string fileName = outputter.generateFileName("rho", count);

			std::map<std::string, real_matrix_type> fields;
			fields.insert({"rho", rho});
			outputter.writeFieldsToFiles(fields, count);
		}

		// SSP RK2
		ins2d::computeRHS(rho, u, v, g, quadNodesProvisioner, RHS1, RHS2, RHS3);
		RHS1 = sum(Filt(ii,kk)*RHS1(kk,jj), kk);
		RHS2 = sum(Filt(ii,kk)*RHS2(kk,jj), kk);
		RHS3 = sum(Filt(ii,kk)*RHS3(kk,jj), kk);

		real_matrix_type rho1(Np,K), u1(Np,K), v1(Np,K);

		rho1 = rho + 0.5*dt*RHS1;
		u1 = u + 0.5*dt*RHS2;
		v1 = v + 0.5*dt*RHS3;

		ins2d::computeRHS(rho1, u1, v1, g, quadNodesProvisioner, RHS1, RHS2, RHS3);
		RHS1 = sum(Filt(ii,kk)*RHS1(kk,jj), kk);
		RHS2 = sum(Filt(ii,kk)*RHS2(kk,jj), kk);
		RHS3 = sum(Filt(ii,kk)*RHS3(kk,jj), kk);

		rho  += dt*RHS1;
		u += dt*RHS2;
		v += dt*RHS3;

		real_type rho_max = normMax(rho);
		if ( std::abs(rho_max) > 1e8  || std::isnan(rho_max) )
			throw std::runtime_error("A numerical instability has occurred!");


		spd = blitz::sqrt(u*u + v*v);
		fullToVector(spd, spdVec, false);
		applyIndexMap(spdVec, vmapM, spdM);

        Fsc = blitz::abs(fsVec)*spdM;
        real_type delta_rho = blitz::max(rhoVec)-blitz::min(rhoVec);
        real_type buoyancy_freq = sqrt(g*delta_rho/H);
        real_type c0_bound = abs(buoyancy_freq*H);
		// dt = std::min(CFL/((N+1)*(N+1)*0.5*Fsc_max), CFL/buoyancy_freq);
        //dt = std::min(blitz::min(CFL/((N+1)*(N+1)*0.5*(Fsc+c0_bound))), CFL/buoyancy_freq);
	    dt = CFL/blitz::max(1/((N+1)*(N+1)*0.5*(Fsc+c0_bound)));
        std::cout << dt << "\n";

		t += dt;
		++count;

	}
	real_matrix_type rhofinal(Np, K);
	rhofinal = rho;

	return 0;
} // end main

namespace blitzdg {
	namespace ins2d {
		void computeRHS(real_matrix_type rho, real_matrix_type u, real_matrix_type v, real_type g, QuadNodesProvisioner& quadNodesProvisioner, real_matrix_type& RHS1, real_matrix_type& RHS2, real_matrix_type& RHS3) {
			// Blitz indices
			firstIndex ii;
			secondIndex jj;
			thirdIndex kk;

			// Differentiation matrices and scaling factors.
			const real_matrix_type& Dr = quadNodesProvisioner.get_Dr();
			const real_matrix_type& Ds = quadNodesProvisioner.get_Ds();

			const real_matrix_type& rx = quadNodesProvisioner.get_rx();
			const real_matrix_type& ry = quadNodesProvisioner.get_ry();
			const real_matrix_type& sx = quadNodesProvisioner.get_sx();
			const real_matrix_type& sy = quadNodesProvisioner.get_sy();

			const real_matrix_type &nx = quadNodesProvisioner.get_nx();
			const real_matrix_type &ny = quadNodesProvisioner.get_ny();

			// Get volume to surface maps.
			const index_vector_type& vmapM = quadNodesProvisioner.get_vmapM();
			const index_vector_type& vmapP = quadNodesProvisioner.get_vmapP();

			// boundary indices.
			const index_hashmap bcHash = quadNodesProvisioner.get_bcMap();
			const std::vector<index_type> mapW = bcHash.at(3);

			index_type numFaces = quadNodesProvisioner.NumFaces;
			index_type Nfp = quadNodesProvisioner.get_NumFacePoints();
			index_type Np = quadNodesProvisioner.get_NumLocalPoints();
			index_type K = quadNodesProvisioner.get_NumElements();

			index_type numFaceNodes = numFaces*Nfp*K;
			index_type numTotalNodes = Np*K;

			real_vector_type drho(numFaceNodes), du(numFaceNodes), dv(numFaceNodes);

			real_vector_type rhoM(numFaceNodes), rhoP(numFaceNodes);
			real_vector_type uM(numFaceNodes), uP(numFaceNodes);
			real_vector_type vM(numFaceNodes), vP(numFaceNodes);

			real_vector_type nxVec(numFaceNodes), nyVec(numFaceNodes);

			real_vector_type rhoVec(numTotalNodes), uVec(numTotalNodes), vVec(numTotalNodes);

			drho = 0.*ii;
			du = 0.*ii;
			dv = 0.*ii;

			uM = 0.*ii;
			uP = 0.*ii;
            vM = 0.*ii;
            vP = 0.*ii;
			nxVec = 0.*ii;
			nyVec = 0.*ii;

			// We want to apply maps to column-wise ordering of the nodes.
			const bool byRowsOpt = false;

			fullToVector(nx, nxVec, byRowsOpt);
			fullToVector(ny, nyVec, byRowsOpt);

			fullToVector(rho, rhoVec, byRowsOpt);
			fullToVector(u, uVec, byRowsOpt);
			fullToVector(v, vVec, byRowsOpt);

			applyIndexMap(rhoVec, vmapM, rhoM);
			applyIndexMap(rhoVec, vmapP, rhoP);

			applyIndexMap(uVec, vmapM, uM);
			applyIndexMap(uVec, vmapP, uP);

			applyIndexMap(vVec, vmapM, vM);
			applyIndexMap(vVec, vmapP, vP);
			

			// Compute jump in state vector
			drho = rhoM - rhoP;
			du = uM - uP;
			dv = vM - vP;

			// calculate flux tensor: '-' traces.
			real_vector_type F1M(numFaceNodes), G1M(numFaceNodes);
            F1M = rhoM*uM; G1M = rhoM*vM;
			real_vector_type F2M(numFaceNodes), G2M(numFaceNodes), G3M(numFaceNodes);

			F2M = uM*uM; G2M = uM*vM; G3M = vM*vM;
			real_vector_type& F3M = G2M;  

			// '+' traces:
            real_vector_type F1P(numFaceNodes), G1P(numFaceNodes);
			F1P = rhoP*uP, G1P = rhoP*vP;
			real_vector_type F2P(numFaceNodes), G2P(numFaceNodes), G3P(numFaceNodes);

			F2P = uP*uP; G2P = uP*vP; G3P = vP*vP;
			real_vector_type& F3P = G2P;

			// Compute required local derivatives for volumetric contributions.
            real_matrix_type rhox(Np, K), rhoy(Np, K), ener(Np, K), vort(Np, K), enerx(Np, K), enery(Np, K);
            rhox = (rx*sum(Dr(ii,kk)*rho(kk,jj), kk) + sx*sum(Ds(ii,kk)*rho(kk,jj), kk));
            rhoy = (ry*sum(Dr(ii,kk)*rho(kk,jj), kk) + sy*sum(Ds(ii,kk)*rho(kk,jj), kk));

            ener = 0.5*(u*u + v*v);
            enerx = (rx*sum(Dr(ii,kk)*ener(kk,jj), kk) + sx*sum(Ds(ii,kk)*ener(kk,jj), kk));
            enery = (ry*sum(Dr(ii,kk)*ener(kk,jj), kk) + sy*sum(Ds(ii,kk)*ener(kk,jj), kk));


            // Add volume contributions
            RHS1 = -u*rhox - v*rhoy;
            RHS2 = -enerx  - v*vort;
            RHS3 = -enery  + u*vort;
            RHS3+= -g*rho;

			// compute maximum linearized wave-speed - for numerical flux stabilization.
			real_vector_type spdM(numFaceNodes), spdP(numFaceNodes);
			real_vector_type spdMax(numFaceNodes);

			spdM = sqrt(uM*uM + vM*vM);
			spdP = sqrt(uP*uP + vP*vP);

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
			dFlux1 = 0.5*((F1M - F1P)*nxVec + (G1M-G1P)*nyVec - spdMax*drho);
			dFlux2 = 0.5*((F2M - F2P)*nxVec + (G2M-G2P)*nyVec - spdMax*du);
			dFlux3 = 0.5*((F3M - F3P)*nxVec + (G3M-G3P)*nyVec - spdMax*dv);

			real_matrix_type dFlux1Mat(Nfp*numFaces, K), dFlux2Mat(Nfp*numFaces, K), dFlux3Mat(Nfp*numFaces, K);

			vectorToFull(dFlux1, dFlux1Mat, byRowsOpt);
			vectorToFull(dFlux2, dFlux2Mat, byRowsOpt);
			vectorToFull(dFlux3, dFlux3Mat, byRowsOpt);

			// Assumes PDE has been left-multiplied by local inverse mass matrix, so all we have left
			// is the differentiation matrix contribution, and the surface integral.
			// RHS == -Flux divergence + surface integral contributions.

			// Surface integral contributions
			const real_matrix_type& Fscale = quadNodesProvisioner.get_Fscale();
			const real_matrix_type& Lift   = quadNodesProvisioner.get_Lift();

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
	} // namespace ins2d
} // namespace blitzdg
