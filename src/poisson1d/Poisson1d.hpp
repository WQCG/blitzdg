// Copyright (C) 2017-2019  Waterloo Quantitative Consulting Group, Inc.
// See COPYING and LICENSE files at project root for more details.

/**
 * @file Poisson1d.hpp
 * @brief Short header file for 1D Poisson solver's main.cpp file.
 */
#pragma once
#include "Nodes1DProvisioner.hpp"
#include "BlitzHelpers.hpp"
#include "LinAlgHelpers.hpp"
#include "Types.hpp"
#include <string>

namespace blitzdg {
	namespace poisson1d {
		class PoissonOperator {
			const Nodes1DProvisioner* nodes1D;
		public:
			explicit PoissonOperator(const Nodes1DProvisioner& _nodes1D)
                : nodes1D{ &_nodes1D }
            {}
		
			bool operator()(const real_vector_type & in, real_vector_type & out) {
				// Blitz indices
				blitz::firstIndex ii;
				blitz::secondIndex jj;
				blitz::thirdIndex kk;
				const real_matrix_type& Dr = nodes1D->get_Dr();
				const real_matrix_type& rx = nodes1D->get_rx();
				const real_matrix_type& J = nodes1D->get_J();
				const real_matrix_type& Lift = nodes1D->get_Lift();
				const real_matrix_type& Fscale = nodes1D->get_Fscale();
				const real_matrix_type& nx = nodes1D->get_nx();
				const real_matrix_type& Vinv = nodes1D->get_Vinv();

				// Get volume to surface maps.
				const index_vector_type& vmapM = nodes1D->get_vmapM();
				const index_vector_type& vmapP = nodes1D->get_vmapP();

				// boundary indices.
				index_type mapO = nodes1D->get_mapO();
				index_type mapI = nodes1D->get_mapI();

				index_type numFaces = nodes1D->NumFaces;
				index_type Nfp = nodes1D->NumFacePoints;
				index_type Np = nodes1D->get_NumLocalPoints();
				index_type K = nodes1D->get_NumElements();

				real_matrix_type VinvTrans(Np, Np);
				VinvTrans = Vinv(jj,ii);

				real_vector_type du(numFaces*Nfp*K);
				real_vector_type uM(numFaces*Nfp*K);
				real_vector_type uP(numFaces*Nfp*K);
				real_vector_type nxVec(numFaces*Nfp*K);
				real_matrix_type u(Np, K);

				const bool byRowsOpt = false;
				blitzdg::vectorToFull(in, u, byRowsOpt);
				blitzdg::fullToVector(nx, nxVec, byRowsOpt);

				real_matrix_type ux(Np, K);
				real_vector_type uVec(Np*K);
				real_vector_type uxVec(Np*K);

				blitzdg::fullToVector(u, uVec, byRowsOpt);

				ux = rx*sum(Dr(ii,kk)*u(kk,jj), kk);
				blitzdg::fullToVector(ux, uxVec, byRowsOpt);

				blitzdg::applyIndexMap(uVec, vmapM, uM);
				blitzdg::applyIndexMap(uVec, vmapP, uP);

				// impose boundary condition -- Dirichlet BC's
				uP(mapI) = -uM(mapI);
				uP(mapO) = -uM(mapO);

				// Compute jump in flux:
				du = (uM - uP);

				real_matrix_type duMat(Nfp*numFaces, K);
				blitzdg::vectorToFull(du, duMat, byRowsOpt);

				// Auxiliary variable to hold first derivative.
				real_matrix_type q(Np, K);
				q = rx*sum(Dr(ii,kk)*u(kk,jj), kk);

				real_matrix_type surfaceRHS(Nfp*numFaces, K);
				real_vector_type qVec(Np*K);

				surfaceRHS = Fscale*(nx*duMat/2.0);

				q  = rx*(sum(Dr(ii,kk)*u(kk,jj), kk));
				q -= sum(Lift(ii,kk)*(surfaceRHS(kk,jj)), kk);

				blitzdg::fullToVector(q, qVec, byRowsOpt);
				real_vector_type dq(numFaces*Nfp*K);
				real_vector_type qM(numFaces*Nfp*K);
				real_vector_type qP(numFaces*Nfp*K);
				real_vector_type uxM(numFaces*Nfp*K);
				real_vector_type uxP(numFaces*Nfp*K);

				blitzdg::applyIndexMap(qVec, vmapM, qM);
				blitzdg::applyIndexMap(uxVec, vmapM, uxM);
				blitzdg::applyIndexMap(uxVec, vmapP, uxP);

				// Neumann BC's on q
				uxP(mapI) = uxM(mapI);
				uxP(mapO) = uxM(mapO);

				qP = 0.5*(uxM + uxP);

				dq = qM -qP;

				real_matrix_type dqMat(Nfp*numFaces, K);
				blitzdg::vectorToFull(dq, dqMat, byRowsOpt);

				const double hmin = 2.0/blitzdg::normMax(rx);
				const double tau = (Np*Np)/hmin;

				surfaceRHS = nx*(dqMat + tau*nx*duMat);

				real_matrix_type MassMatrix(Np, Np);
				MassMatrix = sum(Vinv(kk,ii)*Vinv(kk,jj), kk);

				real_matrix_type tmp(Np, K);
				real_matrix_type result(Np, K);

				tmp  = rx*(sum(Dr(ii,kk)*q(kk,jj), kk));
				tmp -= sum(Lift(ii,kk)*(surfaceRHS(kk,jj)), kk);
				result = J*(sum(MassMatrix(ii,kk)*tmp(kk,jj), kk));

				blitzdg::fullToVector(result, out, byRowsOpt);

				return true;
			} // forward operator
		};

		// identity preconditioner
        class Precon {
        public:
            bool operator()(const real_vector_type& in, real_vector_type& out) const {
                out = in;
                return true;
            }
        };
	} // namespace poisson1d
} // namespace blitzdg
