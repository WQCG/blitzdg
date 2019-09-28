// Copyright (C) 2017-2019  Waterloo Quantitative Consulting Group, Inc.
// See COPYING and LICENSE files at project root for more details.

/**
 * @file Poisson2d.hpp
 * @brief Header file for 2D Poisson solver's main.cpp file.
 */

#pragma once
#include "DGContext2D.hpp"
#include "BlitzHelpers.hpp"
#include "LinAlgHelpers.hpp"
#include "BCtypes.hpp"

namespace blitzdg {
	namespace poisson2d {
		class PoissonOperator {
			const DGContext2D* dgContext;
		public:
			explicit PoissonOperator(const DGContext2D& context)
                : dgContext { &context }
            {}

			bool operator()(const real_vector_type & in, real_vector_type & out) {
				const DGContext2D& dg = *dgContext;

				blitz::firstIndex ii;
				blitz::secondIndex jj;
				blitz::thirdIndex kk;

				const index_type numFaces = dg.numFaces();
				const index_type K = dg.numElements();
				const index_type Nfp = dg.numFacePoints();
				const index_type Np = dg.numLocalPoints();

				real_vector_type du(numFaces*Nfp*K);
				real_vector_type uM(numFaces*Nfp*K);
				real_vector_type uP(numFaces*Nfp*K);
				real_vector_type nxVec(numFaces*Nfp*K);
				real_vector_type nyVec(numFaces*Nfp*K);
				real_matrix_type u(Np, K);

				const bool byRowsOpt = false;
				blitzdg::vectorToFull(in, u, byRowsOpt);
				blitzdg::fullToVector(dg.nx(), nxVec, byRowsOpt);
				blitzdg::fullToVector(dg.ny(), nyVec, byRowsOpt);

				real_matrix_type ux(Np,K), uy(Np,K);
				real_vector_type uVec(Np*K), uxVec(Np*K), uyVec(Np*K);

				blitzdg::fullToVector(u, uVec, byRowsOpt);

				const real_matrix_type& Dr = dg.Dr(), Ds = dg.Ds();
				real_matrix_type dudr(Np,K), duds(Np,K);

				dudr = blitz::sum(Dr(ii,kk)*u(kk,jj), kk);
				duds = blitz::sum(Ds(ii,kk)*u(kk,jj), kk);

				ux = dg.rx()*dudr + dg.sx()*duds;
				uy = dg.ry()*dudr + dg.sy()*duds;

				blitzdg::fullToVector(ux, uxVec, byRowsOpt);
				blitzdg::fullToVector(uy, uyVec, byRowsOpt);

				blitzdg::applyIndexMap(uVec, dg.vmapM(), uM);
				blitzdg::applyIndexMap(uVec, dg.vmapP(), uP);

				// impose boundary condition -- Dirichlet BC's (reflective)
				//const index_hashmap& bcMap = dg.bcmap();
				std::vector<index_type> mapD = dg.bcmap().at(BCTag::Wall); // forcing dirichlet for now.
				for (index_type i=0; i < static_cast<int>(mapD.size()); ++i) {
					index_type d = mapD[i];
					uP(d) = -uM(d);
				}

				// impose boundary condition -- Neuman BC's
				// nothing to do => uP(mapN) = uM(mapN) by default.

				du = uM - uP;
				real_matrix_type duMat(Nfp*numFaces, K);
				blitzdg::vectorToFull(du, duMat, byRowsOpt);

				// Auxiliary vector variable to hold first derivative (gradient).
				const real_matrix_type& Lift = dg.lift();
				real_matrix_type qx(Np, K), qy(Np, K);
				qx = ux;
				qy = uy;

				real_matrix_type surfaceRHSx(Nfp*numFaces, K), surfaceRHSy(Nfp*numFaces, K);
				surfaceRHSx = dg.fscale()*(dg.nx()*duMat/2.0);
				surfaceRHSy = dg.fscale()*(dg.ny()*duMat/2.0);

				qx -= blitz::sum(Lift(ii,kk)*surfaceRHSx(kk,jj), kk);
				qy -= blitz::sum(Lift(ii,kk)*surfaceRHSy(kk,jj), kk);

				real_vector_type qxVec(Np*K), qyVec(Np*K);

				blitzdg::fullToVector(qx, qxVec, byRowsOpt);
				blitzdg::fullToVector(qy, qyVec, byRowsOpt);

				real_vector_type dqx(numFaces*Nfp*K), qxM(numFaces*Nfp*K), qxP(numFaces*Nfp*K);
				real_vector_type dqy(numFaces*Nfp*K), qyM(numFaces*Nfp*K), qyP(numFaces*Nfp*K);
				real_vector_type uxM(numFaces*Nfp*K), uxP(numFaces*Nfp*K);
				real_vector_type uyM(numFaces*Nfp*K), uyP(numFaces*Nfp*K);

				blitzdg::applyIndexMap(qxVec, dg.vmapM(), qxM);
				blitzdg::applyIndexMap(qxVec, dg.vmapP(), qxP);

				blitzdg::applyIndexMap(qyVec, dg.vmapM(), qyM);
				blitzdg::applyIndexMap(qyVec, dg.vmapP(), qyP);
				
				blitzdg::applyIndexMap(uxVec, dg.vmapM(), uxM);
				blitzdg::applyIndexMap(uxVec, dg.vmapP(), uxP);

				blitzdg::applyIndexMap(uyVec, dg.vmapM(), uyM);
				blitzdg::applyIndexMap(uyVec, dg.vmapP(), uyP);

				// Neuman BC's imply auxiliary vector "(qx,qy) = 0" on the boundary
				//std::vector<index_type> mapN = dg.bcmap().at(BCTag::Neuman);
				std::vector<index_type> mapN;
				if (dg.bcmap().count( BCTag::Neuman )) {
					mapN = dg.bcmap().at(BCTag::Neuman);
				}

				for (index_type i=0; i < static_cast<int>(mapN.size()); ++i) {
					index_type n = mapN[i];
					uxP(n) = uxM(n) - 2*nxVec(n)*(uxM(n)*nxVec(n) + uyM(n)*nyVec(n));
					uyP(n) = uyM(n) - 2*nyVec(n)*(uxM(n)*nxVec(n) + uyM(n)*nyVec(n));
				}

				// the weirdness that is interior penalty....
				qxP = 0.5*(uxM + uxP);
				qyP = 0.5*(uyM + uyP);

				dqx = qxM - qxP;
				dqy = qyM - qyP;

				real_matrix_type dqxMat(Nfp*numFaces, K), dqyMat(Nfp*numFaces, K);
				blitzdg::vectorToFull(dqx, dqxMat, byRowsOpt);
				blitzdg::vectorToFull(dqy, dqyMat, byRowsOpt);

				// compute penalty scaling.
				const double hmin = 2.0/ blitzdg::normMax(dg.jacobian());
				const double tau = Np/hmin;

				real_matrix_type surfaceRHS(Nfp*numFaces, K);

				surfaceRHS = (dg.nx()*dqxMat + dg.ny()*dqyMat + tau*duMat)/2.0; // x component

				real_matrix_type MassMatrix(Np,Np), tmp(Np,K), result(Np,K);

				const real_matrix_type& Vinv = dg.Vinv();
				MassMatrix = blitz::sum(Vinv(kk,ii)*Vinv(kk,jj), kk); // local one.

				// Commpute Laplacian.
				tmp  = dg.rx()*blitz::sum(Dr(ii,kk)*qx(kk,jj), kk) + dg.sx()*blitz::sum(Ds(ii,kk)*qx(kk,jj), kk);
				tmp += dg.ry()*blitz::sum(Dr(ii,kk)*qy(kk,jj), kk) + dg.sy()*blitz::sum(Ds(ii,kk)*qy(kk,jj), kk);
				tmp -= sum(Lift(ii,kk)*surfaceRHS(kk,jj), kk);

				// Multiply by mass matrix for a symmetric problem.
				//result = dg.jacobian()*(sum(MassMatrix(ii,kk)*tmp(kk,jj),kk));
				result = tmp;

				blitzdg::fullToVector(result, out, byRowsOpt);

				return true;
			}
		};

		// identity preconditioner
        class Precon {
        public:
            bool operator()(const real_vector_type& in, real_vector_type& out) const {
                out = in;
                return true;
            }
        };
	}
}