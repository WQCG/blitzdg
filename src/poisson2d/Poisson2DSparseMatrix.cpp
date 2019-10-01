#include "Poisson2DSparseMatrix.hpp"
#include "CSCMatrix.hpp"
#include "Nodes1DProvisioner.hpp"
#include "MeshManager.hpp"
#include "VandermondeBuilders.hpp"

using blitz::Range;

namespace blitzdg{
    Poisson2DSparseMatrix::Poisson2DSparseMatrix(DGContext2D& dg, MeshManager& mshManager) {

        const real_matrix_type& Dr = dg.Dr();
        const real_matrix_type& Ds = dg.Ds();

        const real_matrix_type& rx = dg.rx();
        const real_matrix_type& ry = dg.ry();
        const real_matrix_type& sx = dg.sx();
        const real_matrix_type& sy = dg.sy();
        const real_matrix_type& J  = dg.jacobian();
        const real_matrix_type& Fscale = dg.fscale();

        const real_matrix_type& Vinv = dg.Vinv();

        const index_vector_type& EToE = mshManager.get_EToE();
        const index_vector_type& EToF = mshManager.get_EToF();

        blitz::firstIndex ii;
        blitz::secondIndex jj;
        blitz::thirdIndex kk;

        const index_type N = dg.order();
        const index_type Np = dg.numLocalPoints();
        const index_type K = dg.numElements();
        const index_type Nfaces = dg.numFaces();
        const index_type Nfp = dg.numFacePoints();

        const index_vector_type& vmapM = dg.vmapM();
        const index_vector_type& vmapP = dg.vmapP();

        real_matrix_type MassMatrix(Np, Np);

        MassMatrix = blitz::sum(Vinv(kk,ii)*Vinv(kk,jj),kk);

        const index_matrix_type& Fmask = dg.fmask();

        real_tensor3_type massEdge(Nfp, Nfp, Nfaces);

        Nodes1DProvisioner nodes1d(Nfp-1, 1, -1., 1.);

        real_matrix_type V1D(Nfp, Nfp), V1Dinv(Nfp,Nfp);
        VandermondeBuilders vandermonde;

        const real_vector_type& r = dg.r(), s = dg.s();
        const real_vector_type& nx = dg.nx(), ny  = dg.ny();

        real_vector_type edge(Nfp, 1);
        for (index_type i=0; i < Nfp; ++i)
            edge(i) = r(Fmask(i, 1));

        vandermonde.computeVandermondeMatrix(edge, V1D, V1Dinv);
        massEdge(Range::all(), Range::all(), 0) = blitz::sum(V1Dinv(kk,ii)*V1Dinv(kk,jj), kk);

        for (index_type i=0; i < Nfp; ++i)
            edge(i) = r(Fmask(i, 2));

        vandermonde.computeVandermondeMatrix(edge, V1D, V1Dinv);
        massEdge(Range::all(), Range::all(), 1) = blitz::sum(V1Dinv(kk,ii)*V1Dinv(kk,jj), kk);

        for (index_type i=0; i < Nfp; ++i)
            edge(i) = s(Fmask(i, 2));

        vandermonde.computeVandermondeMatrix(edge, V1D, V1Dinv);
        massEdge(Range::all(), Range::all(), 2) = blitz::sum(V1Dinv(kk,ii)*V1Dinv(kk,jj), kk);

        // allocate storage for sparse matrix entries. 
        real_matrix_type MM(K*Np,Np, 3), OP(K*Np*Np*(1+Nfaces), 3);

        index_vector_type entries(Np*Np, 1), entriesMM(Np*Np, 1);

        entries = ii;
        entriesMM = ii;

        CSCMat sparseMat(Np*K, Np*K, K*Np*Np*(1+Nfaces));

        for (index_type k=0; k < K; ++k) {
            //rows1 = ((k1-1)*Np+1:k1*Np)'*ones(1,Np); cols1 = rows1';

            index_matrix_type rows1(Np, Np), cols1(Np, Np);
            rows1 = (((ii+1)*Np*(k+1)) -1) * (0*jj + 1);
            cols1 = rows1(jj, ii);

            real_matrix_type Dx(Np, Np), Dy(Np, Np);

            Dx = rx(1, k)*Dr + sx(1, k)*Ds;
            Dy = ry(1, k)*Dr + sy(1, k)*Ds;

            real_matrix_type OP11(Np, Np), tmp(Np, Np);

            // OP11 = J(1,k1)*(Dx'*MassMatrix*Dx + Dy'*MassMatrix*Dy);
            tmp  = blitz::sum(Dx(kk, ii)*MassMatrix(kk, jj), kk);
            OP11 = blitz::sum(tmp(ii, kk)*Dx(kk, jj), kk);
            tmp  = blitz::sum(Dy(kk, ii)*MassMatrix(kk, jj), kk);
            OP11+= blitz::sum(tmp(ii, kk)*Dy(kk, jj), kk);

            const index_type k1 = k;
            for (index_type f1=0; f1 < Nfaces; ++f1) {
                index_type k2 = EToE(k1, f1);
                index_type f2 = EToF(k1, f1);

                index_matrix_type rows2(Np, Np), cols2(Np, Np);
                rows2 = ((ii + (ii+1)*Np*k) -1) * (0*jj + 1);
                // check this, modified for zero-based indices.
                cols2 = rows2(jj, ii);

                index_vector_type fidM(Nfp), vidM(Nfp), vidP(Nfp),
                    Fm1(Nfp), Fm2(Nfp);

                fidM = (ii+1)*Nfp*Nfaces*(k+1) - 1;
                for (index_type j=0; j < Nfp; ++j) {
                    vidM(j) = vmapM(fidM(j));
                    vidP(j) = vmapP(fidM(j));
                    Fm1(j) = vidM % Np;
                    Fm2(j) = vidP % Np;
                }

                index_type id = f1*Nfp + k1*Nfp*Nfaces;
                real_type lnx = nx(id), lny = ny(id),
                    lsJ = std::sqrt(nx(id)*nx(id)+ny(id)*ny(id));

                real_type hinv = std::max(Fscale(f1*Nfp, k1), Fscale(f2*Nfp, k2))
                real_matrix_type Dx2(Np, Np), Dy2(Np, Np);

                Dx2 = rx(1, k2)*Dr + sx(1, k2)*Ds;
                Dy2 = ry(1, k2)*Dr + sy(1, k2)*Ds;

                real_matrix_type Dn1(Np, Np), Dn2(Np, Np);
                Dn1 = lnx*Dx  + lny*Dy,
                Dn2 = lnx*Dx2 + lny*Dy2;

                real_matrix_type mmE(Nfp, Nfp);
                mmE = lsJ * massEdge(Range::all(), Range::all(), f1);

                real_type gtau = 100*2*(N+1)*(N+1)*hinv; // set penalty scaling

                switch(BCType(k1,f1)) {
                case bcType::Dirichlet:
                    OP11 += ( gtau*mmE - blitz::sum(mmE(ii,kk)*Dn1(kk,jj), kk) - blitz::sum(Dn1(kk,ii)*mmE(kk,jj), kk) ); // ok
                    break;
                case bcType::Neuman:
                case bcType::Wall:
                default:

                }

                // interior face variational terms
                OP11 += + 0.5*( gtau*mmE - blitz::sum(mmE(ii,kk)*Dn1(kk,jj), kk) - blitz::sum(Dn1(kk,ii)*mmE(kk,jj), kk) );


                real_matrix_type OP12(Np, Np);

                for (index_type j=0; j < Np; ++j) {
                    OP12(Range::all(), Fm2(j)) =             - 0.5*( gtau*mmE(Range::all(),Fm1(j)) );  
                    OP12(Fm1(j), Range::all()) += - 0.5*(      mmE(Fm1(j),Fm1(j))*Dn2(Fm2(j), Range::all()));  
                    OP12(Range::all(), Fm2(j)) += - 0.5*(-blitz::sum(Dn1(kk, jj)*mmE(Range::all(), Fm1(j)) (kk), kk) ); /// fix this!!
                }
                OP(entries(:), :) = [rows1(:), cols2(:), OP12(:)]; */

                entries += Np*Np;
  
            }






        }






    } 
}

/*
massEdge = zeros(Np,Np,Nfaces);
Fm = Fmask(:,1); faceR = r(Fm); 
V1D = Vandermonde1D(N, faceR);  massEdge(Fm,Fm,1) = inv(V1D*V1D');
Fm = Fmask(:,2); faceR = r(Fm); 
V1D = Vandermonde1D(N, faceR);  massEdge(Fm,Fm,2) = inv(V1D*V1D');
Fm = Fmask(:,3); faceS = s(Fm); 
V1D = Vandermonde1D(N, faceS);  massEdge(Fm,Fm,3) = inv(V1D*V1D'); */