#include "Poisson2DSparseMatrix.hpp"
#include "Nodes1DProvisioner.hpp"
#include "VandermondeBuilders.hpp"

using blitz::Range;

namespace blitzdg{
    Poisson2DSparseMatrix::Poisson2DSparseMatrix(DGContext2D& dg) {

        const real_matrix_type& Dr = dg.Dr();
        const real_matrix_type& Ds = dg.Ds();

        const real_matrix_type& rx = dg.rx();
        const real_matrix_type& ry = dg.ry();
        const real_matrix_type& sx = dg.sx();
        const real_matrix_type& sy = dg.sy();

        const real_matrix_type& Vinv = dg.Vinv();

        blitz::firstIndex ii;
        blitz::secondIndex jj;
        blitz::thirdIndex kk;

        const index_type Np = dg.numLocalPoints();
        const index_type K = dg.numElements();
        const index_type Nfaces = dg.numFaces();
        const index_type Nfp = dg.numFacePoints();

        real_matrix_type MassMatrix(Np, Np);

        MassMatrix = blitz::sum(Vinv(kk,ii)*Vinv(kk,jj),kk);

        const index_matrix_type& Fmask = dg.fmask();

        real_tensor3_type massEdge(Nfp, Nfp, Nfaces);

        Nodes1DProvisioner nodes1d(Nfp-1, 1, -1., 1.);

        real_matrix_type V1D(Nfp, Nfp), V1Dinv(Nfp,Nfp);
        VandermondeBuilders vandermonde;

        const real_vector_type& r = dg.r(), s = dg.s();

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

        for (index_type k=0; k < K; ++k) {
            //rows1 = ((k1-1)*Np+1:k1*Np)'*ones(1,Np); cols1 = rows1';

            index_matrix_type rows1(Np, Np), cols1(Np, Np);
            rows1 = ((ii*Np)+1) * (0*jj + 1);
            cols1 = rows1(jj, ii);

            real_matrix_type Dx(Np, Np), Dy(Np, Np);

            Dx = rx(1, k)*Dr + sx(1, k)*Ds;
            Dy = ry(1, k)*Dr + sy(1, k)*Ds;




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