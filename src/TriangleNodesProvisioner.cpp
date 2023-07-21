// Copyright (C) 2017-2022  Waterloo Quantitative Consulting Group, Inc.
// See COPYING and LICENSE files at project root for more details.

#include "TriangleNodesProvisioner.hpp"
#include "DenseCholeskyFactorizer.hpp"
#include "DGContext2D.hpp"
#include "GaussFaceContext2D.hpp"
#include "CSCMatrix.hpp"
#include "BlitzHelpers.hpp"
#include "DenseMatrixInverter.hpp"
#include "MeshManager.hpp"
#include "Constants.hpp"
#include "Types.hpp"
#include "BCtypes.hpp"
#include <blitz/array.h>
#include <cmath>
#include <limits>
#include <memory>

using blitz::firstIndex;
using blitz::Range;
using blitz::secondIndex;
using blitz::sum;
using blitz::thirdIndex;
using std::numeric_limits;
using std::unique_ptr;
using std::abs;
using std::sqrt;
using std::log;
using std::vector;
using std::pow;
using std::exp;

namespace blitzdg {
    const index_type TriangleNodesProvisioner::NumFaces = 3;
    const real_type TriangleNodesProvisioner::NodeTol = 1.e-5;
    const real_type pi = blitzdg::constants::pi;

    TriangleNodesProvisioner::TriangleNodesProvisioner(index_type _NOrder, const MeshManager& _MeshManager)
        : NumElements{ _MeshManager.get_NumElements() }, NOrder{ _NOrder },
        NumLocalPoints{ (_NOrder + 2)*(_NOrder+1)/2 },
        NumFacePoints{ _NOrder + 1},
        xGrid{ new real_matrix_type((_NOrder + 2)*(_NOrder+1)/2, _MeshManager.get_NumElements()) },
        yGrid{ new real_matrix_type((_NOrder + 2)*(_NOrder+1)/2, _MeshManager.get_NumElements()) },
        rGrid{ new real_vector_type((_NOrder + 2)*(_NOrder+1)/2) },
        sGrid{ new real_vector_type((_NOrder + 2)*(_NOrder+1)/2) },
        V{ new real_matrix_type((_NOrder + 2)*(_NOrder+1)/2, (_NOrder + 2)*(_NOrder+1)/2) }, 
        Dr{ new real_matrix_type((_NOrder + 2)*(_NOrder+1)/2, (_NOrder + 2)*(_NOrder+1)/2) },
        Ds{ new real_matrix_type((_NOrder + 2)*(_NOrder+1)/2, (_NOrder + 2)*(_NOrder+1)/2) },
        Drw{ new real_matrix_type((_NOrder + 2)*(_NOrder+1)/2, (_NOrder + 2)*(_NOrder+1)/2) },
        Dsw{ new real_matrix_type((_NOrder + 2)*(_NOrder+1)/2, (_NOrder + 2)*(_NOrder+1)/2) },
        Lift{ new real_matrix_type((_NOrder + 2)*(_NOrder+1)/2, (_NOrder+1)*NumFaces) },
        J{ new real_matrix_type((_NOrder + 2)*(_NOrder+1)/2, _MeshManager.get_NumElements()) },
        rx{ new real_matrix_type((_NOrder + 2)*(_NOrder+1)/2, _MeshManager.get_NumElements()) },
        sx{ new real_matrix_type((_NOrder + 2)*(_NOrder+1)/2, _MeshManager.get_NumElements()) },
        ry{ new real_matrix_type((_NOrder + 2)*(_NOrder+1)/2, _MeshManager.get_NumElements()) },
        sy{ new real_matrix_type((_NOrder + 2)*(_NOrder+1)/2, _MeshManager.get_NumElements()) },
        nx{ new real_matrix_type((_NOrder+1)*NumFaces, _MeshManager.get_NumElements()) },
        ny{ new real_matrix_type((_NOrder+1)*NumFaces, _MeshManager.get_NumElements()) },
        Vinv{ new real_matrix_type((_NOrder + 2)*(_NOrder+1)/2, (_NOrder + 2)*(_NOrder+1)/2) },
        Filter{ new real_matrix_type((_NOrder + 2)*(_NOrder+1)/2, (_NOrder + 2)*(_NOrder+1)/2) },
        Fmask{ new index_matrix_type( _NOrder+1, NumFaces) },
        Fscale{ new real_matrix_type((_NOrder+1)*NumFaces, _MeshManager.get_NumElements()) },
        Fx{ new real_matrix_type((_NOrder+1)*NumFaces, _MeshManager.get_NumElements()) },
        Fy{ new real_matrix_type((_NOrder+1)*NumFaces, _MeshManager.get_NumElements()) },
        vmapM{ new index_vector_type((_NOrder+1)*NumFaces*_MeshManager.get_NumElements()) },
        vmapP{ new index_vector_type((_NOrder+1)*NumFaces*_MeshManager.get_NumElements()) },
        mapP{ new index_vector_type((_NOrder+1)*NumFaces*_MeshManager.get_NumElements()) },
        BCmap{ new index_hashmap()},
        Mesh2D { _MeshManager },
        Nodes1D{ new Nodes1DProvisioner(_NOrder, 5, -1.0, 1.0) },
		Jacobi{}, Vandermonde{}, LinSolver{}, Inverter{}, CholeskyFactorizer{}
    {
        // Nodal construction required for physical simulations.
        buildNodes();
        buildLift();
        buildPhysicalGrid();
        buildMaps();
    }

    CubatureContext2D TriangleNodesProvisioner::buildCubatureVolumeMesh(index_type NCubature) {
        firstIndex ii;
        secondIndex jj;
        thirdIndex kk;

        TriangleCubatureRules cubature(NCubature);
        index_type Ncub = cubature.NumCubaturePoints();
        real_vector_type rcub(Ncub), scub(Ncub), wcub(Ncub);

        rcub = cubature.rCoord();
        scub = cubature.sCoord();
        wcub = cubature.weights();

        real_matrix_type Vcub(Ncub, NumLocalPoints), Drcub(Ncub, NumLocalPoints), Dscub(Ncub, NumLocalPoints),
            Drwcub(Ncub, NumLocalPoints), Dswcub(Ncub, NumLocalPoints);
        computeInterpMatrix(rcub, scub, Vcub);
        
        real_matrix_type V2Dr(Ncub, NumLocalPoints), V2Ds(Ncub, NumLocalPoints);

        computeGradVandermondeMatrix(NOrder, rcub, scub, V2Dr, V2Ds);
        computeDifferentiationMatrices(V2Dr, V2Ds, *V, Vcub, Drcub, Dscub, Drwcub, Dswcub);

        real_matrix_type sJcub(Ncub, NumElements), Jcub(Ncub, NumElements),
            xrcub(Ncub, NumElements), yrcub(Ncub, NumElements),
            xscub(Ncub, NumElements), yscub(Ncub, NumElements),
            rxcub(Ncub, NumElements), rycub(Ncub, NumElements),
            sxcub(Ncub, NumElements), sycub(Ncub, NumElements);

        const real_matrix_type& x = *xGrid, y = *yGrid;
        xrcub = sum(Drcub(ii,kk)*x(kk,jj), kk);
        xscub = sum(Dscub(ii,kk)*x(kk,jj), kk);

        yrcub = sum(Drcub(ii,kk)*y(kk,jj), kk);
        yscub = sum(Dscub(ii,kk)*y(kk,jj), kk);

        Jcub = -xscub*yrcub + xrcub*yscub;

        rxcub = yscub/Jcub;
        sxcub =-yrcub/Jcub;
        rycub =-xscub/Jcub;
        sycub = xrcub/Jcub;

        // repair Jacobian components on nodal points as well
        real_matrix_type& Jac = *J, D_r = *Dr, D_s = *Ds;

        real_matrix_type& r_x = *rx;
        real_matrix_type& r_y = *ry;
        real_matrix_type& s_x = *sx;
        real_matrix_type& s_y = *sy;

        index_type Np = NumLocalPoints, K = NumElements;
        real_matrix_type xr(Np, K), yr(Np, K), xs(Np, K), ys(Np, K);

        xr = sum(D_r(ii,kk)*x(kk,jj), kk);
        yr = sum(D_r(ii,kk)*y(kk,jj), kk);
        xs = sum(D_s(ii,kk)*x(kk,jj), kk);
        ys = sum(D_s(ii,kk)*y(kk,jj), kk);

        Jac = xr*ys - xs*yr;

        // Invert the 2x2 mapping matrix.
        r_x = ys/Jac;
        r_y =-xs/Jac;
        s_x =-yr/Jac;
        s_y = xr/Jac;
        

        // compute cholesky factors.
        real_matrix_type W(Ncub, NumElements), xcub(Ncub, NumElements), ycub(Ncub, NumElements);
        real_vector_type ones(NumElements);

        ones = 0*ii + 1;

        W = wcub(ii)*ones(jj);
        W *= Jcub;

        xcub = sum(Vcub(ii,kk)*x(kk,jj), kk);
        ycub = sum(Vcub(ii,kk)*y(kk,jj), kk);

        real_tensor3_type MM(NumLocalPoints, NumLocalPoints, NumElements), MMchol(NumLocalPoints, NumLocalPoints, NumElements);
        MM = 0.;
        MMchol = 0.;


        real_matrix_type Jcubdiag(Ncub, Ncub), tmp(Ncub, NumLocalPoints), MMinvtmp(NumLocalPoints, NumLocalPoints), upperFactor(NumLocalPoints, NumLocalPoints),
            MMtmp(NumLocalPoints,  NumLocalPoints);
        Jcubdiag = 0.0;

        // build custom mass matrix and its pre-factorization on each Element.
        for (index_type k=0; k< NumElements; ++k) {
            for (index_type i=0; i < Ncub; ++i) {
                Jcubdiag(i,i) = Jcub(i, k)*wcub(i);
            }

            tmp = sum(Jcubdiag(ii, kk)*Vcub(kk, jj), kk);
            MMtmp = sum(Vcub(kk, ii)*tmp(kk,jj), kk);
            MM(Range::all(), Range::all(), k) = MMtmp; 
            
            upperFactor = 0.0;
            CholeskyFactorizer.computeCholesky(MMtmp, upperFactor);
            MMchol(Range::all(), Range::all(), k) = upperFactor;
        }

        return CubatureContext2D {
            NCubature,
            Ncub,
            std::make_shared<real_vector_type>(rcub),
            std::make_shared<real_vector_type>(scub),
            std::make_shared<real_vector_type>(wcub),
            std::make_shared<real_matrix_type>(Vcub),
            std::make_shared<real_matrix_type>(rxcub),
            std::make_shared<real_matrix_type>(sxcub),
            std::make_shared<real_matrix_type>(rycub),
            std::make_shared<real_matrix_type>(sycub),
            std::make_shared<real_matrix_type>(Jcub),
            std::make_shared<real_matrix_type>(Drcub),
            std::make_shared<real_matrix_type>(Dscub),
            std::make_shared<real_tensor3_type>(MM),
            std::make_shared<real_tensor3_type>(MMchol),
            std::make_shared<real_matrix_type>(xcub),
            std::make_shared<real_matrix_type>(ycub),
            std::make_shared<real_matrix_type>(W)    
        };
    }


    GaussFaceContext2D TriangleNodesProvisioner::buildGaussFaceNodes(index_type NGauss) {
        blitz::firstIndex ii;
        blitz::secondIndex jj;
        blitz::thirdIndex kk;

        index_type Np = NumLocalPoints;

        real_vector_type z(NGauss), w(NGauss);
        real_vector_type face1r(NGauss), face2r(NGauss), face3r(NGauss);
        real_vector_type face1s(NGauss), face2s(NGauss), face3s(NGauss);

        // alpha = beta = 0, for Legendre polynomials.
        Jacobi.computeJacobiQuadWeights(0, 0, NGauss - 1, z, w);

        face1r = z; face2r = -z; face3r = 0.0*z - 1.0;
        face1s = 0.0*z - 1.0; face2s = z; face3s = -z;

        real_matrix_type V1(NGauss, Np), V2(NGauss, Np), V3(NGauss, Np);
        
        computeVandermondeMatrix(NOrder, face1r, face1s, V1);
        computeVandermondeMatrix(NOrder, face2r, face2s, V2);
        computeVandermondeMatrix(NOrder, face3r, face3s, V3);

        real_matrix_type gaussFaceInterp1(NGauss, Np, ColumnMajorOrder()), 
            gaussFaceInterp2(NGauss, Np, ColumnMajorOrder()), gaussFaceInterp3(NGauss, Np, ColumnMajorOrder());

        real_matrix_type & invV = *Vinv, Drref = *Dr, Dsref = *Ds;
        gaussFaceInterp1 = blitz::sum(V1(ii, kk)*invV(kk, jj), kk);
        gaussFaceInterp2 = blitz::sum(V2(ii, kk)*invV(kk, jj), kk);
        gaussFaceInterp3 = blitz::sum(V3(ii, kk)*invV(kk, jj), kk);

    
        real_matrix_type gaussInterpMatrix(NumFaces*NGauss, Np);
        gaussInterpMatrix(Range(0,          NGauss - 1), Range::all()) = gaussFaceInterp1;
        gaussInterpMatrix(Range(1*NGauss, 2*NGauss - 1), Range::all()) = gaussFaceInterp2;
        gaussInterpMatrix(Range(2*NGauss, 3*NGauss - 1), Range::all()) = gaussFaceInterp3;

        index_hashmap gbcMap = {
            {BCTag::Wall, std::vector<index_type>()},
            {BCTag::Dirichlet, std::vector<index_type>()},
            {BCTag::Neuman, std::vector<index_type>()},
            {BCTag::In, std::vector<index_type>()},
            {BCTag::Out, std::vector<index_type>()},
            {BCTag::Cyl, std::vector<index_type>()},
            {BCTag::Far, std::vector<index_type>()},
            {BCTag::Slip, std::vector<index_type>()}
        };

        real_matrix_type nx(NGauss*NumFaces, NumElements), ny(NGauss*NumFaces, NumElements),
            sJ(NGauss*NumFaces, NumElements), Jac(NGauss*NumFaces, NumElements),
            rx(NGauss*NumFaces, NumElements), ry(NGauss*NumFaces, NumElements),
            sx(NGauss*NumFaces, NumElements), sy(NGauss*NumFaces, NumElements);

        index_vector_type mapM(NGauss*NumFaces*NumElements), mapP(NGauss*NumFaces*NumElements);
        mapM = 0, mapP = 0;

        real_matrix_type& x = *xGrid, y = *yGrid;
        std::vector<index_type> mapB;
        const index_vector_type& bcVec = Mesh2D.get_BCType();

        // global face node id's.
        index_type Nfp = NGauss*NumFaces;
        index_matrix_type nodeIds(Nfp, NumElements, ColumnMajorOrder());
        index_vector_type nodeIdsVec(Nfp*NumElements);
        nodeIds = ii + Nfp*jj;
        fullToVector(nodeIds, nodeIdsVec, false);

        mapM = nodeIdsVec;

        for (index_type f=0; f < NumFaces; ++f) {
            real_matrix_type VM(NGauss, Np), dVMdr(NGauss, Np), dVMds(NGauss, Np);
            VM(Range::all(), Range::all()) = gaussInterpMatrix(Range(f*NGauss, (f+1)*NGauss - 1), Range::all());

            dVMdr = blitz::sum(VM(ii, kk)*Drref(kk, jj), kk);
            dVMds = blitz::sum(VM(ii, kk)*Dsref(kk, jj), kk);

            // get global numbers for local guys.
            index_vector_type ids1(NGauss);
            ids1 = f*NGauss + ii;
            for (index_type k=0; k < NumElements; ++k) {
                real_vector_type gxr(NGauss), gyr(NGauss),
                    gxs(NGauss), gys(NGauss), gJac(NGauss),
                    grx(NGauss), gry(NGauss), gsx(NGauss), gsy(NGauss),
                    gnx(NGauss), gny(NGauss), gsJ(NGauss);

                real_vector_type xk(NumLocalPoints), yk(NumLocalPoints);
                xk = x(Range::all(), k);
                yk = y(Range::all(), k);

                gxr = blitz::sum(dVMdr(ii,jj)*xk(jj), jj);
                gyr = blitz::sum(dVMdr(ii,jj)*yk(jj), jj);
                gxs = blitz::sum(dVMds(ii,jj)*xk(jj), jj);
                gys = blitz::sum(dVMds(ii,jj)*yk(jj), jj);

                gJac = gxr*gys - gxs*gyr;

                // Invert the 2x2 mapping matrix.
                grx = gys/gJac;
                gry =-gxs/gJac;
                gsx =-gyr/gJac;
                gsy = gxr/gJac;

                if (f == 0) { 
                    gnx = -gsx; gny = -gsy; 
                } else if (f == 1) { 
                    gnx = grx + gsx; gny = gry + gsy;
                } else if (f == 2) { 
                    gnx = -grx; gny = -gry;
                }
            
                gsJ = blitz::sqrt(gnx*gnx + gny*gny);
                gnx = gnx / gsJ;
                gny = gny / gsJ;
                gsJ = gsJ * gJac;

                // copy to the global arrays.
                for (index_type ig=0; ig < NGauss; ++ig) {
                    nx(ids1(ig), k)  = gnx(ig);
                    ny(ids1(ig), k)  = gny(ig);
                    sJ(ids1(ig), k)  = gsJ(ig);
                    Jac(ids1(ig), k) = gJac(ig);
                    rx(ids1(ig), k)  = grx(ig);
                    ry(ids1(ig), k)  = gry(ig);
                    sx(ids1(ig), k)  = gsx(ig);
                    sy(ids1(ig), k)  = gsy(ig);
                }
            
                const index_vector_type& E2E = Mesh2D.get_EToE(), E2F = Mesh2D.get_EToF();
                index_type k2 = E2E(NumFaces*k + f), f2 = E2F(NumFaces*k + f);

                index_vector_type ids2(NGauss);
                for (index_type ind=NGauss-1; ind >= 0; --ind) {
                    ids2(ind) = NGauss*(f2+1) - ind - 1;
                }
                if (k != k2) {
                    for (index_type ig=0; ig < NGauss; ++ig) {
                        mapP(ids1(ig) + k*Nfp) = mapM(ids2(ig) + k2*Nfp);
                    }
                } else {
                    for (index_type ig=0; ig < NGauss; ++ig) {
                        index_type bcFaceNode = mapM(ids1(ig) + k*Nfp);
                        mapP(ids1(ig) + k*Nfp) = bcFaceNode;
                        mapB.push_back(bcFaceNode);

                        gbcMap[bcVec(NumFaces*k + f)].push_back(bcFaceNode);
                    }
                }
            }
        }

        // get x & y grids on Gauss face mesh.
        real_matrix_type gx(NumFaces*NGauss, NumElements), gy(NumFaces*NGauss, NumElements);

        gx = blitz::sum(gaussInterpMatrix(ii, kk)*x(kk, jj), kk);
        gy = blitz::sum(gaussInterpMatrix(ii, kk)*y(kk, jj), kk);

        // Form vectorized quadrature weights.
        real_matrix_type W(NumFaces*NGauss, NumElements);
        real_matrix_type wmat(NGauss, NumElements);

        W = 0.0;
        wmat = w(ii) * (0*jj + 1.0);
        W(Range(0, NGauss - 1),        Range::all()) = wmat;
        W(Range(NGauss, 2*NGauss - 1), Range::all()) = wmat;
        W(Range(2*NGauss, 3*NGauss-1), Range::all()) = wmat;

        // Jacobian scaling.
        W *= sJ;

        return GaussFaceContext2D { 
            NGauss, nx, ny, sJ, Jac, rx, ry,
            sx, sy, gbcMap, gx, gy, W, gaussInterpMatrix,
            mapM, mapP
        };
    }

	void TriangleNodesProvisioner::evaluateSimplexPolynomial(const real_vector_type & a, const real_vector_type & b, index_type i, index_type j, real_vector_type & p) const {
		real_vector_type h1(a.length(0));
		real_vector_type h2(b.length(0));

		Jacobi.computeJacobiPolynomial(a, 0.0, 0.0, i, h1);
		Jacobi.computeJacobiPolynomial(b, 2.0*i+1.0, 0.0, j, h2);

		real_vector_type c(b.length(0));

		p = sqrt(2.0)*h1*h2*pow(1.-b, i);
	}

    void TriangleNodesProvisioner::rsToab(const real_vector_type & r, const real_vector_type & s, real_vector_type & a, real_vector_type & b) const {
        index_type Np = r.length(0);
        for( index_type i=0; i < Np; i++) {
            if(s(i) != 1.0) 
                a(i) = 2.0*(1.0+r(i))/(1.0-s(i)) - 1.0;
            else
                a(i) = -1.0;
        }
        b = s;
    }

    void TriangleNodesProvisioner::xyTors(const real_vector_type & x, const real_vector_type & y, real_vector_type & r, real_vector_type & s) const {
        const index_type Np = x.length(0);
        real_vector_type L1(Np), L2(Np), L3(Np);
        L1 = (sqrt(3.0)*y+1.0)/3.0;
        L2 = (-3.0*x - sqrt(3.0)*y + 2.0)/6.0;
        L3 = ( 3.0*x - sqrt(3.0)*y + 2.0)/6.0;

        r =-L2 + L3 - L1;
        s =-L2 - L3 + L1;
    }


    void TriangleNodesProvisioner::computeVandermondeMatrix(index_type N, const real_vector_type & r, const real_vector_type & s, real_matrix_type & V) const {
        const index_type Nr = r.length(0);

        real_vector_type a(Nr), b(Nr);
        rsToab(r, s, a, b);

        // build the Vandermonde matrix
        index_type count = 0;
        for (index_type i=0; i < N + 1; ++i) {
            for (index_type j=0; j < N - i + 1; ++j) {
                real_vector_type p(Nr);
                evaluateSimplexPolynomial(a,b,i,j,p);
                V(Range::all(),count) = p;
                ++count;
            }
        }
    }

    pyarray TriangleNodesProvisioner::computeVandermondeMatrix_numpy(index_type N, const pyarray & r, const pyarray & s) const {
        const index_type Nr = r.shape(0);  // Ns assumed same. // 6 
        const index_type Np = NumLocalPoints;  // 15

        char * rRaw = r.get_data();
        char * sRaw = s.get_data();

        // Target grid.
        real_vector_type rout(Nr), sout(Nr);
        rout = 0.0;
        sout = 0.0;

        std::copy(&rRaw[0], &rRaw[r.shape(0)*sizeof(real_type)], reinterpret_cast<char*>(rout.data()));
        std::copy(&sRaw[0], &sRaw[s.shape(0)*sizeof(real_type)], reinterpret_cast<char*>(sout.data()));

        real_matrix_type V2D(Nr, Np);
        V2D = 0.0;

        computeVandermondeMatrix(N, rout, sout, V2D);

        Py_intptr_t shape[2] = { Nr, Np };
        pyarray Vpy = zeros(2, shape, dtype::get_builtin<real_type>());

        char * raw = reinterpret_cast<char*>(V2D.data());
        std::copy(&raw[0], &raw[Nr*Np*sizeof(real_type)], Vpy.get_data());

        return Vpy;
    }

    void TriangleNodesProvisioner::computeGradVandermondeMatrix(index_type N,  const real_vector_type & r, const real_vector_type & s, real_matrix_type & V2Dr, real_matrix_type & V2Ds) const {
        const index_type Nr = r.length(0);
        real_vector_type a(Nr), b(Nr);

        rsToab(r, s, a, b);

        index_type count = 0;
        for (index_type i=0; i <= N; ++i) {
            for (index_type j=0; j <= N-i; ++j) {

                real_vector_type v2drCol(Nr), v2dsCol(Nr);
                evaluateGradSimplex(a, b, i, j, v2drCol, v2dsCol);

                V2Dr(Range::all(), count) = v2drCol;
                V2Ds(Range::all(), count) = v2dsCol;
                ++count;
            }
        }
    }

    void TriangleNodesProvisioner::computeDifferentiationMatrices(const real_matrix_type & V2Dr, const real_matrix_type & V2Ds, const real_matrix_type & V, const real_matrix_type & Vc, real_matrix_type & Dr, real_matrix_type & Ds, real_matrix_type& Drw, real_matrix_type& Dsw) const {
        firstIndex ii;
        secondIndex jj;
        thirdIndex kk;
 
        const index_type Np = V.rows(); // = V.cols()
        const index_type Nc = V2Dr.rows(); // num cubature points

		// Note: this is not a column major ordering trick. We need these transposes.
        real_matrix_type Vtrans(Np, Np);
        real_matrix_type V2Drtrans(Np, Nc);
        real_matrix_type V2Dstrans(Np, Nc);

        real_matrix_type Drtrans(Np, Nc);
        real_matrix_type Dstrans(Np, Nc);

        Drtrans = 0.*jj;
        Dstrans = 0.*jj;

        Vtrans = V(jj, ii);
        V2Drtrans = V2Dr(jj,ii);
        V2Dstrans = V2Ds(jj,ii);

        // Dr = V2Dr * V^{-1}
        // Drtrans = Vtrans^{-1}* V2Drtrans
        LinSolver.solve(Vtrans, V2Drtrans, Drtrans);

        // LAPACK can mangle our input.
        Vtrans = V(jj,ii);

        // Ds = V2Ds / V;
        LinSolver.solve(Vtrans, V2Dstrans, Dstrans);

        // Take transpose.
        Dr = Drtrans(jj, ii); 
        Ds = Dstrans(jj, ii);

        // Now construct the weak derivatives.
        real_matrix_type VVt(Nc, Nc), VVrt(Nc, Np), VVst(Nc, Np); 
        real_matrix_type VVttrans(Nc, Nc), VVrttrans(Np, Nc), VVsttrans(Np, Nc);
        real_matrix_type Drwtrans(Np, Nc), Dswtrans(Np, Nc), Vctrans(Np, Nc);

        Vctrans = Vc(jj, ii);

        VVt = blitz::sum(Vc(ii,kk)*Vctrans(kk,jj), kk);
        VVrt = blitz::sum(Vc(ii,kk)*V2Drtrans(kk,jj), kk);
        VVst = blitz::sum(Vc(ii,kk)*V2Dstrans(kk,jj), kk);

        VVttrans = VVt(jj,ii);   // Vvttrans must be Nc x Np for this to shake out.
        VVrttrans = VVrt(jj,ii);
        VVsttrans = VVst(jj,ii);

        LinSolver.solve( VVttrans, VVrttrans,  Drwtrans );
        LinSolver.solve( VVttrans, VVsttrans,  Dswtrans );

        Drw = Drwtrans(jj,ii);
        Dsw = Dswtrans(jj,ii);
    }

    void TriangleNodesProvisioner::buildFilter(real_type Nc, index_type s) {
        firstIndex ii;
        secondIndex jj;
        thirdIndex kk;

        real_type alpha = -std::log(std::numeric_limits<real_type>().epsilon());

        real_matrix_type& F = *Filter;
        real_matrix_type& Vref = *V;
        real_matrix_type& Vinvref = *Vinv;

        real_matrix_type Fdiag(NumLocalPoints, NumLocalPoints);
        Fdiag = 0.0*jj;

        // build exponential filter
        index_type count = 0;
        for (index_type i=0; i <= NOrder; ++i) {
            for (index_type j=0; j <= NOrder-i; ++j) {
                if ( (i+j) >= Nc) {
                    real_type k = (static_cast<real_type>(i+j) - Nc) / (static_cast<real_type>(NOrder) - Nc);
                    Fdiag(count, count) = std::exp(-alpha*std::pow(k,s));
                } else {
                    Fdiag(count, count) = 1.0;
                }
                ++count;
            }
        }
        
        real_matrix_type tmp(NumLocalPoints, NumLocalPoints);
        tmp = sum(Fdiag(ii,kk)*Vinvref(kk,jj), kk);

        F = sum(Vref(ii,kk)*tmp(kk,jj), kk);
    }

    void TriangleNodesProvisioner::computeEquilateralNodes(real_vector_type & x, real_vector_type & y) const {
        real_vector_type alphaOptimal(15);

        alphaOptimal = 0.0000,0.0000,1.4152,0.1001,0.2751,0.9800,1.0999,
                       1.2832,1.3648,1.4773,1.4959,1.5743, 1.5770,1.6223,1.6258;

        real_type alpha = 2.0/3.0;
        if (NOrder < 16) {
            alpha = alphaOptimal(NOrder-1);
        }

        const index_type Np = (NOrder+1)*(NOrder+2)/2;

        // Create equidistributed nodes on equilateral triangle.
        real_vector_type L1(Np), L2(Np), L3(Np);

        index_type count = 0;

        for (index_type n=1; n <= NOrder+1; ++n) {
            for (index_type m=1; m <= NOrder+2-n; ++m) {
                L1(count) = (n-1.0)/NOrder;
                L3(count) = (m-1.0)/NOrder;
                ++count;
            }
        }
        L2 = 1.0 - L1 - L3;

        x = -L2+L3; 
        y = (-L2-L3+2*L1)/sqrt(3.0);

        // Compute blending function at each node for each edge.
        real_vector_type blend1(Np), blend2(Np), blend3(Np);
        blend1 = 4*L2*L3; 
        blend2 = 4*L1*L3; 
        blend3 = 4*L1*L2;

        // Amount of warp for each node, for each edge.
        real_vector_type warpf1(Np), warpf2(Np), warpf3(Np), temp(Np);

        temp = L3 - L2; computeWarpFactor(temp, warpf1);
        temp = L1 - L3; computeWarpFactor(temp, warpf2);
        temp = L2 - L1; computeWarpFactor(temp, warpf3);

        // Combine blend & warp
        real_vector_type warp1(Np), warp2(Np), warp3(Np);
		const real_type alphaSquared = alpha*alpha;
        warp1 = blend1*warpf1*(1 + alphaSquared*L1*L1);
        warp2 = blend2*warpf2*(1 + alphaSquared*L2*L2);
        warp3 = blend3*warpf3*(1 + alphaSquared*L3*L3);

        // Accumulate deformations associated with each edge.
        x += + 1*warp1 + cos(2*pi/3)*warp2 + cos(4*pi/3)*warp3;
        y += + 0*warp1 + sin(2*pi/3)*warp2 + sin(4*pi/3)*warp3;
    }

    void TriangleNodesProvisioner::computeWarpFactor(const real_vector_type & r, real_vector_type & warpFactor) const {
        firstIndex ii;
        secondIndex jj;

        const index_type Np = NOrder+1;
        const index_type Nr = r.length(0);

        real_vector_type req(Np), rLGL(Np);
        req = -1.0 + 2*ii/(Np-1.0);

        Jacobi.computeGaussLobottoPoints(0.0, 0.0, NOrder, rLGL);

        real_matrix_type Veq(Np,Np), Veqinv(Np,Np);
        Vandermonde.computeVandermondeMatrix(req, Veq, Veqinv);

        // Evaluate Lagrange polynomial at r.
        real_matrix_type P(Np, Nr), L(Np,Nr);
        for (index_type i=0; i< Np; i++) {
            real_vector_type p(Nr);
            Jacobi.computeJacobiPolynomial(r, 0.0, 0.0, i, p);
            P(i, Range::all()) = p;
        }

        real_matrix_type VeqT(Np,Np);
        VeqT = Veq(jj,ii);

        LinSolver.solve(VeqT, P,  L);

        // warp = L^T * (rLGL - req)
        warpFactor = sum(L(jj,ii)*(rLGL(jj) - req(jj)), jj);

        // Scaling.
        real_vector_type zf(Nr), sf(Nr);
        zf = (abs(r) < 1.0-1e-10);
        sf = 1.0 - (zf*r)*(zf*r); 
        warpFactor = warpFactor / sf + warpFactor*(zf-1.0);
    }

    void TriangleNodesProvisioner::evaluateGradSimplex(const real_vector_type & a, const real_vector_type & b, index_type id, index_type jd, real_vector_type & dpdr, real_vector_type & dpds) const {
        const index_type Np = a.length(0);

        real_vector_type fa(Np), gb(Np), dfa(Np), dgb(Np), tmp(Np);

        Jacobi.computeJacobiPolynomial(a, 0.,       0., id, fa);
        Jacobi.computeJacobiPolynomial(b, 2.*id+1., 0., jd, gb);

        Jacobi.computeGradJacobi(a,       0., 0., id, dfa);
        Jacobi.computeGradJacobi(b, 2.*id+1., 0., jd, dgb);

        // r-derivative
        // d/dr = da/dr d/da + db/dr d/db = (2/(1-s)) d/da = (2/(1-b)) d/da
        dpdr = dfa*gb;
        if (id > 1)
            dpdr *= pow(0.5*(1.-b), id-1);

        // s-derivative
        // d/ds = ((1+a)/2)/((1-b)/2) d/da + d/db
        dpds = dfa*(gb*(0.5*(1+a)));
        if ( id > 1 )
            dpds *= pow(0.5*(1.-b), id-1);

		tmp = dgb*pow(0.5*(1.-b), id);
        if( id > 0 ) {
			tmp -= 0.5*id*gb*pow(0.5*(1-b), id-1);
		}
			

        dpds += fa*tmp;

        // Normalize
        dpdr = pow(2., id+0.5)*dpdr; 
        dpds = pow(2., id+0.5)*dpds;
    }

    void TriangleNodesProvisioner::buildNodes() {
        firstIndex ii;
        secondIndex jj;

        real_vector_type x(NumLocalPoints), y(NumLocalPoints);

        real_vector_type& r = *rGrid.get();
        real_vector_type& s = *sGrid.get();

        computeEquilateralNodes(x, y);
        xyTors(x, y, r, s);

        real_vector_type fmask1(NumFacePoints), fmask2(NumFacePoints), fmask3(NumFacePoints), testField(NumLocalPoints);

        testField = s+1;
        index_type count = 0;
        fmask1 = 0*ii;
        for (index_type i=0; i < NumLocalPoints; i++) {
            if (abs(testField(i)) < NodeTol) {
                fmask1(count) = i;
                ++count;
            }
        }

        testField = r+s;
        count = 0;
        fmask2 = 0*ii;
        for (index_type i=0; i < NumLocalPoints; i++) {
            if (abs(testField(i)) < NodeTol) {
                fmask2(count) = i;
                ++count;
            }
        }

        testField = r+1;
        count = 0;
        fmask3 = 0*ii;
        for (index_type i=0; i < NumLocalPoints; i++) {
            if (abs(testField(i)) < NodeTol) {
                fmask3(count) = i;
                ++count;
            }
        }

        // TODO: Move Fmask computation to a helper method.
        index_matrix_type Fm = *Fmask.get();
        Fm = 0*jj;
        Fm(Range::all(), 0) = fmask1;
        Fm(Range::all(), 1) = fmask2;
        Fm(Range::all(), 2) = fmask3;
    }

    void TriangleNodesProvisioner::buildPhysicalGrid() {
        firstIndex ii;
        secondIndex jj;
        thirdIndex kk;

        // Get element data
        const index_vector_type& EToV = Mesh2D.get_Elements();
        const real_vector_type&  Vert = Mesh2D.get_Vertices();
        index_type NumVertices = Mesh2D.get_NumVerts();

        NumElements = Mesh2D.get_NumElements();

        real_vector_type VX(NumVertices), VY(NumVertices), VZ(NumVertices);
        index_vector_type va(NumElements), vb(NumElements), vc(NumElements);
        index_type count=0;

        // Unpack 1D arrays storing EToV and Vertex coordinates
        for (index_type i=0; i < NumElements; ++i) {
            va(i) = EToV(count);
            vb(i) = EToV(count+1);
            vc(i) = EToV(count+2);
            count += 3;
        }

        count=0;
        for (index_type i=0; i < NumVertices; ++i) {
            VX(i) = Vert(count);
            VY(i) = Vert(count+1);
            VZ(i) = Vert(count+2);
            count += 3;
        }

        real_vector_type VXa(NumElements), VXb(NumElements), VXc(NumElements);
        real_vector_type VYa(NumElements), VYb(NumElements), VYc(NumElements);
        for (index_type i=0; i < NumElements; ++i) {
            VXa(i) = VX(va(i));
            VXb(i) = VX(vb(i));
            VXc(i) = VX(vc(i));

            VYa(i) = VY(va(i));
            VYb(i) = VY(vb(i));
            VYc(i) = VY(vc(i));
        }

        const real_vector_type& r = *rGrid;
        const real_vector_type& s = *sGrid;
        real_matrix_type& V2D = *V;

        index_type Np = NumLocalPoints;
        index_type K = NumElements;

        real_matrix_type V2Dr(Np,Np), V2Ds(Np,Np);
        
        real_matrix_type& D_r = *Dr;
        real_matrix_type& D_s = *Ds;

        real_matrix_type& D_rw = *Drw;
        real_matrix_type& D_sw = *Dsw;

        computeVandermondeMatrix(NOrder, r, s, V2D);
        Inverter.computeInverse(*V, *Vinv);
        
        computeGradVandermondeMatrix(NOrder, r, s, V2Dr, V2Ds);
        computeDifferentiationMatrices(V2Dr, V2Ds, V2D, V2D, D_r, D_s, D_rw, D_sw);

        real_matrix_type& x = *xGrid;
        real_matrix_type& y = *yGrid;

        real_matrix_type& Jac = *J;

        real_matrix_type& r_x = *rx;
        real_matrix_type& r_y = *ry;
        real_matrix_type& s_x = *sx;
        real_matrix_type& s_y = *sy;

        real_matrix_type xr(Np,K), yr(Np,K), xs(Np,K), ys(Np,K);

        x = 0.5*(-(r(ii)+s(ii) + 0.*jj)*VXa(jj) + (1+r(ii) + 0.*jj)*VXb(jj) + (1+s(ii) + 0.*jj)*VXc(jj));
        y = 0.5*(-(r(ii)+s(ii) + 0.*jj)*VYa(jj) + (1+r(ii) + 0.*jj)*VYb(jj) + (1+s(ii) + 0.*jj)*VYc(jj));

        xr = sum(D_r(ii,kk)*x(kk,jj), kk);
        yr = sum(D_r(ii,kk)*y(kk,jj), kk);
        xs = sum(D_s(ii,kk)*x(kk,jj), kk);
        ys = sum(D_s(ii,kk)*y(kk,jj), kk);

        Jac = xr*ys - xs*yr;

        // Invert the 2x2 mapping matrix.
        r_x = ys/Jac;
        r_y =-xs/Jac;
        s_x =-yr/Jac;
        s_y = xr/Jac;

        // Normals and Face scalings:

        const index_matrix_type& Fmsk = *Fmask;
        real_matrix_type& Fscal = *Fscale;

        // interpolate geometric factors to face nodes
        index_type numLocalFaceNodes = NumFaces*NumFacePoints;

        real_matrix_type fxr(numLocalFaceNodes, NumElements), fxs(numLocalFaceNodes, NumElements);
        real_matrix_type fyr(numLocalFaceNodes, NumElements), fys(numLocalFaceNodes, NumElements);

        real_matrix_type& fx = *Fx, fy = *Fy;

        for (index_type k=0; k < NumElements; ++k) {
            index_type count=0;
            for (index_type f=0; f < NumFaces; ++f) {
                for (index_type n=0; n < NumFacePoints; ++n) {
                    fxr(count, k) = xr(Fmsk(n, f), k);
                    fxs(count, k) = xs(Fmsk(n, f), k);
                    fyr(count, k) = yr(Fmsk(n, f), k);
                    fys(count, k) = ys(Fmsk(n, f), k);
                    fx(count,  k) = x(Fmsk(n, f), k);
                    fy(count,  k) = y(Fmsk(n, f), k);
                    ++count;
                }
            }
        }

        // build normals.
        real_matrix_type& n_x = *nx, n_y = *ny;

        real_matrix_type norm(numLocalFaceNodes, NumElements, ColumnMajorOrder());

        index_vector_type fid1(NumFacePoints), fid2(NumFacePoints), fid3(NumFacePoints);

        fid1 = ii;
        fid2 = ii + NumFacePoints;
        fid3 = ii + 2*NumFacePoints;

        for(index_type k=0; k < NumElements; ++k) {
            for (index_type i=0; i < NumFacePoints; i++) {
                // face 1
                n_x(fid1(i), k) =  fyr(fid1(i), k);
                n_y(fid1(i), k) = -fxr(fid1(i), k);

                // face 2
                n_x(fid2(i), k) =  fys(fid2(i), k)-fyr(fid2(i), k);
                n_y(fid2(i), k) = -fxs(fid2(i), k)+fxr(fid2(i), k);

                // face 3
                n_x(fid3(i), k) = -fys(fid3(i), k);
                n_y(fid3(i), k) =  fxs(fid3(i), k);
            }
        }

        // normalise
        norm = sqrt(n_x*n_x + n_y*n_y);

        n_x = n_x/norm;
        n_y = n_y/norm;

        for(index_type k=0; k < NumElements; ++k) {
            index_type count=0;
            for (index_type f=0; f < NumFaces; ++f) {
                for (index_type n=0; n < NumFacePoints; n++) {
                    Fscal(count,k) = norm(count,k)/(Jac(Fmsk(n,f),k));
                    ++count;
                }
            }
        }
    }
    
    void TriangleNodesProvisioner::buildMaps() {
        firstIndex ii;
        secondIndex jj;
        thirdIndex kk;

        index_matrix_type nodeIds(NumLocalPoints, NumElements);
        
        const index_matrix_type& Fmsk = *Fmask;
        index_vector_type& vmM = *vmapM;
        index_vector_type& vmP = *vmapP;
        index_vector_type& mP = *mapP;

        const index_vector_type& E2E = Mesh2D.get_EToE();
        const index_vector_type& E2F = Mesh2D.get_EToF();
        const index_vector_type& E2V = Mesh2D.get_Elements(); 

        real_matrix_type xmat(NumLocalPoints, NumElements, ColumnMajorOrder());
        real_matrix_type ymat(NumLocalPoints, NumElements, ColumnMajorOrder());

        xmat = *xGrid; ymat = *yGrid;

        real_vector_type x(NumElements*NumLocalPoints);
        real_vector_type y(NumElements*NumLocalPoints);
        reshapeMatTo1D(xmat, x.data(), false);
        reshapeMatTo1D(ymat, y.data(), false);

        // Assemble global volume node numbering.
        nodeIds = ii + NumLocalPoints*jj;

        index_tensor3_type vmapM3(NumFacePoints, NumFaces, NumElements);
        index_tensor3_type vmapP3(NumFacePoints, NumFaces, NumElements);
        index_tensor3_type mapP3(NumFacePoints, NumFaces, NumElements);

        vmapM3 = 0*kk; vmapP3 = 0*kk; mapP3 = 0*kk;

        for (index_type k=0; k < NumElements; ++k) {
            for (index_type f=0; f < NumFaces; ++f) {
                for (index_type n=0; n < NumFacePoints; ++n) {
                    vmapM3(n,f,k) = nodeIds(Fmsk(n,f), k);
                }
            }
        }

        // Find index of face nodes with respect to volume node ordering.
        for (index_type n=0; n < NumFacePoints; ++n) {
            for (index_type f=0; f < NumFaces; ++f) {
                for (index_type k=0; k < NumElements; ++k) {
                    // find neighbor.
                    index_type k2 = E2E(NumFaces*k + f);
                    index_type f2 = E2F(NumFaces*k + f);
                    const real_vector_type Vert = Mesh2D.get_Vertices();                    

                    index_type v1 = E2V(k*NumFaces + f); 
                    index_type v2 = E2V(k*NumFaces + ((f+1) % NumFaces));

                    // Compute reference length of edge.
                    real_type vx1 = Vert(3*v1), vy1 = Vert(3*v1 + 1);
                    real_type vx2 = Vert(3*v2), vy2 = Vert(3*v2 + 1);

                    real_type refd = hypot(vx1-vx2, vy1-vy2);

                    // find find volume node numbers of left and right nodes 
                    index_type vidM = vmapM3(n,f,k);
                    real_type x1 = x(vidM), y1 = y(vidM);
                    
                    for (index_type nP=0; nP < NumFacePoints; ++nP) {
                        index_type vidP = vmapM3(nP,f2,k2);
                        real_type x2 = x(vidP), y2 = y(vidP);

                        if (distanceLessThanEps(x1,y1,x2,y2, refd*NodeTol)) {
                            vmapP3(n,f,k) = vidP;
                            mapP3(n,f,k) = nP + f2*NumFacePoints+k2*NumFaces*NumFacePoints;
                        }
                    }
                }
            }
        }

        // reshape to 1D.
        index_type count = 0;
        for (index_type k=0; k < NumElements; ++k) {
            for (index_type f=0; f < NumFaces; ++f) {
                for (index_type n=0; n < NumFacePoints; ++n) {
                    vmM(count) = vmapM3(n,f,k);
                    vmP(count) = vmapP3(n,f,k);
                    mP(count) = mapP3(n,f,k);
                    ++count;
                }
            }
        }


        // identify all boundary nodes
        index_vector_type tmpMapB(NumElements*NumFaces*NumFacePoints);
        index_type numBoundaryNodes = 0;
        for (index_type i=0; i< NumElements*NumFaces*NumFacePoints; ++i) {
            if (vmP(i) == vmM(i)) {
                tmpMapB(numBoundaryNodes) = i;
                ++numBoundaryNodes;
            }
        }
        
        mapB = unique_ptr<index_vector_type>(new index_vector_type(numBoundaryNodes));
        vmapB = unique_ptr<index_vector_type>(new index_vector_type(numBoundaryNodes));

        index_vector_type& mB = *mapB;
        index_vector_type& vmB = *vmapB;
        for (index_type i=0; i < numBoundaryNodes; ++i)  {
            mB(i) = tmpMapB(i);
            vmB(i) = vmM(mB(i));
        }

        buildBCHash();

        real_matrix_type allNodes(NumLocalPoints*NumElements, 2);
        real_vector_type xFlat(NumLocalPoints*NumElements), yFlat(NumLocalPoints*NumElements);
        fullToVector(*xGrid, xFlat, false);
        fullToVector(*yGrid, yFlat, false);

        allNodes(Range::all(), 0) = xFlat;
        allNodes(Range::all(), 1) = yFlat;

        auto foo = uniquetol(allNodes, 1.e-9);
        gather = std::make_unique<std::vector<index_type>>(foo.first);
        scatter = std::make_unique<std::vector<index_type>>(foo.second);
    }

    void TriangleNodesProvisioner::buildBCHash() {
        const index_vector_type& bcVec = Mesh2D.get_BCType();

        buildBCHash(bcVec);
    }

    void TriangleNodesProvisioner::buildBCHash(const index_vector_type& bcType) {
        firstIndex ii;
        secondIndex jj;
    
        index_hashmap& bcMap = *BCmap;

        // create boundary face nodes global numbering.
        index_matrix_type boundaryNodesMat(NumFacePoints, NumFaces*NumElements, ColumnMajorOrder());

        index_vector_type ones(NumFacePoints);
        ones = 0*ii + 1;
        boundaryNodesMat = ones(ii)*bcType(jj);

        index_vector_type boundaryNodes(NumFacePoints*NumFaces*NumElements);
        fullToVector(boundaryNodesMat, boundaryNodes, false);

        index_type count=0;
        for (auto itr = boundaryNodes.begin(); itr != boundaryNodes.end(); ++itr) {
            index_type bct = *itr;

            if (bct != 0) {
                auto search = bcMap.find(bct);
                if (search == bcMap.end())
                    bcMap.insert({ bct, vector<index_type>{ count } });
                else
                    search->second.push_back(count);
            }
            ++count;
        }
    }


    void TriangleNodesProvisioner::buildLift() {
        firstIndex ii;
        secondIndex jj;
        thirdIndex kk;

        real_matrix_type E(NumLocalPoints, NumFaces*NumFacePoints), MassInv(NumLocalPoints, NumLocalPoints);

        real_matrix_type & Liftref = *Lift.get();

        real_vector_type& r = *rGrid.get();
        real_vector_type& s = *sGrid.get();

        index_matrix_type Fm = *Fmask.get();


        real_vector_type faceR(NumFacePoints), faceS(NumFacePoints);
        real_matrix_type V1D(NumFacePoints, NumFacePoints), V1Dinv(NumFacePoints,NumFacePoints);
        real_matrix_type massEdgeInv(NumFacePoints, NumFacePoints);

        real_matrix_type massEdge1(NumFacePoints, NumFacePoints), massEdge2(NumFacePoints,NumFacePoints),
                         massEdge3(NumFacePoints, NumFacePoints);

        // Face 1
        for (index_type i=0; i < NumFacePoints; ++i)
            faceR(i) = r(Fm(i, 0));

        Vandermonde.computeVandermondeMatrix(faceR, V1D, V1Dinv);
        massEdgeInv = sum(V1D(ii,kk)*V1D(jj,kk), kk);

        Inverter.computeInverse(massEdgeInv, massEdge1);

        E = 0.0*jj;
        for (index_type i=0; i < NumFacePoints; ++i) {
            for (index_type j=0; j < NumFacePoints; ++j) {
                E(Fm(i,0),j) = massEdge1(i,j);
            }
        }

        // Face 2
        for (index_type i=0; i < NumFacePoints; ++i)
            faceR(i) = r(Fm(i, 1));

        Vandermonde.computeVandermondeMatrix(faceR, V1D, V1Dinv);
        massEdgeInv = sum(V1D(ii,kk)*V1D(jj,kk), kk);

        Inverter.computeInverse(massEdgeInv, massEdge2);
        for (index_type i=0; i < NumFacePoints; ++i) {
            for (index_type j=NumFacePoints; j < 2*NumFacePoints; ++j) {
                E(Fm(i,1),j) = massEdge2(i,j - NumFacePoints);
            }
        }

        // Face 3.
        for (index_type i=0; i < NumFacePoints; ++i)
            faceS(i) = s(Fm(i, 2));

        Vandermonde.computeVandermondeMatrix(faceS, V1D, V1Dinv);
        massEdgeInv = sum(V1D(ii,kk)*V1D(jj,kk), kk);

        Inverter.computeInverse(massEdgeInv, massEdge3);
        for (index_type i=0; i < NumFacePoints; ++i) {
            for (index_type j=2*NumFacePoints; j < 3*NumFacePoints; ++j) {
                E(Fm(i,2),j) = massEdge3(i,j - 2*NumFacePoints);
            }
        }

        // Compute the Vandermonde guy.
        real_matrix_type& V2D = *V;
        V2D = 0.0*jj;
        computeVandermondeMatrix(NOrder, r, s, V2D);
        
        Inverter.computeInverse(V2D, *Vinv);

        MassInv = 0.0*jj;
        MassInv = sum(V2D(ii,kk)*V2D(jj,kk), kk);

        // Multiply by inverse mass matrix;
        Liftref = sum(MassInv(ii,kk)*E(kk,jj),kk);
    }

    void TriangleNodesProvisioner::computeInterpMatrix(const real_vector_type& rout, const real_vector_type& sout, real_matrix_type& IM) const {
        firstIndex ii;
        secondIndex jj;
        thirdIndex kk;

        index_type length = rout.size();

        real_matrix_type Vout(length, NumLocalPoints);
        computeVandermondeMatrix(NOrder, rout, sout, Vout);

        real_matrix_type& invV = *Vinv;
        IM = sum(Vout(ii,kk)*invV(kk,jj),kk);
    }

    void TriangleNodesProvisioner::splitElements(const real_matrix_type& x, const real_matrix_type& y, const real_matrix_type& field, real_matrix_type& xnew, real_matrix_type& ynew, real_matrix_type& fieldnew) const {
		index_type Np = field.rows();
		index_type K = field.cols();

		real_vector_type rout(Np), sout(Np);

		const index_type N = get_NumFacePoints() - 1;

		index_type count = 0;

		index_matrix_type counter(N+1,N+1);
		counter = -1; // -1 == 'No Value'

		for (index_type n=0; n < N+1; ++n) {
			for (index_type m=0; m < N+2-(n+1); ++m) {
				rout(count) = -1. + 2.*static_cast<double>(m)/static_cast<double>(N);
				sout(count) = -1. + 2.*static_cast<double>(n)/static_cast<double>(N);

				counter(n,m) = count;
				++count;
			}
		}

		real_matrix_type IM(Np,Np);
		IM = 0.;

		computeInterpMatrix(rout, sout, IM);

		vector<index_vector_type> localE2V;

		index_type numLocalElements =0;
		for (index_type n=0; n < N+1; ++n) {
			for (index_type m=0; m < N+1-(n+1); ++m) {
				index_type v1 = counter(n,m), v2 = counter(n,m+1),
					v3 = counter(n+1, m), v4 = counter(n+1,m+1);

				index_vector_type tri123(3);
				tri123 = v1,v2,v3;
				
				localE2V.push_back(tri123);
				if (v4 >= 0) {
					index_vector_type tri243(3);
					tri243 = v2,v4,v3;
					localE2V.push_back(tri243);
					++numLocalElements;
				}
			
				++numLocalElements;
			}
		}

		vector<index_vector_type> E2Vnew;

		for (index_type k=0; k<K; ++k) {
			index_type shift = k*Np;

			for (index_type l=0; l<numLocalElements; ++l) {
				index_vector_type row(3);
				row(0) = localE2V[l](0) + shift;
				row(1) = localE2V[l](1) + shift;
				row(2) = localE2V[l](2) + shift;
				E2Vnew.push_back(row);
			}
		}

		index_type totalNewElements = numLocalElements*K;

		blitz::firstIndex ii;
		blitz::secondIndex jj;
		blitz::thirdIndex kk;

		real_matrix_type resultx(Np,K), resulty(Np,K), resultField(Np,K);
		resultx = blitz::sum(IM(ii,kk)*x(kk,jj),kk);
		resulty = blitz::sum(IM(ii,kk)*y(kk,jj),kk);
		resultField = blitz::sum(IM(ii,kk)*field(kk,jj),kk);
		

		real_vector_type xVec(Np*K), yVec(Np*K), fieldVec(Np*K);
		fullToVector(resultx, xVec, false);
		fullToVector(resulty, yVec, false);
		fullToVector(resultField, fieldVec, false);

		// Unpack 1D arrays storing EToV and Vertex coordinates
		index_vector_type va(totalNewElements), vb(totalNewElements), vc(totalNewElements);
		for (index_type i=0; i < totalNewElements; ++i) {
				va(i) = E2Vnew[i](0);
				vb(i) = E2Vnew[i](1);
				vc(i) = E2Vnew[i](2);
		}

		// resize arrays for the new linear elements.
		xnew.resize(3, totalNewElements);
		ynew.resize(3, totalNewElements);
		fieldnew.resize(3, totalNewElements);

		for (index_type i=0; i < totalNewElements; ++i) {
			 index_type vai = va(i), vbi = vb(i), vci = vc(i);

				xnew(0,i) = xVec(vai);
				xnew(1,i) = xVec(vbi);
				xnew(2,i) = xVec(vci);

				ynew(0,i) = yVec(vai);
				ynew(1,i) = yVec(vbi);
				ynew(2,i) = yVec(vci);
				
				fieldnew(0,i) = fieldVec(vai);
				fieldnew(1,i) = fieldVec(vbi);
				fieldnew(2,i) = fieldVec(vci);
		}
	}

    void TriangleNodesProvisioner::setCoordinates_numpy(const pyarray& x, const pyarray& y) {
		char * xRaw = x.get_data();
        char * yRaw = y.get_data();

        std::copy(&xRaw[0], &xRaw[x.shape(0)*x.shape(1)*sizeof(real_type)] , reinterpret_cast<char*>(xGrid->data()));
        std::copy(&yRaw[0], &yRaw[y.shape(0)*y.shape(1)*sizeof(real_type)] , reinterpret_cast<char*>(yGrid->data()));
    }

    const real_matrix_type & TriangleNodesProvisioner::get_Lift() const {
        return *Lift;
    }

    const real_matrix_type & TriangleNodesProvisioner::get_Dr() const {
        return *Dr;
    }

    const real_matrix_type & TriangleNodesProvisioner::get_Ds() const {
        return *Ds;
    }

    const real_matrix_type & TriangleNodesProvisioner::get_Drw() const {
        return *Drw;
    }

    const real_matrix_type & TriangleNodesProvisioner::get_Dsw() const {
        return *Dsw;
    }

    const real_matrix_type & TriangleNodesProvisioner::get_rx() const {
        return *rx;
    }

    const real_matrix_type & TriangleNodesProvisioner::get_ry() const {
        return *ry;
    }

    const real_matrix_type & TriangleNodesProvisioner::get_sx() const {
        return *sx;
    }

    const real_matrix_type & TriangleNodesProvisioner::get_sy() const {
        return *sy;
    }

    const real_matrix_type & TriangleNodesProvisioner::get_xGrid() const {
        return *xGrid;
    }

    const real_matrix_type & TriangleNodesProvisioner::get_yGrid() const {
        return *yGrid;
    }

    index_type TriangleNodesProvisioner::get_NumLocalPoints() const {
        return NumLocalPoints;
    }

    index_type TriangleNodesProvisioner::get_NumFacePoints() const {
        return NumFacePoints;
    }

    const index_vector_type& TriangleNodesProvisioner::get_vmapM() const {
        return *vmapM;
    }

    const index_vector_type& TriangleNodesProvisioner::get_vmapP() const {
        return *vmapP;
    }

    const index_vector_type& TriangleNodesProvisioner::get_vmapB() const {
        return *vmapB;
    }

    const index_vector_type& TriangleNodesProvisioner::get_mapB() const {
        return *mapB;
    }
    const real_matrix_type& TriangleNodesProvisioner::get_Filter() const{
        return *Filter;
    }
    const index_hashmap& TriangleNodesProvisioner::get_bcMap() const {
        return *BCmap;
    }

    const real_matrix_type& TriangleNodesProvisioner::get_Fscale() const {
        return *Fscale;
    }

    const real_matrix_type& TriangleNodesProvisioner::get_nx() const {
        return *nx;
    }

    const real_matrix_type& TriangleNodesProvisioner::get_ny() const {
        return *ny;
    }

    int TriangleNodesProvisioner::get_NumElements() const {
        return NumElements;
    }

    DGContext2D TriangleNodesProvisioner::get_DGContext() const {
        return DGContext2D {
            NOrder,
            NumLocalPoints,
            NumFacePoints,
            NumElements,
            NumFaces,
            Filter.get(),
            rGrid.get(),
            sGrid.get(),
            xGrid.get(),
            yGrid.get(),
            Fscale.get(),
            Fmask.get(),
            gather.get(),
            scatter.get(),
            V.get(),
            Vinv.get(),
            J.get(),
            rx.get(),
            ry.get(),
            sx.get(),
            sy.get(),
            nx.get(),
            ny.get(),
            Dr.get(),
            Ds.get(),
            Lift.get(),
            vmapM.get(),
            vmapP.get(),
            BCmap.get()
        };
    }
}
