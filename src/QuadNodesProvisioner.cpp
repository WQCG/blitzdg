// Copyright (C) 2017-2019  Waterloo Quantitative Consulting Group, Inc.
// See COPYING and LICENSE files at project root for more details.

#include "QuadNodesProvisioner.hpp"
#include "DGContext2D.hpp"
#include "CSCMatrix.hpp"
#include "BlitzHelpers.hpp"
#include "DenseMatrixInverter.hpp"
#include "MeshManager.hpp"
#include "Constants.hpp"
#include "Types.hpp"
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
    const index_type QuadNodesProvisioner::NumFaces = 4;
    const real_type QuadNodesProvisioner::NodeTol = 1.e-5;
    const real_type pi = blitzdg::constants::pi;

    QuadNodesProvisioner::QuadNodesProvisioner(index_type _NOrder, const MeshManager& _MeshManager)
        : NumElements{ _MeshManager.get_NumElements() }, NOrder{ _NOrder },
        NumLocalPoints{ (_NOrder + 1)*(_NOrder+1) },
        NumFacePoints{ _NOrder + 1},
        xGrid{ new real_matrix_type((_NOrder + 1)*(_NOrder+1), _MeshManager.get_NumElements()) },
        yGrid{ new real_matrix_type((_NOrder + 1)*(_NOrder+1), _MeshManager.get_NumElements()) },
        rGrid{ new real_vector_type((_NOrder + 1)*(_NOrder+1)) },
        sGrid{ new real_vector_type((_NOrder + 1)*(_NOrder+1)) },
        V{ new real_matrix_type((_NOrder + 1)*(_NOrder+1), (_NOrder + 1)*(_NOrder+1)) }, 
        Dr{ new real_matrix_type((_NOrder + 1)*(_NOrder+1), (_NOrder + 1)*(_NOrder+1)) },
        Ds{ new real_matrix_type((_NOrder + 1)*(_NOrder+1), (_NOrder + 1)*(_NOrder+1)) },
        Drw{ new real_matrix_type((_NOrder + 1)*(_NOrder+1), (_NOrder + 1)*(_NOrder+1)) },
        Dsw{ new real_matrix_type((_NOrder + 1)*(_NOrder+1), (_NOrder + 1)*(_NOrder+1)) },
        Lift{ new real_matrix_type((_NOrder + 1)*(_NOrder+1), (_NOrder+1)*NumFaces) },
        J{ new real_matrix_type((_NOrder + 1)*(_NOrder+1), _MeshManager.get_NumElements()) },
        rx{ new real_matrix_type((_NOrder + 1)*(_NOrder+1), _MeshManager.get_NumElements()) },
        sx{ new real_matrix_type((_NOrder + 1)*(_NOrder+1), _MeshManager.get_NumElements()) },
        ry{ new real_matrix_type((_NOrder + 1)*(_NOrder+1), _MeshManager.get_NumElements()) },
        sy{ new real_matrix_type((_NOrder + 1)*(_NOrder+1), _MeshManager.get_NumElements()) },
        nx{ new real_matrix_type((_NOrder+1)*NumFaces, _MeshManager.get_NumElements()) },
        ny{ new real_matrix_type((_NOrder+1)*NumFaces, _MeshManager.get_NumElements()) },
        Vinv{ new real_matrix_type((_NOrder + 1)*(_NOrder+1), (_NOrder + 1)*(_NOrder+1)) },
        Filter{ new real_matrix_type((_NOrder + 1)*(_NOrder+1), (_NOrder + 1)*(_NOrder+1)) },
        Fmask{ new index_matrix_type( _NOrder+1, NumFaces) },
        Fscale{ new real_matrix_type((_NOrder+1)*NumFaces, _MeshManager.get_NumElements()) },
        vmapM{ new index_vector_type((_NOrder+1)*NumFaces*_MeshManager.get_NumElements()) },
        vmapP{ new index_vector_type((_NOrder+1)*NumFaces*_MeshManager.get_NumElements()) },
        mapP{ new index_vector_type((_NOrder+1)*NumFaces*_MeshManager.get_NumElements()) },
        BCmap{ new index_hashmap()},
        Mesh2D { _MeshManager },
        Nodes1D{ new Nodes1DProvisioner(_NOrder, 5, -1.0, 1.0) },
		Jacobi{}, Vandermonde{}, LinSolver{}, Inverter{}
    {
        // Nodal construction required for physical simulations.
        buildNodes();
        buildLift();
        //buildPhysicalGrid();
        //buildMaps();
    }

    void QuadNodesProvisioner::computeVandermondeMatrix(index_type N, const real_vector_type & r, const real_vector_type & s, real_matrix_type & V) const {
        const index_type Nr = r.length(0), Ns = s.length(0);

        // build the Vandermonde matrix
        index_type count = 0;
        for (index_type i=0; i < N + 1; ++i) {
            for (index_type j=0; j < N + 1; ++j) {
                real_vector_type p(Nr), q(Ns);
                Jacobi.computeJacobiPolynomial(s, 0, 0, i, p);
                Jacobi.computeJacobiPolynomial(r, 0, 0, j, q);
                V(Range::all(),count) = p*q;
                ++count;
            }
        }
        real_matrix_type& VinvRef = *Vinv;
        Inverter.computeInverse(V, VinvRef);
    }

    void QuadNodesProvisioner::computeGradVandermondeMatrix(index_type N,  const real_vector_type & r, const real_vector_type & s, real_matrix_type & V2Dr, real_matrix_type & V2Ds) const {

        index_type ind = 0;
        for (index_type i=0; i <= N; ++i) {
            for (index_type j=0; j <= N; ++j) {
                real_vector_type v1dr(NumLocalPoints),
                    v1ds(NumLocalPoints), v1r(NumLocalPoints), v1s(NumLocalPoints);
                Jacobi.computeGradJacobi(r, 0.0, 0.0, j, v1dr);
                Jacobi.computeGradJacobi(s, 0.0, 0.0, i, v1ds);

                Jacobi.computeJacobiPolynomial(s, 0.0, 0.0, i, v1s);
                Jacobi.computeJacobiPolynomial(r, 0.0, 0.0, j, v1r);

                V2Dr(Range::all(), ind) = v1dr*v1s;
                V2Ds(Range::all(), ind) = v1r*v1ds;
                ++ind;
            }
        }
    }

    void QuadNodesProvisioner::computeDifferentiationMatrices(const real_matrix_type & V2Dr, const real_matrix_type & V2Ds, const real_matrix_type & V, real_matrix_type & Dr, real_matrix_type & Ds, real_matrix_type& Drw, real_matrix_type& Dsw) const {
        firstIndex ii;
        secondIndex jj;
        thirdIndex kk;
 
        const index_type numRowsV = V.rows();
        const index_type numColsV = V.cols();


		// Note: this is not a column major ordering trick. We need these transposes.
        real_matrix_type Vtrans(numColsV, numRowsV);
        real_matrix_type V2Drtrans(numColsV, numRowsV);
        real_matrix_type V2Dstrans(numColsV, numRowsV);

        real_matrix_type Drtrans(numColsV, numRowsV);
        real_matrix_type Dstrans(numColsV, numRowsV);

        Drtrans = 0.*jj;
        Dstrans = 0.*jj;

        Vtrans = V(jj, ii);
        V2Drtrans = V2Dr(jj,ii);
        V2Dstrans = V2Ds(jj,ii);

        // Dr = V2Dr * V^{-1}
        LinSolver.solve(Vtrans, V2Drtrans, Drtrans);

        // LAPACK can mangle our input.
        Vtrans = V(jj,ii);

        // Ds = V2Ds / V;
        LinSolver.solve(Vtrans, V2Dstrans, Dstrans);

        // Take transpose.
        Dr = Drtrans(jj, ii); 
        Ds = Dstrans(jj, ii);

        // Now construct the weak derivatives.
        real_matrix_type VVt(numRowsV,numRowsV), VVrt(numRowsV, numRowsV), VVst(numRowsV, numRowsV); 
        real_matrix_type VVttrans(numRowsV,numRowsV), VVrttrans(numRowsV, numRowsV), VVsttrans(numRowsV, numRowsV); 
        real_matrix_type Drwtrans(numRowsV,numRowsV), Dswtrans(numRowsV,numRowsV);
        VVt = blitz::sum(V(ii,kk)*Vtrans(kk,jj), kk);
        VVrt = blitz::sum(V(ii,kk)*V2Drtrans(kk,jj), kk);
        VVst = blitz::sum(V(ii,kk)*V2Dstrans(kk,jj), kk);

        VVttrans = VVt(jj,ii);
        VVrttrans = VVrt(jj,ii);
        VVsttrans = VVst(jj,ii);

        LinSolver.solve( VVttrans, VVrttrans,  Drwtrans );
        LinSolver.solve( VVttrans, VVsttrans,  Dswtrans );

        Drw = Drwtrans(jj,ii);
        Dsw = Dswtrans(jj,ii);
    }

    void QuadNodesProvisioner::buildFilter(real_type Nc, index_type s) {
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

    void QuadNodesProvisioner::buildNodes() {
        firstIndex ii;
        secondIndex jj;

        real_vector_type x(NumLocalPoints), y(NumLocalPoints);

        real_vector_type& r = *rGrid.get();
        real_vector_type& s = *sGrid.get();

        real_vector_type r1d(NOrder+1), s1d(NOrder+1);

        Jacobi.computeGaussLobottoPoints(0.0, 0.0, NOrder, r1d);
        Jacobi.computeGaussLobottoPoints(0.0, 0.0, NOrder, s1d);

        for (index_type j=0; j < NOrder+1; ++j) {
            for (index_type i=0; i < NOrder+1; ++i) {
                r((NOrder+1)*j+i) = r1d(j);
                s((NOrder+1)*j+i) = s1d(i);
            }
        }

        real_vector_type fmask1(NumFacePoints), fmask2(NumFacePoints), fmask3(NumFacePoints),
            fmask4(NumFacePoints), testField(NumLocalPoints);

        testField = s+1;
        index_type count = 0;
        fmask1 = 0*ii;
        for (index_type i=0; i < NumLocalPoints; i++) {
            if (abs(testField(i)) < NodeTol) {
                fmask1(count) = i;
                ++count;
            }
        }

        testField = r-1;
        count = 0;
        fmask2 = 0*ii;
        for (index_type i=0; i < NumLocalPoints; i++) {
            if (abs(testField(i)) < NodeTol) {
                fmask2(count) = i;
                ++count;
            }
        }

        testField = s-1;
        count = 0;
        fmask3 = 0*ii;
        for (index_type i=0; i < NumLocalPoints; i++) {
            if (abs(testField(i)) < NodeTol) {
                fmask3(count) = i;
                ++count;
            }
        }

        testField = r+1;
        count = 0;
        fmask4 = 0*ii;
        for (index_type i=0; i < NumLocalPoints; i++) {
            if (abs(testField(i)) < NodeTol) {
                fmask4(count) = i;
                ++count;
            }
        }

        // TODO: Move Fmask computation to a helper method.
        index_matrix_type Fm = *Fmask.get();
        Fm = 0*jj;
        Fm(Range::all(), 0) = fmask1;
        Fm(Range::all(), 1) = fmask2;
        Fm(Range::all(), 2) = fmask3;
        Fm(Range::all(), 3) = fmask4;
    }

    void QuadNodesProvisioner::buildPhysicalGrid() {
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
        computeGradVandermondeMatrix(NOrder, r, s, V2Dr, V2Ds);
        computeDifferentiationMatrices(V2Dr, V2Ds, V2D, D_r, D_s, D_rw, D_sw);

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

        for (index_type k=0; k < NumElements; ++k) {
            index_type count=0;
            for (index_type f=0; f < NumFaces; ++f) {
                for (index_type n=0; n < NumFacePoints; ++n) {
                    fxr(count, k) = xr(Fmsk(n, f), k);
                    fxs(count, k) = xs(Fmsk(n, f), k);
                    fyr(count, k) = yr(Fmsk(n, f), k);
                    fys(count, k) = ys(Fmsk(n, f), k);
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
    
    void QuadNodesProvisioner::buildMaps() {
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
    }

    void QuadNodesProvisioner::buildBCHash() {
        const index_vector_type& bcVec = Mesh2D.get_BCType();

        buildBCHash(bcVec);
    }

    void QuadNodesProvisioner::buildBCHash(const index_vector_type& bcType) {
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


    void QuadNodesProvisioner::buildLift() {
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

        real_matrix_type massEdge1(NumFacePoints, NumFacePoints), massEdge2(NumFacePoints, NumFacePoints),
                         massEdge3(NumFacePoints, NumFacePoints), massEdge4(NumFacePoints, NumFacePoints);

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
            faceS(i) = s(Fm(i, 1));

        Vandermonde.computeVandermondeMatrix(faceS, V1D, V1Dinv);
        massEdgeInv = sum(V1D(ii,kk)*V1D(jj,kk), kk);

        Inverter.computeInverse(massEdgeInv, massEdge2);
        for (index_type i=0; i < NumFacePoints; ++i) {
            for (index_type j=NumFacePoints; j < 2*NumFacePoints; ++j) {
                E(Fm(i,1),j) = massEdge2(i,j - NumFacePoints);
            }
        }

        // Face 3
        for (index_type i=0; i < NumFacePoints; ++i)
            faceR(i) = r(Fm(i, 2));

        Vandermonde.computeVandermondeMatrix(faceR, V1D, V1Dinv);
        massEdgeInv = sum(V1D(ii,kk)*V1D(jj,kk), kk);

        Inverter.computeInverse(massEdgeInv, massEdge3);
        for (index_type i=0; i < NumFacePoints; ++i) {
            for (index_type j=2*NumFacePoints; j < 3*NumFacePoints; ++j) {
                E(Fm(i,2),j) = massEdge3(i,j - 2*NumFacePoints);
            }
        }

        // Face 4.
        for (index_type i=0; i < NumFacePoints; ++i)
            faceS(i) = s(Fm(i, 3));

        Vandermonde.computeVandermondeMatrix(faceS, V1D, V1Dinv);
        massEdgeInv = sum(V1D(ii,kk)*V1D(jj,kk), kk);

        Inverter.computeInverse(massEdgeInv, massEdge3);
        for (index_type i=0; i < NumFacePoints; ++i) {
            for (index_type j=3*NumFacePoints; j < 4*NumFacePoints; ++j) {
                E(Fm(i,3),j) = massEdge3(i,j - 3*NumFacePoints);
            }
        }

        // Compute the Vandermonde guy.
        real_matrix_type& V2D = *V;
        V2D = 0.0*jj;
        computeVandermondeMatrix(NOrder, r, s, V2D);

        MassInv = 0.0*jj;
        MassInv = sum(V2D(ii,kk)*V2D(jj,kk), kk);

        // Multiply by inverse mass matrix;
        Liftref = sum(MassInv(ii,kk)*E(kk,jj),kk);
    }

    void QuadNodesProvisioner::computeInterpMatrix(const real_vector_type& rout, const real_vector_type& sout, real_matrix_type& IM) const {
        firstIndex ii;
        secondIndex jj;
        thirdIndex kk;

        index_type length = rout.size();

        real_matrix_type Vout(length, NumLocalPoints);
        computeVandermondeMatrix(NOrder, rout, sout, Vout);

        real_matrix_type& invV = *Vinv;
        IM = sum(Vout(ii,kk)*invV(kk,jj),kk);
    }

    DGContext2D QuadNodesProvisioner::get_DGContext() const {
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
