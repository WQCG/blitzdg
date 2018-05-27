// Copyright (C) 2017-2018  Waterloo Quantitative Consulting Group, Inc.
// See COPYING and LICENSE files at project root for more details.

#include "Nodes1DProvisioner.hpp"
#include "DenseMatrixHelpers.hpp"
#include <blitz/array.h>
#include <cmath>
#include <limits>

using blitz::firstIndex;
using blitz::Range;
using blitz::secondIndex;
using blitz::sum;
using blitz::thirdIndex;
using std::numeric_limits;

namespace blitzdg {
    const index_type Nodes1DProvisioner::NumFacePoints = 1;
    const index_type Nodes1DProvisioner::NumFaces = 2;
    const real_type Nodes1DProvisioner::NodeTol = 1.e-5;

    Nodes1DProvisioner::Nodes1DProvisioner(index_type _NOrder, index_type _NumElements, real_type _xmin, real_type _xmax) 
        : Min_x{ _xmin }, Max_x{ _xmax }, NumElements{ _NumElements }, NOrder{ _NOrder },
        NumLocalPoints{ NOrder + 1 }, mapI{ 0 }, mapO{ NumFacePoints*NumFaces*NumElements - 1 },
        vmapI{ 0 }, vmapO{ NumLocalPoints*NumElements - 1 },
        xGrid{ new real_matrix_type(NumLocalPoints, NumElements) },
        rGrid{ new real_vector_type(NumLocalPoints) },
        V{ nullptr }, Dr{ nullptr }, 
        Lift{ new real_matrix_type(NumLocalPoints, NumFacePoints*NumFaces) },
        J{ new real_matrix_type(NumLocalPoints, NumElements) },
        rx{ new real_matrix_type(NumLocalPoints, NumElements) },
        nx{ new real_matrix_type(NumFacePoints*NumFaces, NumElements) },
        Fmask{ new index_vector_type (NumFacePoints*NumFaces) },
        Fx{ new real_matrix_type(NumFacePoints*NumFaces, NumElements) },
        Fscale{ new real_matrix_type(NumFacePoints*NumFaces, NumElements) },
        EToV{ new index_matrix_type(NumElements, NumFaces) },
        EToE{ new index_matrix_type(NumElements, NumFaces) },
        EToF{ new index_matrix_type(NumElements, NumFaces) },
        vmapM{ new index_vector_type(NumFacePoints*NumFaces*NumElements) },
        vmapP{ new index_vector_type(NumFacePoints*NumFaces*NumElements) },
        EigSolver{}, LinSolver{}, Jacobi{}
    {}

    void Nodes1DProvisioner::buildNodes() {
        const real_type alpha = 0.0;
        const real_type beta = 0.0;

        real_vector_type& r = *rGrid;

        Jacobi.computeGaussLobottoPoints(alpha, beta, NOrder, r);
        
        buildVandermondeMatrix();
        buildDr();
        buildLift();

        real_type L = Max_x - Min_x;
        real_type width = L / NumElements;

        real_matrix_type & x = *xGrid;
        for (index_type k=0; k < NumElements; k++) {
            x(Range::all(), k) = Min_x + width*(k + 0.5*(r+1.));
        }

        index_matrix_type & E2V = *EToV;

        // Create Element-to-Vertex connectivity table.
        for (index_type k=0; k < NumElements; k++) {
            E2V(k, 0) = k;
            E2V(k, 1) = k+1;
        }

        buildConnectivityMatrices();
        buildFaceMask();
        buildMaps();
        buildNormals();
    }

    void Nodes1DProvisioner::buildNormals() {
        real_matrix_type & nxref = *nx;

        real_type mult = -1;

        for (index_type k=0; k < NumElements; k++) {
            for (index_type f=0; f < NumFaces*NumFacePoints; f++) {
                nxref(f,k) = mult;
                mult *= -1;
            }
        }
    }

    void Nodes1DProvisioner::buildMaps() {
        firstIndex ii;
        secondIndex jj;

        index_matrix_type nodeIds(NumLocalPoints, NumElements);
        
        index_vector_type & Fmsk = *Fmask;
        index_vector_type & vmM = *vmapM;
        index_vector_type & vmP = *vmapP;

        index_matrix_type & E2E = *EToE;
        index_matrix_type & E2F = *EToF;
        real_matrix_type & xmat = *xGrid;

        real_matrix_type xmatTrans(NumElements, NumLocalPoints);
        xmatTrans = xmat(jj,ii);

        real_type * x = new real_type[NumElements*NumLocalPoints];
        fullToPodArray(xmatTrans, x);

        // Assemble global volume node numbering.
        nodeIds = ii + NumLocalPoints*jj;

        vmM = 0*ii;
        vmP = 0*ii;

        index_type count=0;
        for (index_type k = 0; k < NumElements; k++) {
            for( index_type f = 0; f < NumFaces; f++ ) {
                vmM(count) = nodeIds(Fmsk(f), k);
                count++;
            }
        }

        count = 0;
        for (index_type k1=0; k1 < NumElements; k1++) {
            for (index_type f1=0; f1 < NumFaces; f1++) {
                index_type k2 = E2E(k1, f1);
                index_type f2 = E2F(k1, f1);

                index_type vidM = vmM(k1*NumFaces + f1);
                index_type vidP = vmM(k2*NumFaces + f2);

                real_type dx = x[vidM] - x[vidP];
                real_type dist = sqrt(dx * dx);

                if ( dist < NodeTol) {
                    vmP(count) = vidP;
                }
                count++;
            }
        }

        delete[] x;
    }

    void Nodes1DProvisioner::buildFaceMask() {
        real_matrix_type & x = *xGrid;
        real_matrix_type & Fxref = *Fx;
        index_vector_type & Fmaskref = *Fmask;

        Fmaskref = 0, (NumLocalPoints - 1);

        for (index_type k = 0;  k < NumElements; k++) {
            for (index_type f = 0; f < NumFacePoints*NumFaces; f++) {
                Fxref(f, k) = x(Fmaskref(f), k);
            }
        }
    }

    void Nodes1DProvisioner::buildConnectivityMatrices() {
        firstIndex ii;
        secondIndex jj;
        thirdIndex kk;

        index_type totalFaces = NumFaces*NumElements;
        index_type numVertices = NumElements + 1;

        index_type localVertNum[2];
        localVertNum[0] = 0; localVertNum[1] = 1;

        // Build global face-to-vertex array. (should be sparse matrix in 2D/3D).
        real_matrix_type FToV(totalFaces, numVertices);
        FToV = 0*jj;

        index_matrix_type & E2V = *EToV;

        index_type globalFaceNum = 0;
        for (index_type k=0; k < NumElements; k++) {
            for (index_type f=0; f < NumFaces; f++) {
                index_type v = localVertNum[f];
                index_type vGlobal = E2V(k,v);
                FToV(globalFaceNum, vGlobal) = 1;
                globalFaceNum++;
            }
        }

        real_matrix_type FToF(totalFaces, totalFaces);
        real_matrix_type I(totalFaces, totalFaces);

        FToF = 0*jj;
        I = 0*jj;

        for (index_type f=0; f < totalFaces; f++)
            I(f,f) = 1;

        // Global Face-to-Face connectivity matrix.
        FToF = sum(FToV(ii,kk)*FToV(jj,kk), kk) - I;

        index_vector_type f1(totalFaces - 2); // '- 2' => for physical boundaries.
        index_vector_type f2(totalFaces - 2);

        f1 = 0*ii;
        f2 = 0*ii;

        index_type connectionsCount = 0;
        for (index_type i=0; i < totalFaces; i++) {
            for (index_type j=0; j < totalFaces; j++) {
                if (FToF(i,j) == 1) {
                    f1(connectionsCount) = i;
                    f2(connectionsCount) = j;
                    connectionsCount++;
                }
            }
        }

        index_vector_type e1(totalFaces - 2);
        index_vector_type e2(totalFaces - 2);

        // Convert face global number to local element and face numbers.
        e1 = floor(f1 / NumFaces);
        f1 = (f1 % NumFaces);
        e2 = floor(f2 / NumFaces);
        f2 = (f2 % NumFaces);

        // Build connectivity matrices.
        index_matrix_type & E2E = *EToE;
        index_matrix_type & E2F = *EToF;
        for (index_type k = 0; k < NumElements; k++) {
            for (index_type f = 0; f < NumFaces; f++) {
                E2E(k, f) = k;
                E2F(k, f) = f;
            }
        }

        for (index_type i=0; i < totalFaces - 2; i++) {
            index_type ee1 = e1(i);
            index_type ee2 = e2(i);
            index_type ff1 = f1(i);
            index_type ff2 = f2(i);
            E2E(ee1, ff1) = ee2;
            E2F(ee1, ff1) = ff2;
        }
    }

    void Nodes1DProvisioner::buildLift() {
        index_type Np = NumLocalPoints;
        firstIndex ii;
        secondIndex jj;
        thirdIndex kk;

        real_matrix_type E(Np, NumFaces*NumFacePoints);
        E = 0*jj;
        E(0, 0)  = 1.;
        E(Np-1, 1) = 1.;

        real_matrix_type & Vref = *V;
        real_matrix_type & L = *Lift;

        real_matrix_type Vtrans(Np, Np);
        Vtrans = Vref(jj,ii);
        real_matrix_type temp(Np, NumFaces*NumFacePoints);
        temp = sum(Vtrans(ii,kk)*E(kk,jj), kk);
        L = sum(Vref(ii,kk)*temp(kk,jj), kk);
    }

    void Nodes1DProvisioner::computeJacobian() {
        firstIndex ii;
        secondIndex jj;
        thirdIndex kk;

        real_matrix_type & x = *xGrid;
        real_matrix_type & Drref = *Dr;
        real_matrix_type & Jref = *J;
        real_matrix_type & rxref = *rx;

        real_matrix_type & Fscaleref = *Fscale;

        index_vector_type & Fmsk = *Fmask;

        Jref = sum(Drref(ii,kk)*x(kk,jj), kk);
        rxref = 1/Jref;

        for(index_type f=0; f < NumFaces; f++) {
            Fscaleref(f, Range::all()) = 1/Jref(Fmsk(f), Range::all());
        }
    }

    void Nodes1DProvisioner::buildVandermondeMatrix() {
        V = new real_matrix_type(NOrder+1, NOrder+1);

        real_matrix_type & Vref = *V;

        real_vector_type p(NOrder+1);
        for (index_type j=0; j <= NOrder; j++) {
            Jacobi.computeJacobiPolynomial(*rGrid, 0.0, 0.0, j, p);
            Vref(Range::all(), j) = p;
        }
    }

    void Nodes1DProvisioner::buildDr() {
        firstIndex ii;
        secondIndex jj;

        Dr = new real_matrix_type(NOrder+1, NOrder+1);

        real_matrix_type & Vref = *V;
        real_matrix_type & Drref = *Dr;

        real_matrix_type DVr(NOrder+1, NOrder+1);
        DVr = 0.*jj;

        computeGradVandermonde(DVr);

        // Dr = DVr / V;

        real_matrix_type Vtrans(NOrder+1, NOrder+1);
        real_matrix_type DVrtrans(NOrder+1, NOrder+1);
        real_matrix_type Drtrans(NOrder+1, NOrder+1);

        Vtrans = Vref(jj, ii);
        DVrtrans = DVr(jj, ii);

        LinSolver.solve(Vtrans, DVrtrans,  Drtrans);

        Drref = Drtrans(jj, ii);
    } 

    index_type Nodes1DProvisioner::get_NumElements() const {
        return NumElements;
    }

    const real_matrix_type & Nodes1DProvisioner::get_xGrid() const {
        return *xGrid;
    }

    const real_vector_type & Nodes1DProvisioner::get_rGrid() const {
        return *rGrid;
    }

    const index_matrix_type & Nodes1DProvisioner::get_EToV() const {
        return *EToV;
    }

    const real_matrix_type & Nodes1DProvisioner::get_Lift() const {
        return *Lift;
    }

    const index_matrix_type & Nodes1DProvisioner::get_EToE() const {
        return *EToE;
    }

    const index_matrix_type & Nodes1DProvisioner::get_EToF() const {
        return *EToF;
    }

    const real_matrix_type & Nodes1DProvisioner::get_Dr() const {
        return *Dr;
    }

    const real_matrix_type & Nodes1DProvisioner::get_V() const {
        return *V;
    }

    const real_matrix_type & Nodes1DProvisioner::get_J() const {
        return *J;
    }

    const real_matrix_type & Nodes1DProvisioner::get_rx() const {
        return *rx;
    }

    const real_matrix_type & Nodes1DProvisioner::get_nx() const {
        return *nx;
    } 

    index_type Nodes1DProvisioner::get_NumLocalPoints() const {
        return NumLocalPoints;
    }

    const real_matrix_type & Nodes1DProvisioner::get_Fx() const {
        return *Fx;
    }

    const real_matrix_type & Nodes1DProvisioner::get_Fscale() const {
        return *Fscale;
    }

    const index_vector_type & Nodes1DProvisioner::get_Fmask() const {
        return *Fmask;
    }

    const index_vector_type & Nodes1DProvisioner::get_vmapM() const {
        return *vmapM;
    }

    const index_vector_type & Nodes1DProvisioner::get_vmapP() const {
        return *vmapP;
    }

    index_type Nodes1DProvisioner::get_mapI() const {
        return mapI;
    }

    index_type Nodes1DProvisioner::get_mapO() const {
        return mapO;
    }

    index_type Nodes1DProvisioner::get_vmapI() const {
        return vmapI;
    }

    index_type Nodes1DProvisioner::get_vmapO() const {
        return vmapO;
    }

    Nodes1DProvisioner::~Nodes1DProvisioner() {
        delete V; V = nullptr;
        delete Dr; Dr = nullptr;
        delete rGrid; rGrid = nullptr;
        delete xGrid; xGrid = nullptr;
        delete J; J = nullptr;
        delete rx; rx = nullptr;
        delete Lift; Lift = nullptr;
        delete EToV; EToV = nullptr;
        delete EToE; EToE = nullptr;
        delete EToF; EToF = nullptr;
        delete Fmask; Fmask = nullptr;
        delete Fx; Fx = nullptr;
        delete Fscale; Fscale = nullptr;
        delete nx; nx = nullptr;
        delete vmapM; vmapM = nullptr;
        delete vmapP; vmapP = nullptr;
    }


    void Nodes1DProvisioner::computeGradVandermonde(real_matrix_type & DVr) const {
        firstIndex ii;
        real_vector_type dp(NOrder+1);
        for (index_type i=0; i<=NOrder; i++) {
            dp = 0.*ii;
            Jacobi.computeGradJacobi(*rGrid, 0.0, 0.0, i, dp);
            DVr(Range::all(), i) = dp;
        }
    }
} // namespace blitzdg