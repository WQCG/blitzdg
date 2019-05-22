// Copyright (C) 2017-2019  Waterloo Quantitative Consulting Group, Inc.
// See COPYING and LICENSE files at project root for more details.

#include "Nodes1DProvisioner.hpp"
#include "CSCMatrix.hpp"
#include "BlitzHelpers.hpp"
#include "VandermondeBuilders.hpp"
#include "Types.hpp"
#include <blitz/array.h>
#include <cmath>
#include <limits>
#include <memory>
#include <boost/python/numpy.hpp>

using blitz::firstIndex;
using blitz::Range;
using blitz::secondIndex;
using blitz::sum;
using blitz::thirdIndex;
using std::numeric_limits;
using std::unique_ptr;
using boost::python::numpy::ndarray;
using boost::python::numpy::zeros;
using boost::python::numpy::dtype;



namespace blitzdg {
    const index_type Nodes1DProvisioner::NumFacePoints = 1;
    const index_type Nodes1DProvisioner::NumFaces = 2;
    const real_type Nodes1DProvisioner::NodeTol = 1.e-5;

    Nodes1DProvisioner::Nodes1DProvisioner(index_type _NOrder, index_type _NumElements, real_type _xmin, real_type _xmax) 
        : Min_x{ _xmin }, Max_x{ _xmax }, NumElements{ _NumElements }, NOrder{ _NOrder },
        NumLocalPoints{ _NOrder + 1 }, mapI{ 0 }, mapO{ NumFacePoints*NumFaces*NumElements - 1 },
        vmapI{ 0 }, vmapO{ (_NOrder + 1)*NumElements - 1 },
        xGrid{ new real_matrix_type(_NOrder + 1, NumElements) },
        rGrid{ new real_vector_type(_NOrder + 1) },
        V{ new real_matrix_type(_NOrder + 1, _NOrder + 1) }, 
        Dr{ new real_matrix_type(_NOrder + 1, _NOrder + 1) },
        Lift{ new real_matrix_type(_NOrder + 1, NumFacePoints*NumFaces) },
        J{ new real_matrix_type(_NOrder + 1, NumElements) },
        rx{ new real_matrix_type(_NOrder + 1, NumElements) },
        nx{ new real_matrix_type(NumFacePoints*NumFaces, NumElements) },
        Vinv{ new real_matrix_type(_NOrder + 1, _NOrder + 1) },
        Fmask{ new index_vector_type(NumFacePoints*NumFaces) },
        Fx{ new real_matrix_type(NumFacePoints*NumFaces, NumElements) },
        Fscale{ new real_matrix_type(NumFacePoints*NumFaces, NumElements) },
        EToV{ new index_matrix_type(NumElements, NumFaces) },
        EToE{ new index_matrix_type(NumElements, NumFaces) },
        EToF{ new index_matrix_type(NumElements, NumFaces) },
        vmapM{ new index_vector_type(NumFacePoints*NumFaces*NumElements) },
        vmapP{ new index_vector_type(NumFacePoints*NumFaces*NumElements) },
        LinSolver{}, Jacobi{}, Vandermonde{}
    {}

    void Nodes1DProvisioner::buildNodes() {
        const real_type alpha = 0.0;
        const real_type beta = 0.0;

        real_vector_type& r = *rGrid;

        Jacobi.computeGaussLobottoPoints(alpha, beta, NOrder, r);

        Vandermonde.computeVandermondeMatrix(r, *V, *Vinv);
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

        real_matrix_type xmat(NumLocalPoints, NumElements, ColumnMajorOrder());
        xmat = *xGrid;

        unique_ptr<real_type[]> x(new real_type[NumElements*NumLocalPoints]());
        reshapeMatTo1D(xmat, x.get(), false);

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
        index_type totalFaces = NumFaces * NumElements;
        index_type numVertices = NumElements + 1;
        index_type localVertNum[2] = {0, 1};

        // Build global face-to-vertex array. Note: we actually
        // build its transpose for easy column access.
        CSCMat FToVtrans(numVertices, totalFaces, totalFaces);
        const index_matrix_type& E2V = *EToV;
        index_type globalFaceNum = 0;
        for (index_type k = 0; k < NumElements; ++k) {
            for (index_type f = 0; f < NumFaces; ++f) {
                index_type v = localVertNum[f];
                index_type vGlobal = E2V(k, v);
                FToVtrans.colPtrs(globalFaceNum) = globalFaceNum;
                FToVtrans.rowInds(globalFaceNum) = vGlobal;
                FToVtrans.elems(globalFaceNum) = 1.0;
                ++globalFaceNum;
            }
        }
        FToVtrans.colPtrs(totalFaces) = globalFaceNum;
        CSCMat FToF = multiply(transpose(FToVtrans), FToVtrans);
        
        index_vector_type f1(totalFaces - 2); // '- 2' => for physical boundaries.
        index_vector_type f2(totalFaces - 2);
        f1 = 0;
        f2 = 0;
        
        index_type connectionsCount = 0;
        for (index_type j = 0; j < totalFaces; ++j) {
            for (index_type k = FToF.colPtrs(j); k < FToF.colPtrs(j + 1); ++k) {
                index_type i = FToF.rowInds(k);
                if (i != j && FToF.elems(k) == 1.0) {
                    f1(connectionsCount) = i;
                    f2(connectionsCount) = j;
                    ++connectionsCount;
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
        index_matrix_type& E2E = *EToE;
        index_matrix_type& E2F = *EToF;
        for (index_type k = 0; k < NumElements; ++k) {
            for (index_type f = 0; f < NumFaces; ++f) {
                E2E(k, f) = k;
                E2F(k, f) = f;
            }
        }

        for (index_type i = 0; i < totalFaces - 2; ++i) {
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

    void Nodes1DProvisioner::buildDr() {
        firstIndex ii;
        secondIndex jj;

        real_matrix_type & Vref = *V;
        real_matrix_type & Drref = *Dr;

        real_matrix_type DVr(NOrder+1, NOrder+1);
        DVr = 0.*jj;

        Vandermonde.computeGradVandermonde(*rGrid, DVr);

		// Note: this is not a column major ordering trick. We need these transposes.
        real_matrix_type Vtrans(NOrder+1, NOrder+1);
        real_matrix_type DVrtrans(NOrder+1, NOrder+1);
        real_matrix_type Drtrans(NOrder+1, NOrder+1);

        Vtrans = Vref(jj, ii);
        DVrtrans = DVr(jj, ii);

        // Dr = DVr / V;
        LinSolver.solve(Vtrans, DVrtrans,  Drtrans);

        Drref = Drtrans(jj, ii);
    } 

    index_type Nodes1DProvisioner::get_NumElements() const {
        return NumElements;
    }

    const real_matrix_type & Nodes1DProvisioner::get_xGrid() const {
        return *xGrid;
    }

    ndarray Nodes1DProvisioner::get_xGrid_numpy() const {
        Py_intptr_t shape[2] = { NumLocalPoints, NumElements };
        ndarray result = zeros(2, shape, dtype::get_builtin<real_type>());
        std::copy(xGrid->begin(), xGrid->end(), reinterpret_cast<real_type*>(result.get_data()));
        return result;
    }

    ndarray Nodes1DProvisioner::get_Dr_numpy() const {
        Py_intptr_t shape[2] = { NumLocalPoints, NumLocalPoints };
        ndarray result = zeros(2, shape, dtype::get_builtin<real_type>());
        std::copy(Dr->begin(), Dr->end(), reinterpret_cast<real_type*>(result.get_data()));
        return result;
    }

    ndarray Nodes1DProvisioner::get_Fscale_numpy() const {
        Py_intptr_t shape[2] = { NumFacePoints*NumFaces, NumElements };
        ndarray result = zeros(2, shape, dtype::get_builtin<real_type>());
        std::copy(Fscale->begin(), Fscale->end(), reinterpret_cast<real_type*>(result.get_data()));
        return result;
    }

    ndarray Nodes1DProvisioner::get_rx_numpy() const {
        Py_intptr_t shape[2] = { NumLocalPoints, NumElements };
        ndarray result = zeros(2, shape, dtype::get_builtin<real_type>());
        std::copy(rx->begin(), rx->end(), reinterpret_cast<real_type*>(result.get_data()));
        return result;
    }

    ndarray Nodes1DProvisioner::get_Lift_numpy() const {
        Py_intptr_t shape[2] = { NumLocalPoints, NumFacePoints*NumFaces };
        ndarray result = zeros(2, shape, dtype::get_builtin<real_type>());
        std::copy(Lift->begin(), Lift->end(), reinterpret_cast<real_type*>(result.get_data()));
        return result;
    }

    ndarray Nodes1DProvisioner::get_vmapM_numpy() const {
        Py_intptr_t shape[1] = { NumElements*NumFacePoints*NumFaces };
        ndarray result = zeros(1, shape, dtype::get_builtin<index_type>());
        std::copy(vmapM->begin(), vmapM->end(), reinterpret_cast<index_type*>(result.get_data()));
        return result;
    }

    ndarray Nodes1DProvisioner::get_vmapP_numpy() const {
        Py_intptr_t shape[1] = { NumElements*NumFacePoints*NumFaces };
        ndarray result = zeros(1, shape, dtype::get_builtin<index_type>());
        std::copy(vmapP->begin(), vmapP->end(), reinterpret_cast<index_type*>(result.get_data()));
        return result;
    }

    ndarray Nodes1DProvisioner::get_nx_numpy() const {
        Py_intptr_t shape[2] = { NumFacePoints*NumFaces, NumElements };
        ndarray result = zeros(2, shape, dtype::get_builtin<real_type>());
        std::copy(nx->begin(), nx->end(), reinterpret_cast<real_type*>(result.get_data()));
        return result;
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

    const real_matrix_type & Nodes1DProvisioner::get_Vinv() const {
        return *Vinv;
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

} // namespace blitzdg