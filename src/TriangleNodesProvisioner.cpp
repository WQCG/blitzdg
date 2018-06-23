// Copyright (C) 2017-2018  Waterloo Quantitative Consulting Group, Inc.
// See COPYING and LICENSE files at project root for more details.

#include "TriangleNodesProvisioner.hpp"
#include "CSCMatrix.hpp"
#include "BlitzHelpers.hpp"
#include "DenseMatrixInverter.hpp"
#include "MeshManager.hpp"
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

namespace blitzdg {
    const index_type TriangleNodesProvisioner::NumFaces = 3;
    const real_type TriangleNodesProvisioner::NodeTol = 1.e-5;

    TriangleNodesProvisioner::TriangleNodesProvisioner(index_type _NOrder, index_type _NumElements, const MeshManager * _MeshManager) 
        : NumElements{ _NumElements }, NOrder{ _NOrder },
        NumLocalPoints{ (int)0.5*(_NOrder + 2)*(_NOrder+1) }, 
        xGrid{ new real_matrix_type((int)0.5*(_NOrder + 2)*(_NOrder+1), NumElements) },
        rGrid{ new real_vector_type((int)0.5*(_NOrder + 2)*(_NOrder+1)) },
        V{ new real_matrix_type((int)0.5*(_NOrder + 2)*(_NOrder+1), (int)0.5*(_NOrder + 2)*(_NOrder+1)) }, 
        Dr{ new real_matrix_type((int)0.5*(_NOrder + 2)*(_NOrder+1), (int)0.5*(_NOrder + 2)*(_NOrder+1)) },
        Lift{ new real_matrix_type((int)0.5*(_NOrder + 2)*(_NOrder+1), (_NOrder+1)*NumFaces) },
        J{ new real_matrix_type((int)0.5*(_NOrder + 2)*(_NOrder+1), NumElements) },
        rx{ new real_matrix_type((int)0.5*(_NOrder + 2)*(_NOrder+1), NumElements) },
        nx{ new real_matrix_type((_NOrder+1)*NumFaces, NumElements) },
        Vinv{ new real_matrix_type((int)0.5*(_NOrder + 2)*(_NOrder+1), (int)0.5*(_NOrder + 2)*(_NOrder+1)) },
        Fmask{ new index_vector_type((_NOrder+1)*NumFaces) },
        Fx{ new real_matrix_type((_NOrder+1)*NumFaces, NumElements) },
        Fscale{ new real_matrix_type((_NOrder+1)*NumFaces, NumElements) },
        EToV{ new index_matrix_type(NumElements, NumFaces) },
        EToE{ new index_matrix_type(NumElements, NumFaces) },
        EToF{ new index_matrix_type(NumElements, NumFaces) },
        vmapM{ new index_vector_type((_NOrder+1)*NumFaces*NumElements) },
        vmapP{ new index_vector_type((_NOrder+1)*NumFaces*NumElements) },
        Mesh2D { _MeshManager },
        Nodes1D{ new Nodes1DProvisioner(_NOrder, NumElements, -1.0, 1.0) },
		Jacobi{}, Vandermonde{}, LinSolver{}
    {}


	void TriangleNodesProvisioner::evaluateSimplexPolynomial(const real_vector_type & a, const real_vector_type & b, const index_type i, const index_type j, real_vector_type & p) {
		real_vector_type h1(a.length(0));
		real_vector_type h2(b.length(0));

		Jacobi.computeJacobiPolynomial(a, 0.0, 0.0, i, h1);
		Jacobi.computeJacobiPolynomial(b, 2.0*i+1.0, 0.0, j, h2);

		real_vector_type c(b.length(0));

		p = sqrt(2.0)*h1*h2*pow(1.-b, i);
	}

    void TriangleNodesProvisioner::rsToab(const real_vector_type & r, const real_vector_type & s, real_vector_type & a, real_vector_type & b) {
        int Np = r.length(0);
        for( index_type i=0; i < Np; i++) {
            if(s(i) != 1.0) 
                a(i) = 2.0*(1.0+r(i))/(1.0-s(i)) - 1.0;
            else
                a(i) = -1.0;
        }
        b = s;
    }
    /*
    % function [V2D] = Vandermonde2D(N, r, s);
    % Purpose : Initialize the 2D Vandermonde Matrix,  V_{ij} = phi_j(r_i, s_i);

    V2D = zeros(length(r),(N+1)*(N+2)/2);

    % Transfer to (a,b) coordinates
    [a, b] = rstoab(r, s);

    % build the Vandermonde matrix
    sk = 1;
    for i=0:N
    for j=0:N - i
        V2D(:,sk) = Simplex2DP(a,b,i,j);
        sk = sk+1;
    end
    end
    return;
    */
    void TriangleNodesProvisioner::buildVandermondeMatrix() {
    }

    void TriangleNodesProvisioner::buildNodes() {

    }

    void TriangleNodesProvisioner::computeWarpFactor(const real_vector_type & r, real_vector_type & warpFactor) const {
        firstIndex ii;
        secondIndex jj;

        const int Np = NOrder+1;
        const int Nr = r.length(0);

        real_vector_type req(Np), rLGL(Np);
        req = -1.0 + 2*ii/(Np-1.0);

        Jacobi.computeGaussLobottoPoints(0.0, 0.0, NOrder, rLGL);

        real_matrix_type Veq(Np,Np), Veqinv(Np,Np);
        Vandermonde.computeVandermondeMatrix(req, Veq, Veqinv);

        // Evaluate Lagrange polynomial at rout
    
        real_matrix_type P(Np, Nr), L(Np,Nr);
        for (index_type i=0; i< Np; i++) {
            real_vector_type p(Nr);
            Jacobi.computeJacobiPolynomial(r, 0.0, 0.0, i, p);
            P(i, Range::all()) = p;
        }

        std::cout << "P: " << P << std::endl;

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
}