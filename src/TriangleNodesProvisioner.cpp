// Copyright (C) 2017-2018  Waterloo Quantitative Consulting Group, Inc.
// See COPYING and LICENSE files at project root for more details.

#include "TriangleNodesProvisioner.hpp"
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

namespace blitzdg {
    const index_type TriangleNodesProvisioner::NumFaces = 3;
    const real_type TriangleNodesProvisioner::NodeTol = 1.e-5;
    const real_type pi = blitzdg::constants::pi;

    TriangleNodesProvisioner::TriangleNodesProvisioner(index_type _NOrder, index_type _NumElements, const MeshManager * _MeshManager) 
        : NumElements{ _NumElements }, NOrder{ _NOrder },
        NumLocalPoints{ (_NOrder + 2)*(_NOrder+1)/2 }, 
        xGrid{ new real_matrix_type((_NOrder + 2)*(_NOrder+1)/2, NumElements) },
        rGrid{ new real_vector_type((_NOrder + 2)*(_NOrder+1)/2) },
        V{ new real_matrix_type((_NOrder + 2)*(_NOrder+1)/2, (_NOrder + 2)*(_NOrder+1)/2) }, 
        Dr{ new real_matrix_type((_NOrder + 2)*(_NOrder+1)/2, (_NOrder + 2)*(_NOrder+1)/2) },
        Lift{ new real_matrix_type((_NOrder + 2)*(_NOrder+1)/2, (_NOrder+1)*NumFaces) },
        J{ new real_matrix_type((_NOrder + 2)*(_NOrder+1)/2, NumElements) },
        rx{ new real_matrix_type((_NOrder + 2)*(_NOrder+1)/2, NumElements) },
        nx{ new real_matrix_type((_NOrder+1)*NumFaces, NumElements) },
        Vinv{ new real_matrix_type((_NOrder + 2)*(_NOrder+1)/2, (_NOrder + 2)*(_NOrder+1)/2) },
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


	void TriangleNodesProvisioner::evaluateSimplexPolynomial(const real_vector_type & a, const real_vector_type & b, index_type i, index_type j, real_vector_type & p) const {
		real_vector_type h1(a.length(0));
		real_vector_type h2(b.length(0));

		Jacobi.computeJacobiPolynomial(a, 0.0, 0.0, i, h1);
		Jacobi.computeJacobiPolynomial(b, 2.0*i+1.0, 0.0, j, h2);

		real_vector_type c(b.length(0));

		p = sqrt(2.0)*h1*h2*pow(1.-b, i);
	}

    void TriangleNodesProvisioner::rsToab(const real_vector_type & r, const real_vector_type & s, real_vector_type & a, real_vector_type & b) const {
        int Np = r.length(0);
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


    void TriangleNodesProvisioner::computeVandermondeMatrix(int N, const real_vector_type & r, const real_vector_type & s, const real_matrix_type & V) const {
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
        x = x + 1*warp1 + cos(2*pi/3)*warp2 + cos(4*pi/3)*warp3;
        y = y + 0*warp1 + sin(2*pi/3)*warp2 + sin(4*pi/3)*warp3;
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
}