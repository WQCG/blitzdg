// Copyright (C) 2017-2022  Waterloo Quantitative Consulting Group, Inc.
// See COPYING and LICENSE files at project root for more details.

#include "JacobiBuilders.hpp"
#include <cmath>
#include <limits>

using blitz::firstIndex;
using blitz::Range;
using blitz::secondIndex;
using std::numeric_limits;
using std::pow;
using std::sqrt;
using std::tgamma;

namespace blitzdg {

    void JacobiBuilders::computeJacobiPolynomial(real_vector_type const & x, real_type alpha, real_type beta, index_type N, real_vector_type & p) const {
        Range all = Range::all();
        index_type Np = (x.length())(0);

        real_matrix_type pStorage(N+1, Np);
        pStorage = 0.0;

        real_type gamma0 = pow(2,(alpha+beta+1))/(alpha+beta+1)*tgamma(alpha+1)*tgamma(beta+1)/tgamma(alpha+beta+1);

        p = 1/sqrt(gamma0);
        pStorage(0, all) = p;

        if (N==0) 
            return;

        real_type gamma1 = (alpha+1)*(beta+1)/(alpha+beta+3)*gamma0;
        p = ((alpha+beta+2)*x/2 + (alpha-beta)/2)/sqrt(gamma1);

        pStorage(1, all) = p;

        if (N==1) 
            return;

        real_type aold = 2/(2+alpha+beta)*sqrt((alpha+1)*(beta+1)/(alpha+beta+3));

        // Forward recurrence using the symmetry of the recurrence.
        for(index_type i=1; i <= N-1; i++) {
            real_type h1 = 2*i+alpha+beta;
            real_type anew = 2/(h1+2)*sqrt( (i+1)*(i+1+alpha+beta)*(i+1+alpha)*(i+1+beta)/(h1+1)/(h1+3));
            real_type bnew = - (alpha*alpha-beta*beta)/h1/(h1+2);
            pStorage(i+1,all) = 1/anew*( -aold*pStorage(i-1,all) + (x-bnew)*pStorage(i,all));
            aold = anew;
        }
        p = pStorage(N, all);
    }

    void JacobiBuilders::computeJacobiQuadWeights(real_type alpha, real_type beta, index_type N, real_vector_type& x, real_vector_type& w) const {

        if ( N == 0) {
            x(0) = -(alpha-beta)/(alpha+beta+2);
            w(0) = 2.0;
            return;
        }

        const real_type eps = numeric_limits<real_type>::epsilon();

        firstIndex ii;
        secondIndex jj;

        // Form symmetric matrix.
        real_matrix_type J(N+1,N+1);
        J = 0.;
        for (index_type i=0; i < N+1; i++) {
            real_type h1 = 2.*i+alpha+beta;
            J(i,i)   = -0.5*(alpha*alpha-beta*beta)/(h1+2.)/h1;
            if (i < N) {
                J(i,i+1) = 2./(h1+2)*sqrt((i+1)*((i+1)+alpha+beta)*((i+1)+alpha)*((i+1)+beta)/(h1+1)/(h1+3));
            }
        }

        if ((alpha + beta) < 10*eps ) {
            J(0,0) = 0.0;
        }
        J = J(ii,jj) + J(jj,ii);
        
        real_matrix_type eigenvectors(N+1, N+1);

        EigSolver.solve(J, x, eigenvectors);

        // The eigenvalues give the x points.
        
        // The weights are given by:
        real_vector_type v1(N+1);
        v1 = eigenvectors( 0, Range::all() ); 

        real_type gamma0 = pow(2,(alpha+beta+1))/(alpha+beta+1)*tgamma(alpha+1)*tgamma(beta+1)/tgamma(alpha+beta+1);

        w = (v1*v1)*gamma0;
    }

    void JacobiBuilders::computeGaussLobottoPoints(real_type alpha, real_type beta, index_type N, real_vector_type& x) const {
        if (N==1) {
            x(0) = -1.0;
            x(1) = 1.0;
            return;
        }

        x(0) = -1.0;
        x(N) = 1.0;

        real_vector_type xJG(N-1);
        real_vector_type w(N-1);

        computeJacobiQuadWeights(alpha+1., beta+1., N-2, xJG, w);
        
        for(index_type i=1; i < N; i++)
            x(i) = xJG(i-1);
    }

    void JacobiBuilders::computeGradJacobi(real_vector_type const & x, real_type alpha, real_type beta, index_type N, real_vector_type& dp) const {
        if (N == 0) {
            dp = 0.0;
            return;
        }

        real_vector_type p(x.length());

        computeJacobiPolynomial(x, alpha+1, beta+1, N-1, p);
        dp = sqrt(N*(N+alpha+beta+1))*p;
    }
}