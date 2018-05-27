// Copyright (C) 2017-2018  Waterloo Quantitative Consulting Group, Inc.
// See COPYING and LICENSE files at project root for more details.

#include "DirectSolver.hpp"
#include "DenseMatrixHelpers.hpp"
#include <blitz/array.h>
#include <string>
#include <stdexcept>
#include <iomanip>

using blitz::firstIndex;
using blitz::secondIndex;
using std::runtime_error;
using std::stringstream;

namespace blitzdg {
    extern "C" {
        void dsgesv_( int* n, int* nrhs, double* a, int* lda,
                    int* ipiv, double* b, int* ldb, double* x, int* ldx, 
                    double* work, float* swork, int* iter, int* info );
    }

    void DirectSolver::solve(const real_matrix_type& A, const real_matrix_type& B, real_matrix_type& X) const {

        firstIndex ii;
        secondIndex jj;

        index_type sz = A.rows();
        index_type Nrhs = B.cols();

        index_type dim = sz*Nrhs;

        index_type lda = sz;
        index_type ldb = sz; 
        index_type ldx = sz;

        index_type ipiv[sz];

        real_type work[sz*Nrhs];
        float swork[sz*(sz+Nrhs)];

        index_type info;
        index_type iter;

        real_type Apod[sz*lda];
        real_type Bpod[dim];
        real_type Xpod[dim];

        real_matrix_type Atrans(sz, sz);
        real_matrix_type Btrans(Nrhs, sz);
        real_matrix_type Xtrans(Nrhs, sz);

        Atrans = A(jj,ii);
        Btrans = B(jj,ii);


        fullToPodArray(Atrans, Apod);
        fullToPodArray(Btrans, Bpod);

        dsgesv_(&sz, &Nrhs, Apod, &lda,
                ipiv, Bpod, &ldb, Xpod, &ldx, 
                work, swork, &iter, &info);

        podArrayToFull(Xpod, Xtrans);

        X = Xtrans(jj,ii);
    }
} // namespace blitzdg