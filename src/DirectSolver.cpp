// Copyright (C) 2017-2018  Waterloo Quantitative Consulting Group, Inc.
// See COPYING and LICENSE files at project root for more details.

#include "DirectSolver.hpp"
#include <blitz/array.h>

using blitz::firstIndex;
using blitz::secondIndex;

namespace blitzdg {
    extern "C" {
        void dsgesv_( int* n, int* nrhs, double* a, int* lda,
                    int* ipiv, double* b, int* ldb, double* x, int* ldx, 
                    double* work, float* swork, int* iter, int* info );
    }

    void DirectSolver::solve(const matrix_type& A, const matrix_type& B, matrix_type& X) const {

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

        matrix_type Atrans(sz, sz);
        matrix_type Btrans(Nrhs, sz);
        matrix_type Xtrans(Nrhs, sz);

        Atrans = A(jj,ii);
        Btrans = B(jj,ii);

        MatrixConverter.fullToPodArray(Atrans, Apod);
        MatrixConverter.fullToPodArray(Btrans, Bpod);

        dsgesv_(&sz, &Nrhs, Apod, &lda,
                ipiv, Bpod, &ldb, Xpod, &ldx, 
                work, swork, &iter, &info);

        MatrixConverter.podArrayToFull(Xpod, Xtrans);

        X = Xtrans(jj,ii);
    }
} // namespace blitzdg