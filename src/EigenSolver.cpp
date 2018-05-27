// Copyright (C) 2017-2018  Waterloo Quantitative Consulting Group, Inc.
// See COPYING and LICENSE files at project root for more details.

#include "EigenSolver.hpp"
#include "DenseMatrixHelpers.hpp"
#include <blitz/array.h>

using blitz::firstIndex;
using blitz::secondIndex;

namespace blitzdg {
    extern "C" {
        void dsyevd_( char* jobz, char* uplo, int* n, double* a, int* lda,
                    double* w, double* work, int* lwork, int* iwork, int* liwork, int* info );
    }

    void EigenSolver::solve(const real_matrix_type& A, real_vector_type& eigenvalues, real_matrix_type& eigenvectors) const {
        index_type sz = A.rows();
        index_type lda = sz;
        index_type iwkopt;

        real_type ww[sz];
        real_type wkopt;
        index_type lwork = -1;
        index_type liwork = -1;
        index_type info;

        char JOBZ = 'V';
        char UPLO[] = "UP";

        real_type* Apod = new real_type[sz*lda];

        fullToPodArray(A, Apod);

        /* Determining optimal workspace parameters */
        dsyevd_( &JOBZ, UPLO, &sz, Apod, &lda, ww, &wkopt, &lwork, &iwkopt, &liwork, &info );

        lwork = static_cast<index_type>(wkopt);
        real_type* work = new real_type[lwork];
        liwork = iwkopt;
        index_type * iwork = new index_type[liwork];

        /* Solve eigenproblem */
        dsyevd_( &JOBZ, UPLO, &sz, Apod, &lda, ww, work, &lwork, iwork, &liwork, &info );

        podArrayToFull(Apod, eigenvectors);

        for (index_type i=0; i < sz; i++)
            eigenvalues(i) = ww[i];

        firstIndex ii;
        secondIndex jj;
        eigenvectors = eigenvectors(jj,ii);

        delete [] Apod;
        delete [] work;
        delete [] iwork;
    }
} // namespace blitzdg

