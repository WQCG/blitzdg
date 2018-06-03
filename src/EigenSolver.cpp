// Copyright (C) 2017-2018  Waterloo Quantitative Consulting Group, Inc.
// See COPYING and LICENSE files at project root for more details.

#include "EigenSolver.hpp"
#include "BlitzHelpers.hpp"
#include <blitz/array.h>
#include <iomanip>
#include <stdexcept>
#include <memory>

using blitz::firstIndex;
using blitz::secondIndex;
using std::stringstream;
using std::endl;
using std::runtime_error;
using std::unique_ptr;

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

        unique_ptr<real_type[]> Apod(new real_type[sz*lda]());

        reshapeMatTo1D(A, Apod.get());

        /* Determining optimal workspace parameters */
        dsyevd_( &JOBZ, UPLO, &sz, Apod.get(), &lda, ww, &wkopt, &lwork, &iwkopt, &liwork, &info );
        stringstream strm;
        if (info < 0) {
            strm << "Error calling DSYEVD to determine workspace parameters. Error was in Argument " << info*(-1) << "." << endl;
            throw runtime_error(strm.str());
        } else if (info  > 0) {
            strm << "Error calling DSYEVD to determine workspace parameters. Error code: " << info << "." << endl;
            throw runtime_error(strm.str());
        }

        lwork = static_cast<index_type>(wkopt);
        unique_ptr<real_type[]> work(new real_type[lwork]());
        liwork = iwkopt;
        unique_ptr<index_type[]> iwork(new index_type[liwork]());

        /* Solve eigenproblem */
        dsyevd_( &JOBZ, UPLO, &sz, Apod.get(), &lda, ww, work.get(), &lwork, iwork.get(), &liwork, &info );
        if (info < 0) {
            strm << "Error calling DSYEVD. Error was in Argument " << info*(-1) << "." << endl;
            throw runtime_error(strm.str());
        } else if (info > 0) {
            strm << "The algorithm failed to converge; i off-diagonal elements of an intermediate tridiagonal form did not converge to zero. i=" << info << "." << endl;
            throw runtime_error(strm.str());
        }

        reshape1DToMat(Apod.get(), eigenvectors, false);

        for (index_type i=0; i < sz; i++)
            eigenvalues(i) = ww[i];

    }
} // namespace blitzdg

