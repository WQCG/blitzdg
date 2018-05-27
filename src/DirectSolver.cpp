// Copyright (C) 2017-2018  Waterloo Quantitative Consulting Group, Inc.
// See COPYING and LICENSE files at project root for more details.

#include "DirectSolver.hpp"
#include "DenseMatrixHelpers.hpp"
#include "Types.hpp"
#include <string>
#include <stdexcept>
#include <iomanip>

using std::runtime_error;
using std::stringstream;
using std::endl;

namespace blitzdg {
    extern "C" {
        void dsgesv_( int* n, int* nrhs, double* a, int* lda,
                    int* ipiv, double* b, int* ldb, double* x, int* ldx, 
                    double* work, float* swork, int* iter, int* info );
    }

    void DirectSolver::solve(const real_matrix_type& A, const real_matrix_type& B, real_matrix_type& X) const {

        if ( !(isRowMajor(A) == isRowMajor(B) ) )
            throw runtime_error("Matrices A and B must have same storage order!");

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

        fullToPodArray(A, Apod, false);
        fullToPodArray(B, Bpod, false);

        dsgesv_(&sz, &Nrhs, Apod, &lda,
                ipiv, Bpod, &ldb, Xpod, &ldx, 
                work, swork, &iter, &info);

        stringstream strm;
        if (info < 0) {
            strm << "Error calling DSGESV. Error was in Argument " << info*(-1) << "." << endl;
            throw runtime_error(strm.str());
        } else if (info > 0) {
            strm << "Solution is singular. Factor U contains a diagonal element U(i,i) that is exactly zero, with i=" << info << "." << endl;
            throw runtime_error(strm.str());
        }

        podArrayToFull(Xpod, X, false);
    }
} // namespace blitzdg