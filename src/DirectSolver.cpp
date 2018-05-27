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
    
        real_matrix_type Acol = real_matrix_type(A.rows(), A.cols(), ColumnMajorOrder());
        real_matrix_type Bcol = real_matrix_type(B.rows(), B.cols(), ColumnMajorOrder());
        real_matrix_type Xcol = real_matrix_type(X.rows(), X.cols(), ColumnMajorOrder());
        
        Acol = A;
        Bcol = B;
        Xcol = X;

        index_type sz = Acol.rows();
        index_type Nrhs = Bcol.cols();

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

        fullToPodArray(Acol, Apod, false);
        fullToPodArray(Bcol, Bpod, false);

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

        podArrayToFull(Xpod, Xcol, false);

        X = Xcol;
    }
} // namespace blitzdg