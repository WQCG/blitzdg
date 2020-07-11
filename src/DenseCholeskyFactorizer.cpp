// Copyright (C) 2017-2020  Waterloo Quantitative Consulting Group, Inc.
// See COPYING and LICENSE files at project root for more details.

#include "DenseCholeskyFactorizer.hpp"
#include "BlitzHelpers.hpp"
#include "Types.hpp"
#include <string>
#include <stdexcept>
#include <iomanip>
#include <memory>

using std::runtime_error;
using std::stringstream;
using std::endl;
using std::unique_ptr;

namespace blitzdg {
    extern "C" {
        // Cholesky decomoposition of a a real symmetric positive definite matrix A.
        void dpotrf_(char* UPLO, int *N, double* A, int* LDA, int* INFO);
    }

    void DenseCholeskyFactorizer::computeCholesky(const real_matrix_type& A, real_matrix_type& R) const {
        // Assumes a square matrix so that A is NxN (M=N).
        index_type N = A.rows();
        index_type info;

        unique_ptr<real_type[]> Apod(new real_type[N*N]());
        char UPLO[3] = "UP";

        reshapeMatTo1D(A, Apod.get(), false);

        dpotrf_(UPLO, &N, Apod.get(), &N, &info);
        
        stringstream strm;
        if (info < 0) {
            strm << "Error calling DPOTRF. Error was in Argument " << info*(-1) << "." << endl;
            throw runtime_error(strm.str());
        } else if (info > 0) {
            // verify this message.
            strm << "The leading minor order of i is not positive definite, with i=" << info << ". The Cholesky factorization could not be completed." << endl;
            throw runtime_error(strm.str());
        }


        reshape1DToMat(Apod.get(), R, false);

        // zero out the junk below the main diagonal
        for (index_type i=1; i < R.rows(); ++i) {
            for (index_type j=0; j < i; ++j) {
                R(i,j) = 0.0;
            }
        }
    }
} // namespace blitzdg