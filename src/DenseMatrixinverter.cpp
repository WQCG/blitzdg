// Copyright (C) 2017-2018  Waterloo Quantitative Consulting Group, Inc.
// See COPYING and LICENSE files at project root for more details.

#include "DenseMatrixInverter.hpp"
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
        // LU decomoposition of a general matrix.
        void dgetrf_(int* M, int *N, double* A, int* lda, int* IPIV, int* INFO);

        // generate inverse of a matrix given its LU decomposition
        void dgetri_(int* N, double* A, int* lda, int* IPIV, double* WORK, int* lwork, int* INFO);

    }

    void DenseMatrixInverter::computeInverse(const real_matrix_type& A, real_matrix_type& Ainv) const {

        // Assumes a square matrix so that A is NxN (M=N).
        index_type N = A.rows();
        index_type lwork = N*N;
        index_type info;

        unique_ptr<index_type[]> ipiv(new index_type[N+1]());
        unique_ptr<real_type[]> work(new real_type[lwork]());
        unique_ptr<real_type[]> Apod(new real_type[N*N]());

        reshapeMatTo1D(A, Apod.get(), false);

        dgetrf_(&N, &N, Apod.get(), &N, ipiv.get(), &info);
        
        stringstream strm;
        if (info < 0) {
            strm << "Error calling DGETRF. Error was in Argument " << info*(-1) << "." << endl;
            throw runtime_error(strm.str());
        } else if (info > 0) {
            // verify this message.
            strm << "Solution is singular. Factor U contains a diagonal element U(i,i) that is exactly zero, with i=" << info << "." << endl;
            throw runtime_error(strm.str());
        }

        dgetri_(&N, Apod.get(), &N, ipiv.get(), work.get(), &lwork, &info);
        if (info < 0) {
            strm << "Error calling DGETRI. Error was in Argument " << info*(-1) << "." << endl;
            throw runtime_error(strm.str());
        } else if (info > 0) {
            // verify this message.
            strm << "Unable to compute inverse from LU factors with i=" << info << "." << endl;
            throw runtime_error(strm.str());
        }

        reshape1DToMat(Apod.get(), Ainv, false);
    }
} // namespace blitzdg