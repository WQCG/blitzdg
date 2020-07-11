// Copyright (C) 2017-2020  Waterloo Quantitative Consulting Group, Inc.
// See COPYING and LICENSE files at project root for more details.

#include "DirectSolver.hpp"
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
        void dsgesv_( int* n, int* nrhs, double* a, int* lda,
                    int* ipiv, double* b, int* ldb, double* x, int* ldx, 
                    double* work, float* swork, int* iter, int* info );
    }

    void DirectSolver::solve(const real_matrix_type& A, const real_matrix_type& B, real_matrix_type& X) const {

        index_type sz = A.rows();
        index_type Nrhs = B.cols();

        index_type dim = sz*Nrhs;

        index_type lda = sz;
        index_type ldb = sz; 
        index_type ldx = sz;

        unique_ptr<index_type[]> ipiv(new index_type[sz]());
        unique_ptr<real_type[]> work(new real_type[sz*Nrhs]());
        unique_ptr<float[]> swork(new float[sz*(sz+Nrhs)]());

        unique_ptr<real_type[]> Apod(new real_type[sz*lda]());
        unique_ptr<real_type[]> Bpod(new real_type[dim]());
        unique_ptr<real_type[]> Xpod(new real_type[dim]());

        index_type info;
        index_type iter;

        reshapeMatTo1D(A, Apod.get(), false);
        reshapeMatTo1D(B, Bpod.get(), false);

        dsgesv_(&sz, &Nrhs, Apod.get(), &lda,
                ipiv.get(), Bpod.get(), &ldb, Xpod.get(), &ldx, 
                work.get(), swork.get(), &iter, &info);

        stringstream strm;
        if (info < 0) {
            strm << "Error calling DSGESV. Error was in Argument " << info*(-1) << "." << endl;
            throw runtime_error(strm.str());
        } else if (info > 0) {
            strm << "Solution is singular. Factor U contains a diagonal element U(i,i) that is exactly zero, with i=" << info << "." << endl;
            throw runtime_error(strm.str());
        }

        reshape1DToMat(Xpod.get(), X, false);
    }
} // namespace blitzdg