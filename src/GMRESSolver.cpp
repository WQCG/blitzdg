// Copyright (C) 2017-2022  Waterloo Quantitative Consulting Group, Inc.
// See COPYING and LICENSE files at project root for more details.

#include "GMRESSolver.hpp"
#include <stdexcept>

using std::invalid_argument;
using std::ostream;
using std::setprecision;
using std::scientific;
using std::string;

namespace blitzdg {
    namespace details {
        extern "C" {
            void dtrsv_(char* uplo, char* trans, char* diag, int* n,
                double* A, int* lda, double* x, int* incx);
            
            void dgemv_(char* trans, int* m, int* n, double* alpha,
                double* A, int* lda, double* x, int* incx,
                double* beta, double* y, int* incy);
        }

        void backSolve(index_type n, real_matrix_type& A, real_vector_type& x) {
            char uplo = 'U';
            char trans = 'N';
            char diag = 'N';
            index_type lda = A.rows();
            index_type incx = 1;
            real_type* Aptr = A.data();
            real_type* xptr = x.data();
            dtrsv_(&uplo, &trans, &diag, &n, Aptr, &lda, xptr, &incx);
        }

        void matTimesVec(index_type n, real_matrix_type& A, real_vector_type& x, real_vector_type& result) {
            char trans = 'N';
            index_type m = A.rows();
            index_type lda = A.rows();
            index_type incx = 1;
            index_type incy = 1;
            real_type alpha = real_type(1);
            real_type beta = real_type(0);
            real_type* Aptr = A.data();
            real_type* xptr = x.data();
            real_type* yptr = result.data();
            dgemv_(&trans, &m, &n, &alpha, Aptr, &lda, xptr, &incx, &beta, yptr, &incy);
        }
    } // namespace details

    string ConvFlagToStr(ConvFlag flag) {
        switch(flag) {
        case ConvFlag::success:
            return "success -> convergence critera satisfied";
        case ConvFlag::unconverged:
            return "unconverged -> neither converged nor diverged";
        case ConvFlag::diverged:
            return "diverged -> excessive growth in residual norm";
        case ConvFlag::maxits:
            return "maximum iterations reached";
        case ConvFlag::true_rnrm:
            return "true residual norm failed convergence test";
        case ConvFlag::inf_or_nan:
            return "residual norm is either inf or nan";
        case ConvFlag::precon_fail:
            return "preconditioner application failed";
        case ConvFlag::matvec_fail:
            return "matrix-vector product failed";
        case ConvFlag::breakdown:
            return "input matrix or preconditioner are likely singular";
        case ConvFlag::stagnation:
            return "stagnation";
        }
        return string();
    }

    void checkGMRESParams(const GMRESParams& p) {
        if (p.kspaceSz < 1)
            throw invalid_argument("GMRESParams: kspaceSz < 1");
        if (p.maxits < 1)
            throw invalid_argument("GMRESParams: maxits < 1");
        if (p.relTol < 0)
            throw invalid_argument("GMRESParams: relTol < 0");
        if (p.absTol < 0)
            throw invalid_argument("GMRESParams: absTol < 0");
        if (p.divTol <= 0)
            throw invalid_argument("GMRESParams: divTol <= 0");
        if (p.stgTol < 0)
            throw invalid_argument("GMRESParams: stgTol < 0");
    }

    ostream& operator<<(ostream& strm, const GMRESOut& out) {
        strm << "outcome: " << out.flag << "\n";
        if (!out.msg.empty())
            strm << "info: " << out.msg.empty() << "\n";
        strm << "outer iter: " << out.outerIts << "\n";
        strm << "inner iter: " << out.innerIts << "\n";
        strm << "relative residual: " << std::scientific
            << std::setprecision(2) << out.relres << "\n";
        return strm;
    }
} // namespace blitzdg
