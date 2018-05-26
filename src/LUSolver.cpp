// Copyright (C) 2017-2018  Waterloo Quantitative Consulting Group, Inc.
// See COPYING and LICENSE files at project root for more details.

#include "LUSolver.hpp"
#include <suitesparse/umfpack.h>
#include <stdexcept>

using std::runtime_error;

namespace blitzdg {
    bool LUSolver::symbolicFactorize() {
        index_type flag = umfpack_di_symbolic(mat_->rows(), mat_->cols(), mat_->colPtrs(), 
            mat_->rowInds(), mat_->elems(), &symbolic_, (double*)NULL, (double*)NULL);
        return flag == UMFPACK_OK;
    }

    bool LUSolver::numericFactorize() {
        index_type flag = umfpack_di_numeric(mat_->colPtrs(), mat_->rowInds(), 
            mat_->elems(), symbolic_, &numeric_, (double*)NULL, (double*)NULL);
        return flag == UMFPACK_OK;
    }

    void LUSolver::freeMem() {
        if (symbolic_) {
            umfpack_di_free_symbolic(&symbolic_);
            symbolic_ = nullptr;
        }
        if (numeric_) {
            umfpack_di_free_numeric(&numeric_); 
            numeric_ = nullptr;
        }
    }

    void LUSolver::factorize(const CSCMat& mat) {
        if (mat.rows() != mat.cols())
            throw runtime_error("LUSolver::factorize: input matrix is not square");
        order_ = mat.rows();
        mat_ = &mat;
        freeMem(); // free any existing memory
        if (!symbolicFactorize())
            throw runtime_error("LUSolver::factorize: symbolic factorization failed");
        else if (!numericFactorize())
            throw runtime_error("LUSolver::factorize: numeric factorization failed");
        // free the symbolic factorization since we don't need it any longer
        umfpack_di_free_symbolic(&symbolic_);
    }

    void LUSolver::solve(const vector_type& rhs, vector_type& soln) const {
        if (!numeric_ || !mat_) // check that factorize has been called
            throw runtime_error("LUSolver::solve: call factorize before calling solve");
        if (rhs.length(0) < order_) // check length of rhs
            throw runtime_error("LUSolver::solve: rhs length is less than matrix order");
        if (soln.length(0) < order_) // check length of soln
            throw runtime_error("LUSolver::solve: soln length is less than matrix order");
        index_type flag = umfpack_di_solve (UMFPACK_A, mat_->colPtrs(), mat_->rowInds(), 
            mat_->elems(), soln.data(), rhs.data(), numeric_, (double*)NULL, (double*)NULL);
        if (flag != UMFPACK_OK)
            throw runtime_error("LUSolver::solve: failed");
    }
} // namespace blitzdg