// Copyright (C) 2017-2018  Waterloo Quantitative Consulting Group, Inc.
// See COPYING and LICENSE files at project root for more details.

#include "LUSolver.hpp"
#include <suitesparse/umfpack.h>
#include <stdexcept>

using std::runtime_error;

namespace blitzdg {
    void LUSolver::factorize() {
        index_type flag = umfpack_di_symbolic(mat_.rows(), mat_.cols(), mat_.colPtrs(), 
            mat_.rowInds(), mat_.elems(), &symbolic_, (double*)NULL, (double*)NULL);

        if (flag != UMFPACK_OK)
             throw runtime_error("LUSolver::factorize: symbolic factorization failed");
        
        flag = umfpack_di_numeric (mat_.colPtrs(), mat_.rowInds(), mat_.elems(), 
            symbolic_, &numeric_, (double*)NULL, (double*)NULL);

        if (flag != UMFPACK_OK)
            throw runtime_error("LUSolver::factorize: numeric factorization failed");
        
        // delete the symbolic factorization since we don't need it any longer
        if (symbolic_) {
            umfpack_di_free_symbolic(&symbolic_);
            symbolic_ = nullptr;
        }
    }

    void LUSolver::solve(vector_type const & rhs, vector_type& soln) const {
        index_type flag = umfpack_di_solve (UMFPACK_A, mat_.colPtrs(), mat_.rowInds(), 
            mat_.elems(), soln.data(), rhs.data(), numeric_, (double*)NULL, (double*)NULL);
        if (flag != UMFPACK_OK)
            throw runtime_error("LUSolver::solve: failed");
    }

    LUSolver::~LUSolver() {
        if (numeric_) {
            umfpack_di_free_numeric (&numeric_); 
            numeric_ = nullptr;
        }
    }
} // namespace blitzdg