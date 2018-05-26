// Copyright (C) 2017-2018  Waterloo Quantitative Consulting Group, Inc. 
// See COPYING and LICENSE files at project root for more details. 

/**
 * @file LUSolver.hpp
 * @brief Defines the LUSolver class that implements UMFPACK LU factorization
 * (umfpack_di_numeric, umfpack_di_solve) for sparse matrices stored in compressed 
 * sparse column (CSC) format. UMFPACK is part of the SuiteSparse package: 
 * http://faculty.cse.tamu.edu/davis/suitesparse.html.
 */

#pragma once
#include "CSCMatrix.hpp"
#include "Types.hpp"

namespace blitzdg {
  class LUSolver {   
  public:
    explicit LUSolver(const CSCMat& mat)
      : mat_{ mat }, symbolic_{ nullptr }, numeric_{ nullptr }
    {}

    const CSCMat& getMatrix() const {
      return mat_;
    }

    void factorize();

    void solve(const vector_type& rhs, vector_type& soln) const;

    ~LUSolver();
  private:
    const CSCMat& mat_;
    void* symbolic_;
    void* numeric_;
  };
} // namespace blitzdg

