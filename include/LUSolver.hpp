// Copyright (C) 2017-2018  Waterloo Quantitative Consulting Group, Inc. 
// See COPYING and LICENSE files at project root for more details. 

/**
 * @file LUSolver.hpp
 * @brief Defines the LUSolver class that implements UMFPACK LU factorization.
 */

#pragma once
#include "SparseMatrixConverter.hpp"
#include "Types.hpp"

namespace blitzdg {
  class LUSolver {
      const matrix_type* A;

      // Umfpack-specific fields
      index_type * Ap;
      index_type * Ai;
      real_type * Ax;
      index_type * Map;

      void * Symbolic;
      void * Numeric;
      
      SparseMatrixConverter MatrixConverter;
    
    public:
      explicit LUSolver(const matrix_type* Ain);
      
      const matrix_type& get_A() const;

      void factorize();

      void solve(vector_type const &, vector_type&);

      ~LUSolver();
  };
} // namespace blitzdg

