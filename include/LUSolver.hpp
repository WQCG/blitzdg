// Copyright (C) 2017-2018  Derek Steinmoeller. 
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
      index_type N;
      matrix_type* A;

      // Umfpack-specific fields
      index_type * Ap;
      index_type * Ai;
      real_type * Ax;
      index_type * Map;
      real_type * null;

      void * Symbolic;
      void * Numeric;

      SparseTriplet Triplet;

      SparseMatrixConverter MatrixConverter;
    
    public:
      LUSolver(matrix_type* const &, SparseMatrixConverter const &);
      
      matrix_type& get_A();

      void factorize();

      void solve(vector_type const &, vector_type&);

      ~LUSolver();
  };
} // namespace blitzdg

