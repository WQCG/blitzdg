// Copyright (C) 2017-2018  Derek Steinmoeller. 
// See COPYING and LICENSE files at project root for more details. 

/**
 * @file LUSolver.hpp
 * @brief Defines the LUSolver class that implements UMFPACK LU factorization.
 */

#pragma once
#include <blitz/array.h>
#include "SparseMatrixConverter.hpp"

using namespace blitz;

namespace blitzdg {
  class LUSolver {
      int N;
      Array<double, 2> * A;

      // Umfpack-specific fields
      int * Ap;
      int * Ai;
      double * Ax;
      int * Map;
      double * null;

      void * Symbolic;
      void * Numeric;

      SparseTriplet Triplet;

      SparseMatrixConverter MatrixConverter;
    
    public:
      LUSolver(Array<double, 2> * const &, SparseMatrixConverter const &);
      
      Array<double, 2> & get_A();

      void factorize();

      void solve(Array<double, 1> const &, Array<double, 1> &);

      ~LUSolver();
  };
} // namespace blitzdg

