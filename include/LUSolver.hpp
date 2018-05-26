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
  /**
   * Implements a class for solving linear systems via LU factorization
   * in which the coefficient matrix is stored as a compressed sparse
   * column (CSC) matrix.
   */
  class LUSolver {   
  public:
    /**
     * Constructor that builds an LUSolver from a CSC matrix.
     * @param[in] mat The CSC matrix to be factorized.
     */
    explicit LUSolver(const CSCMat& mat)
      : mat_ { mat }, symbolic_{ nullptr }, numeric_{ nullptr }
    {}

    /**
     * Returns a reference to the underlying CSC matrix.
     */
    const CSCMat& getMatrix() const {
      return mat_;
    }

    /**
     * Computes an LU factorization of the CSC matrix.
     */
    void factorize();

    /**
     * Solves the linear system A*soln = rhs.
     * @param[in] rhs The right-hand side.
     * @param[out] soln The solution of the linear system.
     */
    void solve(const vector_type& rhs, vector_type& soln) const;

    /**
     * Destructor.
     */
    ~LUSolver() {
      freeMem();
    }
  private:
    CSCMat mat_;
    void* symbolic_;
    void* numeric_;

    /**
     * Computes the symbolic factorization and
     * returns true if successful. 
     */
    bool symbolicFactorize();

    /**
     * Computes the numeric factorization and
     * returns true if successful.
     */
    bool numericFactorize();

    /**
     * Frees any memory allocated for either the symbolic
     * or numeric factorization.
     */
    void freeMem();
  };
} // namespace blitzdg

