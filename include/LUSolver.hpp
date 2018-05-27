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
     * Default constructor.
     */
    LUSolver()
      : order_{ 0 }, mat_{ nullptr }, 
      symbolic_{ nullptr }, numeric_{ nullptr }
    {}

    /**
     * Computes an LU factorization of the input CSC matrix.
     * Overwrites any existing factorization.
     * @param[in] mat The CSC matrix.
     */
    void factorize(const CSCMat& mat);

    /**
     * Solves the linear system A*soln = rhs.
     * @param[in] rhs The right-hand side.
     * @param[out] soln The solution of the linear system.
     */
    void solve(const real_vector_type& rhs, real_vector_type& soln) const;

    /**
     * Destructor.
     */
    ~LUSolver() {
      freeMem();
    }
  private:
    index_type order_; // order of the current linear system
    const CSCMat* mat_; // non-owning pointer to input matrix
    void* symbolic_; // pointer to symbolic factorization structure
    void* numeric_; // pointer to numeric factorization structure

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

