// Copyright (C) 2017-2018  Derek Steinmoeller. 
// See COPYING and LICENSE files at project root for more details. 

#include <suitesparse/umfpack.h>
#include "LUSolver.hpp"
#include <iostream>
#include <cmath>

using std::cout;
using std::endl;

namespace blitzdg {
    /**
     * Constructor. Takes a pointer reference to a blitz 2D array (The matrix A to be used by the solver in Ax=b).
     */
    LUSolver::LUSolver(matrix_type* const & Ain, SparseMatrixConverter const & _matrixConverter) {
        A = Ain;
        MatrixConverter = _matrixConverter;
        Triplet.row = nullptr;
        Triplet.col = nullptr;
        Triplet.val = nullptr;
        Numeric = nullptr;
        Symbolic = nullptr;
    }

    /**
     * Factorize the matrix A with UMFPACK. Computes L,U factors and permutation matrices P,Q such that P*A*Q=LU.
     */
    void LUSolver::factorize() {
        const matrix_type& Aref = *A;
        const index_type n_rows = Aref.rows();
        const index_type n_cols = Aref.cols();
        const index_type nz = MatrixConverter.getNumNonZeros(*A);

        Ap = new index_type[n_rows+1];
        Ai = new index_type[nz];
        Ax = new real_type[nz];
        Map = new index_type[nz];

        cout << "Computing LU factorization!" << endl;

        // convert sparse Triplet to compressed column format
        MatrixConverter.fullToCompressedColumn(*A, Ap, Ai, Ax);

        umfpack_di_symbolic(n_rows, n_cols, Ap, Ai, Ax, &Symbolic, (double *) NULL, (double *) NULL);
        umfpack_di_numeric (Ap, Ai, Ax, Symbolic, &Numeric, (double *) NULL, (double *) NULL) ;
        cout << "Done!" << endl;
    }

    /**
     * Solve Ax=b using UMFPACK. Requires LUSolver.factorize() to be called first. 'x' is returned in 'soln' reference.
     * 'b' is specified by 'rhs' reference.
     */
    void LUSolver::solve(vector_type const & rhs, vector_type& soln) {
        index_type n = rhs.length(0);
        real_type* b = new real_type[n];

        for(index_type i=0; i<n; i++) {
            b[i] = rhs(i);
        }

        real_type* x = new real_type[n];

        cout << "Solving Ax = b.." << endl;
        umfpack_di_solve (UMFPACK_A, Ap, Ai, Ax, x, b, Numeric, (double *)NULL, (double  *)NULL);
        cout << "Done." << endl;

        for(index_type i =0; i<n; i++) {
            soln(i) = x[i];
        }
    }

    /**
     * Returns a reference to the matrix A.
     */
    matrix_type& LUSolver::get_A() {
        return *A;
    }

    LUSolver::~LUSolver() {
        if (Triplet.row != nullptr) delete[] Triplet.row;
        if (Triplet.col != nullptr) delete[] Triplet.col;
        if (Triplet.val != nullptr) delete[] Triplet.val;
        if (Numeric != nullptr) umfpack_di_free_numeric (&Numeric);
        if (Symbolic != nullptr) umfpack_di_free_symbolic (&Symbolic);
    }
} // namespace blitzdg