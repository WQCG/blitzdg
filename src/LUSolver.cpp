#include <suitesparse/umfpack.h>
#include <LUSolver.hpp>
#include <iostream>
#include <math.h>

using namespace std;

/**
 * Constructor. Takes a pointer reference to a blitz 2D array (The matrix A to be used by the solver in Ax=b).
 */
LUSolver::LUSolver(Array<double, 2> * const & Ain, SparseMatrixConverter const & _matrixConverter) {
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
    const Array<double, 2> & Aref = *A;
    const int n_rows = Aref.rows();
    const int n_cols = Aref.cols();

    const int nz = MatrixConverter.getNumNonZeros(*A);

    Ap = new int[n_rows+1];
    Ai = new int[nz];
    Ax = new double[nz];
    Map = new int[nz];

    cout << "Computing LU factorization!" << endl;

    // convert sparse Triplet to compressed column format
    MatrixConverter.fullToCompressedColumn(*A, Ap, Ai, Ax);

    umfpack_di_symbolic(n_rows, n_cols, Ap, Ai, Ax, &Symbolic, null, null);
    umfpack_di_numeric (Ap, Ai, Ax, Symbolic, &Numeric, null, null) ;
    cout << "Done!" << endl;
}

/**
 * Solve Ax=b using UMFPACK. Requires LUSolver.factorize() to be called first. 'x' is returned in 'soln' reference.
 * 'b' is specified by 'rhs' reference.
 */
void LUSolver::solve(Array<double, 1> const & rhs, Array<double,1> & soln) {
    int n = rhs.length(0);
    double * b = new double[n];

    for(int i=0; i<n; i++) {
        b[i] = rhs(i);
    }

    double * x = new double[n];

    cout << "Solving Ax = b.." << endl;
    umfpack_di_solve (UMFPACK_A, Ap, Ai, Ax, x, b, Numeric, null, null);
    cout << "Done." << endl;

    for(int i =0; i<n; i++) {
        soln(i) = x[i];
    }
}

/**
 * Returns a reference to the matrix A.
 */
Array<double, 2> & LUSolver::get_A() {
    return *A;
}

LUSolver::~LUSolver() {
    if (Triplet.row != nullptr) delete[] Triplet.row;
    if (Triplet.col != nullptr) delete[] Triplet.col;
    if (Triplet.val != nullptr) delete[] Triplet.val;
    if (Numeric != nullptr) umfpack_di_free_numeric (&Numeric);
    if (Symbolic != nullptr) umfpack_di_free_symbolic (&Symbolic) ;
}