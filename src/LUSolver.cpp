#include <suitesparse/umfpack.h>
#include <LUSolver.hpp>
#include <iostream>
#include <math.h>

using namespace std;

LUSolver::LUSolver(Array<double, 2> * const & Ain) {
    A = Ain;
    Triplet.row = nullptr;
    Triplet.col = nullptr;
    Triplet.val = nullptr;
    Numeric = nullptr;
    Symbolic = nullptr;
}

void LUSolver::factorize() {
    const Array<double, 2> & Aref = *A;
    const int n_rows = Aref.rows();
    const int n_cols = Aref.cols();

    // Convert full matrix A to sparse Triplet.
    toSparseTriplet();

    const int nz = Triplet.nz;

    Ap = new int[n_rows+1];
    Ai = new int[nz];
    Ax = new double[nz];
    Map = new int[nz];

    cout << "Computing LU factorization!" << endl;
    // convert sparse Triplet to compressed column format
    umfpack_di_triplet_to_col(n_rows, n_cols, Triplet.nz, Triplet.row, Triplet.col,
        Triplet.val, Ap, Ai, Ax, Map);

    umfpack_di_symbolic(n_rows, n_cols, Ap, Ai, Ax, &Symbolic, null, null);
    umfpack_di_numeric (Ap, Ai, Ax, Symbolic, &Numeric, null, null) ;
    cout << "Done!" << endl;
}

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

void LUSolver::toSparseTriplet() {
    const Array<double, 2> & Aref = *A;

    const int n_rows = Aref.rows();
    const int n_cols = Aref.cols();

    Triplet.nz = n_rows*n_cols;  // find a better bound, or don't.
    Triplet.col = new int[Triplet.nz];
    Triplet.row = new int[Triplet.nz];
    Triplet.val = new double[Triplet.nz];

    int nz = 0;

    for( int i=0; i < n_rows; i++ ) {
        for ( int j=0; j < n_cols; j++ ) {
            double val = Aref(i,j);
            if ( abs(val) < 1.e-15 ) {
                continue;
            }

            Triplet.row[nz] = i;
            Triplet.col[nz] = j;
            Triplet.val[nz] = val;

            nz++;
        }
    }

    Triplet.nz = nz;
}

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