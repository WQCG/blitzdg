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
}

void LUSolver::factorize() {
    const Array<double, 2> & Aref = *A;
    const int n_rows = Aref.rows();
    const int n_cols = Aref.cols();

    // Convert full matrix A to sparse Triplet.
    toSparseTriplet();

    const int nz = Triplet.nz;

    int *Ap = new int[n_rows+1];
    int *Ai = new int[nz];
    double *Ax = new double[nz];
    int *map = new int[nz];

    void *symbolic;
    double *null = (double *) NULL ;

    double b[] = {8., 45., -3., 3., 19.};
    double x[5];

    cout << "Computing LU factorization!" << endl;
    // convert sparse Triplet to compressed column format
    umfpack_di_triplet_to_col(n_rows, n_cols, Triplet.nz, Triplet.row, Triplet.col,
        Triplet.val, Ap, Ai, Ax, map);

    umfpack_di_symbolic(n_rows, n_cols, Ap, Ai, Ax, &symbolic, null, null);
    umfpack_di_numeric (Ap, Ai, Ax, symbolic, &Numeric, null, null) ;
    cout << "Done!" << endl;

    umfpack_di_free_symbolic (&symbolic) ;

    umfpack_di_solve (UMFPACK_A, Ap, Ai, Ax, x, b, Numeric, null, null) ;
    umfpack_di_free_numeric (&Numeric) ;

    Array<double, 2> xa(5,1);
    xa(0) = x[0];
    xa(1) = x[1];
    xa(2) = x[2];
    xa(3) = x[3];
    xa(4) = x[4];

    cout << "x: " << xa << endl;
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
            if ( abs(Aref(i,j)) < 1.e-15 ) {
                continue;
            }

            Triplet.row[nz] = i;
            Triplet.col[nz] = j;
            Triplet.val[nz] = Aref(i,j);

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
}