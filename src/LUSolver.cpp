#include <suitesparse/umfpack.h>
#include <LUSolver.hpp>
#include <iostream>
#include <math.h>

using namespace std;

LUSolver::LUSolver(Array<double, 2> * const & Ain) {
    A = Ain;
    triplet.row = nullptr;
    triplet.col = nullptr;
    triplet.val = nullptr;
}

void LUSolver::factorize() {
    const Array<double, 2> & Aref = *A;
    const int n_rows = Aref.rows();
    const int n_cols = Aref.cols();

    // Convert full matrix A to sparse triplet.
    toSparseTriplet();

    const int nz = triplet.nz;

    int *Ap = new int[n_rows+1];
    int *Ai = new int[nz];
    double *Ax = new double[nz];
    int *map = new int[nz];

    void *Symbolic, *Numeric;
    double *null = (double *) NULL ;

    double b[] = {8., 45., -3., 3., 19.};
    double x[5];

    cout << "Computing LU factorization!" << endl;
    // convert sparse triplet to compressed column format
    umfpack_di_triplet_to_col(n_rows, n_cols, triplet.nz, triplet.row, triplet.col,
        triplet.val, Ap, Ai, Ax, map);

    umfpack_di_symbolic(n_rows, n_cols, Ap, Ai, Ax, &Symbolic, null, null);
    umfpack_di_numeric (Ap, Ai, Ax, Symbolic, &Numeric, null, null) ;
    cout << "Done!" << endl;

    umfpack_di_free_symbolic (&Symbolic) ;

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

    triplet.nz = n_rows*n_cols;  // find a better bound, or don't.
    triplet.col = new int[triplet.nz];
    triplet.row = new int[triplet.nz];
    triplet.val = new double[triplet.nz];

    int nz = 0;

    for( int i=0; i < n_rows; i++ ) {
        for ( int j=0; j < n_cols; j++ ) {
            if ( abs(Aref(i,j)) < 1.e-15 ) {
                continue;
            }

            triplet.row[nz] = i;
            triplet.col[nz] = j;
            triplet.val[nz] = Aref(i,j);

            nz++;
        }
    }

    triplet.nz = nz;
}

Array<double, 2> & LUSolver::get_A() {
    return *A;
}


LUSolver::~LUSolver() {
    if (triplet.row != nullptr) delete[] triplet.row;
    if (triplet.col != nullptr) delete[] triplet.col;
    if (triplet.val != nullptr) delete[] triplet.val;
}