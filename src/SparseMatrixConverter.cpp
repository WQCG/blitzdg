#include <SparseMatrixConverter.hpp>

using namespace std;

SparseMatrixConverter::SparseMatrixConverter() {

}

void SparseMatrixConverter::fullToSparseTriplet(const Array<double, 2> & A, SparseTriplet & triplet) {

    const double eps = numeric_limits<double>::epsilon();
    const int n_rows = A.rows();
    const int n_cols = A.cols();

    triplet.nz = n_rows*n_cols;  // find a better bound, or don't.
    triplet.col = new int[triplet.nz];
    triplet.row = new int[triplet.nz];
    triplet.val = new double[triplet.nz];

    int nz = 0;

    for( int i=0; i < n_rows; i++ ) {
        for ( int j=0; j < n_cols; j++ ) {
            double val = A(i,j);
            if ( abs(val) < 2*eps ) {
                continue;
            }

            triplet.row[nz] = i;
            triplet.col[nz] = j;
            triplet.val[nz] = val;

            nz++;
        }
    }

    triplet.nz = nz;

}

SparseMatrixConverter::~SparseMatrixConverter() {

}