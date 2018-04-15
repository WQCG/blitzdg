// Copyright (C) 2017-2018  Derek Steinmoeller. 
// See COPYING and LICENSE files at project root for more details. 

#include <SparseMatrixConverter.hpp>
#include <suitesparse/umfpack.h>

using namespace std;

SparseMatrixConverter::SparseMatrixConverter() {

}

void SparseMatrixConverter::podArrayToFull(const double * Apod, Array<double, 2> & A) {
    int ind = 0;
    for (int i=0; i<A.rows(); i++) {
        for (int j=0; j<A.cols(); j++) {
            A(i,j) = Apod[ind];
            ind++;
        }
    }
}

void SparseMatrixConverter::fullToPodArray(const Array<double, 2> & A, double * Apod) {
    int ind = 0;
    for (int i =0; i< A.rows(); i++) {
        for (int j=0; j< A.cols(); j++) {
            Apod[ind] = A(i,j);
            ind++;
        }
    }
}

void SparseMatrixConverter::fullToCompressedColumn(const Array<double, 2> & A,
                                                int * Aptr, int * Aind, double * Avalues) {
    SparseTriplet triplet;

    fullToSparseTriplet(A, triplet);

    sparseTripletToCompressedColumn(A.rows(), A.cols(), triplet, Aptr, Aind, Avalues);
}

void SparseMatrixConverter::sparseTripletToCompressedColumn(const int numRows, const int numCols, const SparseTriplet & triplet,
                                                        int * Aptr, int * Aind, double * Avalues) {

    const int nz = triplet.nz;
    int * map = new int[nz];

    umfpack_di_triplet_to_col(numRows, numCols, nz, triplet.row, triplet.col,
        triplet.val, Aptr, Aind, Avalues, map);
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

int SparseMatrixConverter::getNumNonZeros(const Array<double, 2> & A) {
    const double eps = numeric_limits<double>::epsilon();
    const int n_rows = A.rows();
    const int n_cols = A.cols();

    int nz = 0;

    for( int i=0; i < n_rows; i++ ) {
        for ( int j=0; j < n_cols; j++ ) {
            double val = A(i,j);
            if ( abs(val) < 2*eps ) {
                continue;
            }

            nz++;
        }
    }

    return nz;
}

SparseMatrixConverter::~SparseMatrixConverter() {

}