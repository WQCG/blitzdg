// Copyright (C) 2017-2018  Waterloo Quantitative Consulting Group, Inc.
// See COPYING and LICENSE files at project root for more details.

#include "SparseMatrixConverter.hpp"
#include <suitesparse/umfpack.h>
#include <limits>

using std::numeric_limits;

namespace blitzdg {
    SparseMatrixConverter::SparseMatrixConverter() {

    }

    void SparseMatrixConverter::podArrayToFull(const real_type * Apod, matrix_type & A) {
        index_type ind = 0;
        for (index_type i=0; i<A.rows(); i++) {
            for (index_type j=0; j<A.cols(); j++) {
                A(i,j) = Apod[ind];
                ind++;
            }
        }
    }

    void SparseMatrixConverter::fullToPodArray(const matrix_type & A, real_type * Apod) {
        index_type ind = 0;
        for (index_type i =0; i< A.rows(); i++) {
            for (index_type j=0; j< A.cols(); j++) {
                Apod[ind] = A(i,j);
                ind++;
            }
        }
    }

    void SparseMatrixConverter::fullToCompressedColumn(const matrix_type & A,
                                                    index_type * Aptr, index_type * Aind, real_type * Avalues) {
        SparseTriplet triplet;

        fullToSparseTriplet(A, triplet);

        sparseTripletToCompressedColumn(A.rows(), A.cols(), triplet, Aptr, Aind, Avalues);
    }

    void SparseMatrixConverter::sparseTripletToCompressedColumn(const index_type numRows, const index_type numCols, const SparseTriplet & triplet,
                                                            index_type * Aptr, index_type * Aind, real_type * Avalues) {

        const index_type nz = triplet.nz;
        index_type * map = new index_type[nz];

        umfpack_di_triplet_to_col(numRows, numCols, nz, triplet.row, triplet.col,
            triplet.val, Aptr, Aind, Avalues, map);
    }

    void SparseMatrixConverter::fullToSparseTriplet(const matrix_type & A, SparseTriplet & triplet) {

        const real_type eps = numeric_limits<real_type>::epsilon();
        const index_type n_rows = A.rows();
        const index_type n_cols = A.cols();

        triplet.nz = n_rows*n_cols;  // find a better bound, or don't.
        triplet.col = new index_type[triplet.nz];
        triplet.row = new index_type[triplet.nz];
        triplet.val = new real_type[triplet.nz];

        index_type nz = 0;

        for( index_type i=0; i < n_rows; i++ ) {
            for ( index_type j=0; j < n_cols; j++ ) {
                real_type val = A(i,j);
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

    index_type SparseMatrixConverter::getNumNonZeros(const matrix_type & A) {
        const real_type eps = numeric_limits<real_type>::epsilon();
        const index_type n_rows = A.rows();
        const index_type n_cols = A.cols();

        index_type nz = 0;

        for( index_type i=0; i < n_rows; i++ ) {
            for ( index_type j=0; j < n_cols; j++ ) {
                real_type val = A(i,j);
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
} // namespace blitzdg