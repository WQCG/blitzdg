// Copyright (C) 2017-2018  Waterloo Quantitative Consulting Group, Inc.
// See COPYING and LICENSE files at project root for more details.

/**
 * @file SparseMatrixConverter.hpp
 * @brief Defines the SparseMatrixConverter class that converts between blitz++
 * 2D arrays, sparse triplet matrices, compressed sparse column (CSC) matrices, and contiguous 
 * C-style (POD) arrays.
 * @note When converting from a dense matrix to a sparse representation, an element \f$a_{ij}\f$ of
 * the dense matrix is considered nonzero if \f$|a_{ij}| \geq 2\varepsilon\f$, where \f$\varepsilon\f$
 * is machine precision for the type real_type. This approach is also employed when counting the
 * the number of nonzero elements in a dense matrix.
 */
#pragma once
#include "SparseTriplet.hpp"
#include "Types.hpp"

namespace blitzdg {
  class SparseMatrixConverter {
  public:
      /**
       * Converts a dense matrix to a sparse triplet matrix.
       * @param[in] A Dense matrix.
       * @param[out] triplet Sparse triplet matrix.
       */
      void fullToSparseTriplet(const matrix_type& A, SparseTriplet & triplet) const;

      /**
       * Converts a sparse triplet matrix to a compressed sparse column matrix with numRows rows and numCols columns.
       * @param[in] numRows Number of rows in CSC matrix.
       * @param[in] numCols Number of columns in CSC matrix.
       * @param[in] triplet Sparse triplet matrix.
       * @param[out] Aptr Pointer to array of the locations in the Avalues array that start a column.
       * @param[out] Aind Pointer to array of row indices of elements in Avalues.
       * @param[out] Avalues Pointer to array of nonzero elements.
       */
      void sparseTripletToCompressedColumn(index_type numRows, index_type numCols, const SparseTriplet & triplet, index_type * Aptr, index_type * Aind, real_type * Avalues) const;

      /**
       * Converts a dense matrix to a compressed sparse column matrix.
       * @param[in] A Dense matrix.
       * @param[out] Aptr Pointer to array of the locations in the Avalues array that start a column.
       * @param[out] Aind Pointer to array of row indices of elements in Avalues.
       * @param[out] Avalues Pointer to array of nonzero elements.
       */
      void fullToCompressedColumn(const matrix_type& A, index_type * Aptr, index_type * Aind, real_type * Avalues) const;

      /**
       * Vectorizes a dense matrix. The rows of A are stored contiguously and in ascending order of row index.
       * @param[in] A Dense matrix.
       * @param[out] Apod Pointer to array that stores the vectorized matrix.
       * @note The array pointed to by Apod must have space allocated for at 
       * least m*n elements, where m is the the number of rows of A and n is 
       * the number of columns.
       */
      void fullToPodArray(const matrix_type& A, real_type * Apod) const;

      /**
       * Reshapes an array into a dense matrix. The dense matrix is filled 
       * rowwise in ascending order of row index.
       * @param[in] Apod Pointer to an array.
       * @param[out] A Dense matrix.
       * @note The array pointed to by Apod must have at least m*n elements, 
       * where m is the the number of rows of A and n is the number of columns.
       */
      void podArrayToFull(const real_type * Apod, matrix_type& A) const;

      /**
       * Returns the number of nonzero elements in a dense matrix.
       * @param[in] A Dense matrix.
       * @return The number of nonzeros.
       */
      index_type getNumNonZeros(const matrix_type& A) const;
  };
} // namespace blitzdg
