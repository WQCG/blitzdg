#include "DenseMatrixHelpers.hpp"
#include <cmath>
#include <limits>
#include <stdexcept>

using std::abs;
using std::numeric_limits;
using std::runtime_error;

namespace blitzdg {
    index_type countNonzeros(const matrix_type& mat, real_type dropTol) {
        real_type nnz = 0;
        for (matrix_type::const_iterator itr = mat.begin(); itr != mat.end(); ++itr) {
            if (abs(*itr) > dropTol)
                ++nnz;
        }
        if (nnz > numeric_limits<index_type>::max())
            throw runtime_error("countNonzeros: number of nonzero elements exceeds maximum allowable");
        return static_cast<index_type>(nnz);
    }

    void fullToPodArray(const matrix_type& mat, real_type* arr, bool byRows) {
        index_type nnz = 0;
        if (byRows) {
            for (index_type i = 0; i < mat.rows(); ++i) {
                for (index_type j = 0; j < mat.cols(); ++j)
                    arr[nnz++] = mat(i, j);
            }
        }
        else {
            for (index_type j = 0; j < mat.cols(); ++j) {
                for (index_type i = 0; i < mat.rows(); ++i)
                    arr[nnz++] = mat(i, j);
            }
        }
    }

    void podArrayToFull(const real_type* arr, matrix_type& mat, bool byRows) {
        index_type nnz = 0;
        if (byRows) {
            for (index_type i = 0; i < mat.rows(); ++i) {
                for (index_type j = 0; j < mat.cols(); ++j)
                    mat(i, j) = arr[nnz++];
            }
        }
        else {
            for (index_type j = 0; j < mat.cols(); ++j) {
                for (index_type i = 0; i < mat.rows(); ++i)
                   mat(i, j) = arr[nnz++];
            }
        }
    }
} // namespace blitzdg