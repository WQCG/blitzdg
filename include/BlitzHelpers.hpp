// Copyright (C) 2017-2022  Waterloo Quantitative Consulting Group, Inc.
// See COPYING and LICENSE files at project root for more details.

/**
 * @file BlitzHelpers.hpp
 * @brief Defines a set of utility functions for working with blitz arrays
 * and their matrix and vector aliases.
 */
#pragma once
#include "Traits.hpp"
#include "Types.hpp"
#include <blitz/array.h>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <algorithm>
#include <cmath>
#include <numeric>
#include <utility>
#include <vector>


namespace blitzdg {
    namespace details {
        // Exact lexicographic comparison for a given ordering.
        template <typename T>
        class CompareExact {
        public:
            explicit CompareExact(const matrix_type<T>& A)
                : ptr_{ &A }, ord_(A.cols())
            {
                std::iota(ord_.begin(), ord_.end(), index_type(0));
            }

            CompareExact(const matrix_type<T>& A, const std::vector<index_type>& order)
                : ptr_{ &A }, ord_{ order }
            {}

            bool operator()(index_type lhs, index_type rhs) const;
        private:
            const matrix_type<T>* ptr_;
            std::vector<index_type> ord_;
        };

        template <typename T>
        bool CompareExact<T>::operator()(index_type lhs, index_type rhs) const {
            for (auto j : ord_) {
                if ((*ptr_)(lhs, j) < (*ptr_)(rhs, j)) {
                    return true;
                }
                else if ((*ptr_)(lhs, j) > (*ptr_)(rhs, j)) {
                    return false;
                }
            }
            return false;
        }

        // Equality comparison with tolerance.
        template <typename T>
        class CompareEQ {
        public:
            explicit CompareEQ(const matrix_type<T>& A, real_type tol = real_type(0));

            bool operator()(index_type lhs, index_type rhs) const;
        private:
            const matrix_type<T>* ptr_;
            std::vector<real_type> tols_;
        };

        template <typename T>
        CompareEQ<T>::CompareEQ(const matrix_type<T>& A, real_type tol)
            : ptr_{ &A }, tols_(A.cols(), real_type(0))
        {
            if (tol > real_type(0)) {
                for (index_type i = 0; i < A.rows(); ++i) {
                    for (index_type j = 0; j < A.cols(); ++j) {
                        tols_[j] = std::max(tols_[j], std::abs(A(i, j)));
                    }
                }
                for (auto& t : tols_) {
                    t *= tol;
                }
            }
        }

        template <typename T>
        bool CompareEQ<T>::operator()(index_type lhs, index_type rhs) const {
            for (index_type j = 0; j < ptr_->cols(); ++j) {
                if (std::abs((*ptr_)(lhs, j) - (*ptr_)(rhs, j)) > tols_[j]) {
                    return false;
                }
            }
            return true;
        }

        // Equality comparison with tolerance over a single dimension.
        template <typename T>
        class CompareEQByDim {
        public:
            CompareEQByDim(const matrix_type<T>& A, index_type dim, real_type tol = real_type(0));

            bool operator()(index_type lhs, index_type rhs) const {
                return std::abs((*ptr_)(lhs, dim_) - (*ptr_)(rhs, dim_)) <= tol_;
            }
        private:
            const matrix_type<T>* ptr_;
            index_type dim_;
            real_type tol_;
        };

        template <typename T>
        CompareEQByDim<T>::CompareEQByDim(const matrix_type<T>& A, index_type dim, real_type tol)
            : ptr_{ &A }, dim_{ dim }, tol_{ real_type(0) }
        {
            if (tol > real_type(0)) {
                for (index_type i = 0; i < A.rows(); ++i) {
                    tol_ = std::max(tol_, std::abs(A(i, dim)));
                }
                tol_ *= tol;
            }
        }

        // Order the columns of A by decreasing size of their variance.
        template <typename T>
        std::vector<index_type> getOrdering(const matrix_type<T>& A) {
            if (A.cols() < 2) {
                return std::vector<index_type>{ index_type(0) };
            }
            std::vector<real_type> mn(A.cols(), real_type(0)), var(A.cols(), real_type(0));
            for (index_type i = 0; i < A.rows(); ++i) {
                for (index_type j = 0; j < A.cols(); ++j) {
                    auto delta = A(i, j) - mn[j];
                    mn[j] += delta / (i + 1);
                    auto delta2 = A(i, j) - mn[j];
                    var[j] += delta * delta2;
                }
            }
            std::vector<index_type> order(A.cols());
            std::iota(order.begin(), order.end(), index_type(0));
            std::sort(order.begin(), order.end(),
                [&var](index_type i, index_type j) { return var[i] > var[j]; });
            return order;
        }
    } // namespace details

    /**
     * Used for constructing a matrix with column-major ordering.
     * 
     * The following example shows how to build a 5 x 5
     * real matrix with column-major ordering.
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     * real_matrix_type mat(5, 5, ColumnMajorOrder());
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     */
    using ColumnMajorOrder = blitz::ColumnMajorArray<2>;

    /**
     * Used for constructing a matrix with row-major ordering.
     * 
     * The following example shows how to build a 5 x 5
     * real matrix with row-major ordering.
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     * real_matrix_type mat(5, 5, RowMajorOrder());
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     * @note: Row-major ordering is the default ordering
     * and does not need to be explicitly specified.
     */
    using RowMajorOrder = blitz::GeneralArrayStorage<2>;

    /**
     * Returns true if the matrix use row-major storage.
     */
    template <typename T>
    bool isRowMajor(const matrix_type<T>& mat) {
        return mat.isMajorRank(0);
    }

    /**
     * Returns true if the matrix uses column-major storage.
     */
    template <typename T>
    bool isColumnMajor(const matrix_type<T>& mat) {
        return mat.isMajorRank(1);
    }

    /**
     * Returns the number of nonzero elements in the dense matrix.
     * @param[in] mat The dense matrix.
     * @param[in] dropTol The drop tolerance. Defaults to 0.
     * @note An element \f$a_{ij}\f$ of the input matrix is considered
     * nonzero if \f$|a_{ij}| > \mathrm{dropTol}\f$. 
     */
    template <typename T>
    index_type countNonzeros(const matrix_type<T>& mat, real_type dropTol = real_type(0)) {
        std::size_t nnz = 0;
        for (typename matrix_type<T>::const_iterator itr = mat.begin(); itr != mat.end(); ++itr) {
            if (std::abs(*itr) > dropTol)
                ++nnz;
        }
        if (nnz > std::numeric_limits<index_type>::max())
            throw std::runtime_error("countNonzeros: number of nonzero elements exceeds maximum allowable");
        return static_cast<index_type>(nnz);
    }

    /**
     * Vectorizes a dense matrix by writing it to an array.
     * @param[in] mat The Dense matrix.
     * @param[in] byRows If true, the matrix is copied to the array by rows. Defaults to true.
     * @param[out] arrItr An output iterator to an array.
     * @note The array must have space allocated for at 
     * least m*n elements, where m is the the number of rows of mat and n is 
     * the number of columns. 
     */
    template <typename T, typename OutputItr>
    void reshapeMatTo1D(const matrix_type<T>& mat, OutputItr arrItr, bool byRows = true) {
        // Fail at compile time if the type T is not
        // the same as the value type of OutputItr.
        static_assert(iteratorValueTypeSameAs<OutputItr, T>(),
        "Matrix value type differs from array value type");
        if (byRows) {
            for (index_type i = 0; i < mat.rows(); ++i) {
                for (index_type j = 0; j < mat.cols(); ++j)
                    *arrItr++ = mat(i, j);
            }
        }
        else {
            for (index_type j = 0; j < mat.cols(); ++j) {
                for (index_type i = 0; i < mat.rows(); ++i)
                    *arrItr++ = mat(i, j);
            }
        }
    }

    /**
     * Reshapes an array to a dense matrix.
     * @param[in] arrItr An input iterator to the array.
     * @param[in] byRows If true, the array is stored in mat rowwise. Defaults to true.
     * @param[out] mat The dense matrix.
     * @note The array must have space allocated for at 
     * least m*n elements, where m is the the number of rows of mat and n is 
     * the number of columns. 
     */
    template <typename T, typename InputItr>
    void reshape1DToMat(InputItr arrItr, matrix_type<T>& mat, bool byRows = true) {
        // Fail at compile time if the type T is not
        // the same as the value type of InputItr.
        static_assert(iteratorValueTypeSameAs<InputItr, T>(),
        "Matrix value type differs from array value type");
        if (byRows) {
            for (index_type i = 0; i < mat.rows(); ++i) {
                for (index_type j = 0; j < mat.cols(); ++j)
                    mat(i, j) = *arrItr++;
            }
        }
        else {
            for (index_type j = 0; j < mat.cols(); ++j) {
                for (index_type i = 0; i < mat.rows(); ++i)
                   mat(i, j) = *arrItr++;
            }
        }
    }

	/**
	 * Evaluate vector at an array of indices.
	 * @param[in] vec The input vector to be evaluated at a list of indices (or map).
	 * @param[in] map The map or list of indices.
	 * @param[out] out The output vector containing the result of applying the map.
	 */
	template <typename T, typename U>
	void applyIndexMap(const vector_type<T>& vec, const vector_type<U>& map, vector_type<T>& out) {
		static_assert(isIntegral<U>(), "map value_type is not integral");
        for (index_type k = 0; k < map.length(0); ++k)
            out(k) = vec(map(k));
    }

	/**
	 * Convert a vector to a blitz matrix.
	 * @param[in] vec The vector.
	 * @param[out] mat The output dense matrix.
	 * @param[in] byRows whether to convert row-wise (default) or column-wise.
	 */
	template <typename T>
	void vectorToFull(const vector_type<T>& vec, matrix_type<T>& mat, bool byRows = true) {
		reshape1DToMat(vec.begin(), mat, byRows);
	}

	/**
     * Convert a blitz matrix to a vector.
	 * @param[in] mat The dense matrix.
	 * @param[out] vec The output vector.
	 * @param[in] byRows whether to convert row-wise (default) or column-wise.
     */
	template <typename T>
	void fullToVector(const matrix_type<T>& mat, vector_type<T>& vec, bool byRows = true) {
		reshapeMatTo1D(mat, vec.begin(), byRows);
	}

    /**
     *  Computes Kronecker product of matrices A and B.
     *  @param[in] A The left matrix.
     *  @param[in] B The right matrix.
     *  @returns C = kron(A,B)
     */
    template <typename T>
    matrix_type<T> kron(const matrix_type<T>& A, const matrix_type<T>&B) {
        using range = blitz::Range;
        const index_type Br = B.rows(), Bc = B.cols();
        const index_type Ar = A.rows(), Ac = A.cols();
        matrix_type<T> C(Ar*Br, Ac*Bc);

        for (index_type i = 0; i < Ar; ++i) {
            for (index_type j = 0; j < Ac; ++j) {
                C(range(j * Br, (j + 1) * Br - 1), 
                range(i * Bc, (i + 1) * Bc - 1)) = A(i,j) * B;
            }
        }

        return C;
    }

    template <typename T>
    std::pair<std::vector<index_type>, std::vector<index_type>>
    unique(const matrix_type<T>& A) {
        using details::CompareExact;

        std::vector<index_type> gather(A.rows()), scatter(A.rows());
        std::iota(gather.begin(), gather.end(), index_type(0));
        CompareExact<T> comp(A);
        std::sort(gather.begin(), gather.end(), comp);
        index_type i = 0;
        for (index_type k = 0; k < A.rows(); ++k) {
            auto curr = gather[i];
            auto next = gather[k];
            if (comp(curr, next) || comp(next, curr)) { // A(next,:) != A(curr,:)
                ++i;
                gather[i] = next;
            }
            scatter[next] = i;
        }
        gather.resize(i + 1);
        return std::make_pair(std::move(gather), std::move(scatter));
    }

    template <typename T>
    std::pair<std::vector<index_type>, std::vector<index_type>>
    uniquetol(const matrix_type<T>& A, real_type tol = real_type(0)) {
        using details::CompareEQ;
        using details::CompareEQByDim;
        using details::CompareExact;
        using details::getOrdering;

        if (tol <= real_type(0)) {
            return unique(A);
        }
        std::vector<index_type> gather(A.rows()), scatter(A.rows());
        std::iota(gather.begin(), gather.end(), index_type(0));
        auto order = getOrdering(A);
        std::sort(gather.begin(), gather.end(), CompareExact<T>(A, order));
        CompareEQ<T> comp(A, tol);
        CompareEQByDim<T> compByDim(A, order[0], tol);
        index_type nunique = 0;
        auto curr = gather.begin();
        while (curr != gather.end()) {
            auto seed = *curr;
            auto first = curr;
            auto last = ++first;
            while (last != gather.end() && compByDim(seed, *last)) {
                if (!comp(seed, *last)) {
                    *first = *last;
                    ++first;
                }
                else {
                    scatter[*last] = nunique;
                }
                ++last;
            }
            scatter[seed] = nunique++;
            gather.erase(first, last);
            ++curr;
        }
        return std::make_pair(std::move(gather), std::move(scatter));
    }

    // Order the columns of A by decreasing size of their 2-norm.
    template <typename T>
    std::vector<index_type> orderColsByNorm(const matrix_type<T>& A) {
        if (A.cols() < 2) {
            return std::vector<index_type>{ index_type(0) };
        }
        std::vector<T> nrm(A.cols());
        for (index_type j = 0; j < A.cols(); ++j) {
            auto s = T(0);
            for (index_type i = 0; i < A.rows(); ++i) {
                s += A(i, j) * A(i, j);
            }
            nrm[j] = s;
        }
        std::vector<index_type> order(A.cols());
        std::iota(order.begin(), order.end(), index_type(0));
        std::sort(order.begin(), order.end(),
            [&nrm](index_type i, index_type j) { return nrm[i] > nrm[j]; });
        return order;
    }
    template <typename T>
    std::pair<std::vector<index_type>, std::vector<index_type>>
    uniquetolMatlab(const matrix_type<T>& A, real_type tol = real_type(0)) {
        using details::CompareEQ;
        using details::CompareExact;
        
        if (tol <= real_type(0)) {
            return unique(A);
        }
        std::vector<index_type> gather(A.rows()), scatter(A.rows());
        std::iota(gather.begin(), gather.end(), index_type(0));
        std::sort(gather.begin(), gather.end(), CompareExact<T>(A, orderColsByNorm(A)));
        CompareEQ<T> comp(A, tol);
        index_type i = 0;
        scatter[gather[0]] = 0;
        for (index_type k = 1; k < A.rows(); ++k) {
            auto curr = gather[k - 1];
            auto next = gather[k];
            if (!comp(curr, next)) {
                ++i;
                gather[i] = next;
            }
            scatter[next] = i;
        }
        gather.resize(i + 1);
        return std::make_pair(std::move(gather), std::move(scatter));
    }

} // namespace blitzdg