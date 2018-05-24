#include "CSCMatrix.hpp"
#include <cmath>
#include <iomanip>
#include <limits>
#include <stdexcept>

using std::abs;
using std::move;
using std::numeric_limits;
using std::ostream;
using std::runtime_error;
using std::setw;

namespace blitzdg {
    namespace {
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

        index_type numDigits(index_type n) {
            if (n == 0) return 1;
            index_type ret = 0;
            while (n > 0) {
                ++ret;
                n /= 10;
            }
            return ret;
        }
    } // anonymous namespace

    /////////////////////////////////////
    // MEMBER FUNCTIONS
    /////////////////////////////////////
    CSCMat::CSCMat(index_type rows, index_type cols, index_type nnz)
		: mat_{ cs_spalloc(rows, cols, (rows * cols == 0) ? 0 : nnz, 1, 0) }
	{
		if (!mat_)
			throw runtime_error("CSCMat::CSCMat: matrix construction failed");
        if (rows == 0 || cols == 0) {
            mat_->m = 0;
            mat_->n = 0;
        }
	}

    CSCMat::CSCMat(index_type rows, index_type cols, const SparseTriplet& triplet) 
        : mat_{ nullptr }
    {
        // create a CXSparse triplet that is a copy of triplet
        cs_di* tmp = cs_di_spalloc(rows, cols, triplet.nz, 1, 1);
        if (!tmp)
            throw runtime_error("CSCMat::CSCMat: unable to create matrix from sparse triplet");
        for (index_type k = 0; k < triplet.nz; ++k) {
            if (!cs_di_entry(tmp, triplet.row[k], triplet.col[k], triplet.val[k])) {
                cs_di_spfree(tmp);
                throw runtime_error("CSCMat::CSCMat: unable to create matrix from sparse triplet");
            }
        }
        // create mat_ by compressing the triplet to a CSC matrix
        mat_.reset(cs_di_compress(tmp));
        // free the temporary storage
        cs_di_spfree(tmp);
        if (!mat_)
            throw runtime_error("CSCMat::CSCMat: unable to create matrix from sparse triplet");
    }

    CSCMat::CSCMat(const matrix_type& mat, real_type dropTol) 
        : mat_{ nullptr }
    {
        // allocate memory
        index_type nnz = countNonzeros(mat, dropTol);
        mat_.reset(cs_di_spalloc(mat.rows(), mat.cols(), nnz, 1, 0));
        if (!mat_)
            throw runtime_error("CSCMat::CSCMat: unable to create matrix from dense matrix");
        // assign values
        // Note: this may be slow if mat is stored by rows in memory.
        nnz = 0;
        for (index_type j = 0; j < mat.cols(); ++j) {
            mat_->p[j] = nnz;
            for (index_type i = 0; i < mat.rows(); ++i) {
                real_type elem = mat(i, j);
                if (abs(elem) > dropTol) {
                    mat_->i[nnz] = i;
                    mat_->x[nnz++] = elem;
                }
            }
        }
        mat_->p[mat.cols()] = nnz;
    }
	
	CSCMat::CSCMat(cs_di_smart_ptr mat)
		: mat_{ move(mat) }
	{
		if (!mat_)
			throw runtime_error("CSCMat::CSCMat: input matrix is null");
		if (!mat_->x)
			throw runtime_error("CSCMat::CSCMat: values array in input matrix is null");
        // if mat_ is a triplet, then convert to csc
        if (mat_->nz >= 0) {
            mat_.reset(cs_di_compress(mat_.get()));
            if (!mat_)
                throw runtime_error("CSCMat::CSCMat: failed to build from smart pointer");
        }
	}

	CSCMat::CSCMat(const CSCMat& other)
		: CSCMat(other.rows(), other.cols(), other.nnz())
	{
		for (index_type k = 0; k < other.nnz(); ++k) {
			mat_->i[k] = other.mat_->i[k];
			mat_->x[k] = other.mat_->x[k];
		}
		for (index_type k = 0; k <= other.cols(); ++k)
			mat_->p[k] = other.mat_->p[k];
	}

	CSCMat& CSCMat::operator=(CSCMat other) {
		swap(*this, other);
		return *this;
	}

    void CSCMat::removeDuplicates() {
        if (!cs_di_dupl(mat_.get()))
            throw runtime_error("CSCMat::removeDuplicates: failed");
    }

    void CSCMat::prune(real_type dropTol) {
        if (!cs_di_droptol(mat_.get(), dropTol))
            throw runtime_error("CSCMat::prune: failed");
    }

    void CSCMat::transpose() {
        cs_di* trsp = cs_di_transpose(mat_.get(), 1);
        if (!trsp)
            throw runtime_error("CSCMat::tranpose: failed");
        mat_.reset(trsp);
    }
		
    /////////////////////////////////////
    // FREE FUNCTIONS
    /////////////////////////////////////
	void swap(CSCMat& lhs, CSCMat& rhs) {
		using std::swap;
		swap(lhs.mat_->m, rhs.mat_->m);
		swap(lhs.mat_->n, rhs.mat_->n);
		swap(lhs.mat_->nzmax, rhs.mat_->nzmax);
		swap(lhs.mat_->nz, rhs.mat_->nz);
		swap(lhs.mat_->i, rhs.mat_->i);
		swap(lhs.mat_->p, rhs.mat_->p);
		swap(lhs.mat_->x, rhs.mat_->x);
	}

    std::ostream& operator<<(std::ostream& strm, const CSCMat& mat) {
        const index_type ndr = numDigits(mat.rows());
        const index_type ndc = numDigits(mat.cols());
        strm << "rows = " << mat.rows() << ", " 
            << "cols = " << mat.cols() << ", " 
            << "nnz = " << mat.nnz() << "\n\n";
        for (index_type j = 0; j < mat.cols(); ++j) {
            for (index_type k = mat.colPtrs(j); k < mat.colPtrs(j + 1); ++k) {
                strm << setw(ndr) << mat.rowInds(k) << " " 
                    << setw(ndc) << j << " " 
                    << mat.elems(k) << "\n";
            }
        }
        return strm;
    }
	
	CSCMat transpose(const CSCMat& mat) {
		cs_di* tmp = cs_di_transpose(mat.matPtr(), 1);
		if (!tmp)
			throw runtime_error("CSCMat matrix transpose failed");
		return CSCMat(CSCMat::cs_di_smart_ptr{ tmp });
	}
	
	CSCMat multiply(const CSCMat& lhs, const CSCMat& rhs) {
		cs_di* tmp = cs_di_multiply(lhs.matPtr(), rhs.matPtr());
		if (!tmp)
			throw runtime_error("CSCMat matrix-matrix multiplication failed");
		return CSCMat(CSCMat::cs_di_smart_ptr{ tmp });
	}
} // namespace blitzdg