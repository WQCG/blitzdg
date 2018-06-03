#include "SparseTriplet.hpp" 
#include "BlitzHelpers.hpp"
#include <cmath>
#include <iomanip>
#include <limits>
#include <stdexcept>

using std::abs;
using std::numeric_limits;
using std::ostream;
using std::runtime_error;
using std::setw;

namespace blitzdg {
	namespace {
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

	index_type SparseTriplet::newSize() const {
		real_type nzmax_new = nzmax_ < 2 ? 2.0 : 1.5 * nzmax_;
		if (nzmax_new > numeric_limits<index_type>::max())
			throw runtime_error("SparseTriplet::newSize: matrix capacity exceeds maximum allowable");
		return static_cast<index_type>(nzmax_new);
	}

	void SparseTriplet::grow(index_type nzmax_new) {
		if (nzmax_new <= nzmax_)
			return;
		row_.resize(nzmax_new);
		col_.resize(nzmax_new);
		elems_.resize(nzmax_new);
		nzmax_ = nzmax_new;
	}

	SparseTriplet::SparseTriplet(const real_matrix_type& mat, real_type dropTol) 
		: SparseTriplet(mat.rows(), mat.cols(), countNonzeros(mat, dropTol))
	{
		for (real_matrix_type::const_iterator itr = mat.begin(); itr != mat.end(); ++itr) {
			if (abs(*itr) > dropTol) {
				row_[nnz_] = itr.position()[0];
				col_[nnz_] = itr.position()[1];
				elems_[nnz_++] = *itr;
			}
		}
	}

	void swap(SparseTriplet& lhs, SparseTriplet& rhs) {
		using std::swap;
		swap(lhs.rows_, rhs.rows_);
		swap(lhs.cols_, rhs.cols_);
		swap(lhs.nnz_, rhs.nnz_);
		swap(lhs.nzmax_, rhs.nzmax_);
		swap(lhs.row_, rhs.row_);
		swap(lhs.col_, rhs.col_);
		swap(lhs.elems_, rhs.elems_);
	}

	ostream& operator<<(ostream& strm, const SparseTriplet& mat) {
		const index_type ndr = numDigits(mat.rows());
		const index_type ndc = numDigits(mat.cols());
		strm << "rows = " << mat.rows() << ", " 
			<< "cols = " << mat.cols() << ", " 
			<< "nnz = " << mat.nnz() << "\n\n";
		for (index_type k = 0; k < mat.nnz(); ++k) {
			strm << setw(ndr) << mat.row(k) << " " 
				<< setw(ndc) << mat.col(k) << " " 
				<< mat.elem(k) << "\n";
		}
		return strm;
	}
} // namespace blitzdg