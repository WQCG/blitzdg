#include "CSCMatrix.hpp"
#include "LinAlgHelpers.hpp"
#include "SparseMatrixConverter.hpp"
#include <blitz/array.h>
#include <igloo/igloo_alt.h>
#include <iostream>

using blitz::firstIndex;
using blitz::secondIndex;
using blitz::ColumnMajorArray;
using std::cout;
using std::endl;
using std::move;

namespace blitzdg {
    namespace CSCMatrixTests {
        using namespace igloo;

        matrix_type cscToFull(const CSCMat& mat) {
            matrix_type ret(mat.rows(), mat.cols(), ColumnMajorArray<2>());
            ret = real_type(0);
            for (index_type j = 0; j < mat.cols(); ++j) {
                for (index_type k = mat.colPtrs(j); k < mat.colPtrs(j + 1); ++k)
                    ret(mat.rowInds(k), j) = mat.elems(k);
            }
            return ret;
        }

        Describe(CSCMat_Object) {
            It(Should_Build_From_A_Dense_Matrix) {
                matrix_type full(5,5);
                full = 1, 0, 0, 3, 2,
                       0, 0, 1, 1, 0,
                       0, 5, 0, 0, 0,
                       0, 2, 0, 0, 4,
                       1, 1, 1, 0, 0;
                CSCMat csc(full);
                matrix_type ret(5,5);
                ret = cscToFull(csc);
                ret -= full;
                real_type diff = normMax(ret);
                cout << "CSCMat: build from full matrix" << endl;
                Assert::That(csc.rows(), Equals(5));
                Assert::That(csc.cols(), Equals(5));
                Assert::That(csc.nnz(), Equals(11));
                Assert::That(csc.empty(), Equals(false));
                Assert::That(diff, Equals(0.0));
            }

            It(Should_Build_From_A_Sparse_Triplet) {
                matrix_type full(5,5);
                full = 1, 0, 0, 3, 2,
                       0, 0, 1, 1, 0,
                       0, 5, 0, 0, 0,
                       0, 2, 0, 0, 4,
                       1, 1, 1, 0, 0;
                SparseTriplet trip;
                SparseMatrixConverter conv;
                conv.fullToSparseTriplet(full, trip);
                CSCMat csc(5, 5, trip);
                matrix_type ret(5,5);
                ret = cscToFull(csc);
                ret -= full;
                real_type diff = normMax(ret);
                cout << "CSCMat: build from triplet matrix" << endl;
                Assert::That(diff, Equals(0.0));
            }

            It(Should_Build_From_A_Cs_Di_Pointer) {
                cs_di* tmp = cs_di_spalloc(5, 5, 25, 1, 1);
                index_type k = 0;
                for (index_type j = 0; j < 5; ++j) {
                    tmp->p[j] = k;
                    for (index_type i = 0; i < 5; ++i) {
                        tmp->x[k] = 1;
                        tmp->i[k++] = i;
                    }
                }
                tmp->p[5] = k;
                tmp->nz = -1;
                CSCMat::cs_di_smart_ptr ptr(tmp);
                CSCMat csc(move(ptr));
                matrix_type full(5,5), ret(5,5);
                full = 1,1,1,1,1,
                       1,1,1,1,1,
                       1,1,1,1,1,
                       1,1,1,1,1,
                       1,1,1,1,1;
                ret = cscToFull(csc);
                ret -= full;
                real_type diff = normMax(ret);
                cout << "CSCMat: build from cs_di_smart_ptr" << endl;
                Assert::That(diff, Equals(0.0));
            }

            It(Should_Prune_Any_Explicit_Zeros) {
                matrix_type full(5,5);
                full = 1, 0, 0, 3, 2,
                       0, 0, 1, 1, 0,
                       0, 5, 0, 0, 0,
                       0, 2, 0, 0, 4,
                       1, 1, 1, 0, 0;
                CSCMat csc(full, -1);
                Assert::That(csc.nnz(), Equals(25));
                csc.prune(-1);
                matrix_type ret(5,5);
                ret = cscToFull(csc);
                ret -= full;
                real_type diff = normMax(ret);
                cout << "CSCMat: prune zero elements from matrix" << endl;
                Assert::That(diff, Equals(0.0));
            }

            It(Should_Remove_Duplicate_Elements) {
                matrix_type full(5,5);
                full = 1, 0, 0, 3, 2,
                       0, 0, 1, 1, 0,
                       0, 5, 0, 0, 0,
                       0, 2, 0, 0, 4,
                       1, 1, 1, 0, 0;
                SparseTriplet trip;
                SparseMatrixConverter conv;
                conv.fullToSparseTriplet(full, trip);
                trip.row[0] = 0; trip.col[0] = 0; trip.val[0] = 1;
                trip.row[1] = 0; trip.col[1] = 0; trip.val[1] = 2;
                trip.row[5] = 0; trip.col[5] = 0; trip.val[5] = 3;
                trip.row[9] = 0; trip.col[9] = 0; trip.val[9] = 4;
                CSCMat csc(5, 5, trip);
                csc.removeDuplicates();
                matrix_type ret(5,5);
                ret = cscToFull(csc);
                cout << "CSCMat: remove duplicate elements from matrix" << endl;
                Assert::That(csc.nnz(), Equals(8));
                Assert::That(ret(0,0), Equals(10));
            }

            It(Should_Compute_The_Transpose) {
                matrix_type full(5,5);
                full = 1, 0, 0, 3, 2,
                       0, 0, 1, 1, 0,
                       0, 5, 0, 0, 0,
                       0, 2, 0, 0, 4,
                       1, 1, 1, 0, 0;
                CSCMat csc(full);
                CSCMat csct = transpose(csc);
                firstIndex ii;
                secondIndex jj;
                matrix_type trans(5,5);
                trans = full(jj,ii);
                matrix_type ret(5,5);
                ret = cscToFull(csct);
                ret -= trans;
                real_type diff = normMax(ret);
                cout << "CSCMat: tranpose" << endl;
                Assert::That(diff, Equals(0.0));
                csct.transpose();
                ret = cscToFull(csct);
                ret -= full;
                diff = normMax(ret);
                Assert::That(diff, Equals(0.0));
            }

            It(Should_Multiply_Two_Matrices) {
                matrix_type full(5,5), fullp(5,5);
                full = 1, 0, 0, 3, 2,
                       0, 0, 1, 1, 0,
                       0, 5, 0, 0, 0,
                       0, 2, 0, 0, 4,
                       1, 1, 1, 0, 0;
                fullp =   2,    1,    1,    3,    2,
                          1,   30,    1,    0,    8,
                          1,    1,    2,    1,    0,
                          3,    0,    1,   10,    6,
                          2,    8,    0,    6,   20;
                CSCMat csc(full);
                CSCMat csct = transpose(csc);
                CSCMat prod = multiply(csct, csc);
                matrix_type ret(5,5);
                ret = cscToFull(prod);
                ret -= fullp;
                real_type diff = normMax(ret);
                cout << "CSCMat: matrix-matrix multiplication" << endl;
                Assert::That(diff, Equals(0.0));
            }

            It(Should_Copy_Construct) {
                matrix_type full(5,5);
                full = 1, 0, 0, 3, 2,
                       0, 0, 1, 1, 0,
                       0, 5, 0, 0, 0,
                       0, 2, 0, 0, 4,
                       1, 1, 1, 0, 0;
                CSCMat csc(full);
                CSCMat cscnew = csc;
                matrix_type ret(5,5);
                ret = cscToFull(cscnew);
                ret -= full;
                real_type diff = normMax(ret);
                cout << "CSCMat: copy constructor" << endl;
                Assert::That(diff, Equals(0.0));
            }

            It(Should_Copy_Assign) {
                matrix_type full(5,5), dens(3,3);
                full = 1, 0, 0, 3, 2,
                       0, 0, 1, 1, 0,
                       0, 5, 0, 0, 0,
                       0, 2, 0, 0, 4,
                       1, 1, 1, 0, 0;
                dens = 1, 2, 3,
                       4, 5, 6,
                       7, 8, 9;
                CSCMat cscf(full);
                CSCMat cscd(dens);
                cscf = cscd;
                matrix_type ret(3,3);
                ret = cscToFull(cscf);
                ret -= dens;
                real_type diff = normMax(ret);
                cout << "CSCMat: copy assignment" << endl;
                Assert::That(diff, Equals(0.0));
            }
        };
    } // namespace CSCMatrixTests
} // namespace blitzdg
