#include "BlitzHelpers.hpp"
#include "LinAlgHelpers.hpp"
#include "Types.hpp"
#include <blitz/array.h>
#include <igloo/igloo_alt.h>
#include <iostream>

using blitz::firstIndex;
using blitz::secondIndex;
using blitz::ColumnMajorArray;
using std::cout;
using std::endl;

namespace blitzdg {
    namespace DenseMatrixHelpersTests {
        using namespace igloo;

        Describe(DenseMatrixHelpers_Tests) {
            It(Should_Return_Matrix_Storage_Order) {
                real_matrix_type rmat(5,3);
                real_matrix_type cmat(4,2,ColumnMajorArray<2>());
                cout << "BlitzHelpers: storage order" << endl;
                Assert::That(isRowMajor(rmat), Equals(true));
                Assert::That(isColumnMajor(cmat), Equals(true));
            }

            It(Should_Count_The_Number_Of_Nonzeros) {
                real_matrix_type mat(5, 5);
                mat =   2,3,0,0,0,
                        3,0,4,0,6,
                        0,-1,-3,2,0,
                        0,0,1,0,0,
                        0,4,2,0,1;
                cout << "BlitzHelpers: count nonzeros" << endl;
                Assert::That(countNonzeros(mat), Equals(12));
            }

            It(Should_Convert_Full_Matrix_To_POD) {
                real_matrix_type mat(5, 5);
                mat =   2, 3, 0,0,0,
                        3, 0, 4,0,6,
                        0,-1,-3,2,0,
                        0, 0, 1,0,0,
                        0, 4, 2,0,1;
                real_vector_type pod(25);
                real_vector_type podRow(25), podCol(25);
                podRow = 2, 3, 0, 0, 0, 
                3, 0, 4, 0, 6,
                0, -1, -3, 2, 0,
                0, 0, 1, 0, 0,
                0, 4, 2, 0, 1;
                podCol = 2,3,0,0,0,3,0,-1,0,4,0,4,-3,1,2,0,0,2,0,0,0,6,0,0,1;
                reshapeMatTo1D(mat, pod.data());
                podRow -= pod;
                cout << "BlitzHelpers: full to pod (rowwise)" << endl;
                real_type diff = normInf(podRow);
                Assert::That(diff, Equals(0.0));

                reshapeMatTo1D(mat, pod.data(), false);
                podCol -= pod;
                cout << "BlitzHelpers: full to pod (colwise)" << endl;
                diff = normInf(podCol);
                Assert::That(diff, Equals(0.0));
            }

            It(Should_Convert_POD_To_Full_Matrix) {
                real_matrix_type mat(5, 5), matT(5,5);
                mat =   2,3,0,0,0,
                        3,0,4,0,6,
                        0,-1,-3,2,0,
                        0,0,1,0,0,
                        0,4,2,0,1;
                matT = 2,3,0,0,0,
                       3,0,-1,0,4,
                       0,4,-3,1,2,
                       0,0,2,0,0,
                       0,6,0,0,1;
                real_vector_type pod(25);
                pod = 2, 3, 0, 0, 0, 
                3, 0, 4, 0, 6,
                0, -1, -3, 2, 0,
                0, 0, 1, 0, 0,
                0, 4, 2, 0, 1;
                real_matrix_type ret(5,5);
                reshape1DToMat(pod.data(), ret);
                mat -= ret;
                cout << "BlitzHelpers: pod to full (rowwise)" << endl;
                real_type diff = normMax(mat);
                Assert::That(diff, Equals(0.0));

                reshape1DToMat(pod.data(), ret, false);
                
                matT -= ret;
                cout << "BlitzHelpers: pod to full (colwise)" << endl;
                diff = normMax(matT);
                Assert::That(diff, Equals(0.0));
            }

			It(Should_Convert_Dense_Matrix_To_Vector) {
				cout << "Should_Convert_Dense_Matrix_To_Vector" << endl;
				real_matrix_type mat(5, 5);
				real_vector_type vec(25), expectedvec(25);

                mat =   2,3,0,0,0,
                        3,0,4,0,6,
                        0,-1,-3,2,0,
                        0,0,1,0,0,
                        0,4,2,0,1;
				
				// Ask for column-wise vector ordering.
				const bool byRowsOpt = false;
				
				expectedvec = 2,
				3,
				0,
				0,
				0,
				3,
				0,
				-1,
				0,
				4,
				0,
				4,
				-3,
				1,
				2,
				0,
				0,
				2,
				0,
				0,
				0,
				6,
				0,
				0,
				1;
				
				fullToVector(mat, vec, byRowsOpt);

				vec -= expectedvec;
				real_type diff = normInf(vec);
				Assert::That(diff, Equals(0.0));
			}

			It(Should_Convert_Vector_To_Dense_Matrix) {
				cout << "Should_Convert_Vector_To_Dense_Matrix" << endl;
				real_vector_type vec(25);
				real_matrix_type mat(5,5);

				vec = 2,
				3,
				0,
				0,
				0,
				3,
				0,
				-1,
				0,
				4,
				0,
				4,
				-3,
				1,
				2,
				0,
				0,
				2,
				0,
				0,
				0,
				6,
				0,
				0,
				1;

				real_matrix_type expectedmat(5,5);
				expectedmat =   2,3,0,0,0,
								3,0,4,0,6,
								0,-1,-3,2,0,
								0,0,1,0,0,
								0,4,2,0,1;

				// Ask for column-wise reshaping.
				const bool byRowsOpt = false;

				vectorToFull(vec, mat, byRowsOpt);
				mat -= expectedmat;

				real_type diff = normMax(mat);
                Assert::That(diff, Equals(0.0));
			}

			It(Should_Apply_Index_Map) {
				real_vector_type vec(5);
				vec = 1.1,1.2,1.3,1.4,1.5;
				
				index_vector_type map(3);
				map = 0,2,4;

				real_vector_type expectedout(3);
				real_vector_type out(3);

				expectedout = 1.1,1.3,1.5;

				applyIndexMap(vec, map, out);

				out -= expectedout;
				real_type diff = normInf(out);
                Assert::That(diff, Equals(0.0));
			}
			It(Should_Compute_Kronecker_Product) {
				real_matrix_type A(3,3), B(3,3), C(9,9), Cexpected(9,9);
				A = 1.,0,0,
					0,1.,0,
					0,0,1.;
				B = 3.,3.,3.,
					3.,3.,3.,
					3.,3.,3.;
				
				C = kron(A,B);

				Cexpected = 3.,3.,3.,0.,0.,0.,0.,0.,0.,
					3.,3.,3.,0.,0.,0.,0.,0.,0.,
					3.,3.,3.,0.,0.,0.,0.,0.,0.,
					0.,0.,0.,3.,3.,3.,0.,0.,0.,
					0.,0.,0.,3.,3.,3.,0.,0.,0.,
					0.,0.,0.,3.,3.,3.,0.,0.,0.,
					0.,0.,0.,0.,0.,0.,3.,3.,3.,
					0.,0.,0.,0.,0.,0.,3.,3.,3.,
					0.,0.,0.,0.,0.,0.,3.,3.,3.;

				C -= Cexpected;
				real_type diff = normMax(C);
                Assert::That(diff, Equals(0.0));
			}

			It(Should_Find_Unique_Rows_With_Tol) {
				real_matrix_type A(6, 3);
				A = 1.,1.,1,
					1.,1.,1,
					2., 0., 3.,
					2.00000001, 0., 3,
					2.,2.,2.,
					2.,2.,2.1;

				std::pair<std::vector<index_type>,std::vector<index_type>> uniqueIndices;
				uniqueIndices = uniquetol(A, 1e-4);
				
				auto gather = uniqueIndices.first;
				auto scatter = uniqueIndices.second;


				real_matrix_type AuniqExpected(4, 3);
				AuniqExpected = 1.,1.,1,
								2.,2.,2.,
								2.,2.,2.1,
								2.,0.,3.;

				real_matrix_type Auniq(gather.size(), 3);

				real_matrix_type Ascatter(scatter.size(), 3);

				for (index_type i=0; i < static_cast<index_type>(gather.size()); ++i) {
					Auniq(i, blitz::Range::all()) = A(gather[i], blitz::Range::all());
				}

				AuniqExpected -= Auniq;
				Assert::That(normFro(AuniqExpected), IsLessThan(1e-5));

				for (index_type i=0; i < static_cast<index_type>(scatter.size()); ++i) {
					Ascatter(i, blitz::Range::all()) = Auniq(scatter[i], blitz::Range::all());
				}

				Ascatter -= A;
				Assert::That(normFro(Ascatter), IsLessThan(1e-5));
			}
			It(Should_Find_Unique_Rows_With_Tol_Matlab_Style_Impl) {
				real_matrix_type A(6, 3);
				A = 1.,1.,1,
					1.,1.,1,
					2., 0., 3.,
					2.00000001, 0., 3,
					2.,2.,2.,
					2.,2.,2.1;

				std::pair<std::vector<index_type>,std::vector<index_type>> uniqueIndices;
				uniqueIndices = uniquetolMatlab(A, 1e-4);
				
				auto gather = uniqueIndices.first;
				auto scatter = uniqueIndices.second;


				real_matrix_type AuniqExpected(4, 3);
				AuniqExpected = 1.,1.,1,
								2.,2.,2.,
								2.,2.,2.1,
								2.,0.,3.;

				real_matrix_type Auniq(gather.size(), 3);

				real_matrix_type Ascatter(scatter.size(), 3);

				for (index_type i=0; i < static_cast<index_type>(gather.size()); ++i) {
					Auniq(i, blitz::Range::all()) = A(gather[i], blitz::Range::all());
				}

				AuniqExpected -= Auniq;
				Assert::That(normFro(AuniqExpected), IsLessThan(1e-5));

				for (index_type i=0; i < static_cast<index_type>(scatter.size()); ++i) {
					Ascatter(i, blitz::Range::all()) = Auniq(scatter[i], blitz::Range::all());
				}

				Ascatter -= A;
				Assert::That(normFro(Ascatter), IsLessThan(1e-5));
			}
        };
    } // namespace DenseMatrixHelpersTests
} // namespace blitzdg