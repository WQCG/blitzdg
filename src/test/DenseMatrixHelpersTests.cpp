#include "DenseMatrixHelpers.hpp"
#include "LinAlgHelpers.hpp"
#include "Types.hpp"
#include <igloo/igloo_alt.h>
#include <iostream>

using blitz::ColumnMajorArray;
using std::cout;
using std::endl;

namespace blitzdg {
    namespace DenseMatrixHelpersTests {
        using namespace igloo;

        Describe(DenseMatrixHelpers_Tests) {
            It(Should_Return_Matrix_Storage_Order) {
                matrix_type rmat(5,3);
                matrix_type cmat(4,2,ColumnMajorArray<2>());
                cout << "DenseMatrixHelpers: storage order" << endl;
                Assert::That(isRowMajor(rmat), Equals(true));
                Assert::That(isColumnMajor(cmat), Equals(true));
            }

            It(Should_Count_The_Number_Of_Nonzeros) {
                matrix_type mat(5, 5);
                mat =   2,3,0,0,0,
                        3,0,4,0,6,
                        0,-1,-3,2,0,
                        0,0,1,0,0,
                        0,4,2,0,1;
                cout << "DenseMatrixHelpers: count nonzeros" << endl;
                Assert::That(countNonzeros(mat), Equals(12));
            }

            It(Should_Convert_Full_Matrix_To_POD) {
                matrix_type mat(5, 5);
                mat =   2,3,0,0,0,
                        3,0,4,0,6,
                        0,-1,-3,2,0,
                        0,0,1,0,0,
                        0,4,2,0,1;
                vector_type pod(25);
                vector_type podTrue(25);
                podTrue = 2, 3, 0, 0, 0, 
                3, 0, 4, 0, 6,
                0, -1, -3, 2, 0,
                0, 0, 1, 0, 0,
                0, 4, 2, 0, 1;
                fullToPodArray(mat, pod.data());
                podTrue -= pod;
                cout << "DenseMatrixHelpers: full to pod" << endl;
                real_type diff = normInf(podTrue);
                Assert::That(diff, Equals(0.0));
            }

            It(Should_Convert_POD_To_Full_Matrix) {
                matrix_type mat(5, 5);
                mat =   2,3,0,0,0,
                        3,0,4,0,6,
                        0,-1,-3,2,0,
                        0,0,1,0,0,
                        0,4,2,0,1;
                vector_type pod(25);
                pod = 2, 3, 0, 0, 0, 
                3, 0, 4, 0, 6,
                0, -1, -3, 2, 0,
                0, 0, 1, 0, 0,
                0, 4, 2, 0, 1;
                matrix_type ret(5,5);
                podArrayToFull(pod.data(), ret);
                mat -= ret;
                cout << "DenseMatrixHelpers: pod to full" << endl;
                real_type diff = normMax(mat);
                Assert::That(diff, Equals(0.0));
            }
        };
    } // namespace DenseMatrixHelpersTests
} // namespace blitzdg