#include "DenseMatrixInverter.hpp"
#include "LinAlgHelpers.hpp"
#include "Types.hpp"
#include <blitz/array.h>
#include <igloo/igloo_alt.h>
#include <iostream>

using blitz::firstIndex;
using blitz::secondIndex;
using blitz::ColumnMajorArray;
using blitz::sum;
using blitzdg::DenseMatrixInverter;
using std::cout;
using std::endl;

namespace blitzdg {
    namespace DenseMatrixHelpersTests {
        using namespace igloo;

        Describe(DenseMatrixInverter_Object) {
            const real_type epsf = 5.e-6;
            It(Should_Invert_Matrix) {
                cout << "DenseMatrixInverter: should invert " << endl;

                DenseMatrixInverter inverter;
                real_matrix_type mat(4,4), matinv(4,4), expectedmatinv(4,4);

                mat = -0.0030118,-0.1275669,-1.6166021,-2.5901948,
                       1.2757077,-1.1244873,-0.3668178,1.0018966,
                      -0.1412610,1.7840822,-0.4695502,-0.2113200,
                      -0.6717781,-0.4816790,-0.1003273,1.3764271;

                inverter.computeInverse(mat, matinv);

                expectedmatinv = -0.122272,0.499965,0.152316,-0.570633,
                                 -0.129222,-0.0241404,0.495695,-0.149498,
                                 -0.394064,-0.337091,-0.390806,-0.556192,
                                 -0.13362,0.210994,0.219321,0.355159;

                real_matrix_type res(4,4);
                res = matinv - expectedmatinv;
                Assert::That(sqrt(sum(res*res)), IsLessThan(epsf));
            }
        };
    } // namespace DenseMatrixInverterTests
} // namespace blitzdg