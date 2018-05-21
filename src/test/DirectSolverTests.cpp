// Copyright (C) 2017-2018  Waterloo Quantitative Consulting Group, Inc.
// See COPYING and LICENSE files at project root for more details.


#include "DirectSolver.hpp"
#include "SparseMatrixConverter.hpp"
#include "Types.hpp"
#include <igloo/igloo_alt.h>
#include <blitz/array.h>
#include <iostream>
#include <limits>

using blitz::firstIndex;
using blitz::Range;
using blitz::secondIndex;
using std::cout;
using std::endl;
using std::numeric_limits;

namespace blitzdg {
    namespace DirectSolverTests {
        using namespace igloo;
        const index_type N=5;
        const real_type eps=10*numeric_limits<real_type>::epsilon();
        const float epsf = 1.e-5;
        
        firstIndex ii;
        secondIndex jj;

        matrix_type b(N,1), x(N,1);
        matrix_type A(N,N);
        vector_type expectedx(N);

        Describe(DirectSolver_Object) {
            void SetUp() {
                cout << "DirectSolver" << endl;

                A = -1.550762,0.524638,1.483705,-2.354095,1.397928,
                0.370182,-0.747854,-0.020957,-1.120856,0.121897,
                0.343633,0.855609,-0.206901,0.713489,1.442008,
                0.825662,-1.416542,0.246512,-0.884317,1.672087,
                0.014445,1.065063,0.151608,1.221347,-0.280884;

                b(Range::all(), 0) = 2.36882,
                                    0.29768,
                                    0.35453,
                                    0.88904,
                                    -1.22226;
                expectedx = -1.32806,
                            -0.38880,
                            -0.96358,
                            -0.33740,
                            0.82171;
            }

            It(Should_Solve_Random_Linear_System) {
                cout << "It Should Solve Random Linear System" << endl;
                
                DirectSolver solver;
                solver.solve(A, b, x);

                Assert::That(fabs(x(0)-expectedx(0)), IsLessThan(epsf));
                Assert::That(fabs(x(1)-expectedx(1)), IsLessThan(epsf));
                Assert::That(fabs(x(2)-expectedx(2)), IsLessThan(epsf));
                Assert::That(fabs(x(3)-expectedx(3)), IsLessThan(epsf));
                Assert::That(fabs(x(4)-expectedx(4)), IsLessThan(epsf));
            }
        };
    } // namespace DirectSolverTests
} // namespace blitzdg