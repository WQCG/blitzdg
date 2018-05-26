// Copyright (C) 2017-2018  Waterloo Quantitative Consulting Group, Inc.
// See COPYING and LICENSE files at project root for more details.

#include "LUSolver.hpp"
#include "Types.hpp"
#include <igloo/igloo_alt.h>
#include <blitz/array.h>
#include <iostream>
#include <limits>
#include <stdexcept>

using std::cout;
using std::endl;
using std::exception;
using std::numeric_limits;
using std::abs;

namespace blitzdg {
    namespace LUSolverTests {
        using namespace igloo;
        Describe(LUSolver_Object) {
            It(Solves_Ax_equals_b)  {
                const real_type eps=10*numeric_limits<real_type>::epsilon();
                const index_type N=5;
                vector_type b(N), x(N);
                matrix_type A(N,N);
                 A = 2,3,0,0,0,
                        3,0,4,0,6,
                        0,-1,-3,2,0,
                        0,0,1,0,0,
                        0,4,2,0,1;

                b =  8,
                    45,
                    -3,
                    3,
                    19;

                x = 1,
                    2,
                    3,
                    4,
                    5;
                cout << "LUSolver" << endl;
                CSCMat csc(A, 2*numeric_limits<real_type>::epsilon());
                LUSolver solver(csc);
                cout << "Done building LUSolver" << endl;
                vector_type soln(N);

                // Compute LU factors.
                try {
                    cout << "Starting factorization" << endl;
                    solver.factorize();
                    cout << "Factorization successful" << endl;
                    solver.solve(b, soln);
                    cout << "Solve successful" << endl;
                }
                catch (exception& e) {
                    cout << "Excpetion caught: " << e.what() << endl;
                }

                Assert::That(abs(soln(0)-x(0)), IsLessThan(eps));
                Assert::That(abs(soln(1)-x(1)), IsLessThan(eps));
                Assert::That(abs(soln(2)-x(2)), IsLessThan(eps));
                Assert::That(abs(soln(3)-x(3)), IsLessThan(eps));
                Assert::That(abs(soln(4)-x(4)), IsLessThan(eps));
            }
        };
    } // namespace LUSolverTests
} // namespace blitzdg
