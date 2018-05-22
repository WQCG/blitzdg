// Copyright (C) 2017-2018  Waterloo Quantitative Consulting Group, Inc.
// See COPYING and LICENSE files at project root for more details.

#include "LUSolver.hpp"
#include "SparseMatrixConverter.hpp"
#include "Types.hpp"
#include <igloo/igloo_alt.h>
#include <blitz/array.h>
#include <iostream>
#include <limits>

using std::cout;
using std::endl;
using std::numeric_limits;
using std::abs;

namespace blitzdg {
    namespace LUSolverTests {
        using namespace igloo;
        const index_type N=5;
        const real_type eps=10*numeric_limits<real_type>::epsilon();

        vector_type b(N), x(N);
        matrix_type A(N,N);

        SparseMatrixConverter * matrixConverter = nullptr;
        LUSolver * luSolver = nullptr;

        Describe(LUSolver_Object) {
            void SetUp() {

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

                matrixConverter = new SparseMatrixConverter();
                luSolver = new LUSolver(&A);
            }

            It(Solves_Ax_equals_b)  {
                cout << "LUSolver" << endl;
                LUSolver & solver = *luSolver;
                vector_type soln(N);

                // Compute LU factors.
                solver.factorize();
                solver.solve(b, soln);

                Assert::That(abs(soln(0)-x(0)), IsLessThan(eps));
                Assert::That(abs(soln(1)-x(1)), IsLessThan(eps));
                Assert::That(abs(soln(2)-x(2)), IsLessThan(eps));
                Assert::That(abs(soln(3)-x(3)), IsLessThan(eps));
                Assert::That(abs(soln(4)-x(4)), IsLessThan(eps));
            }

            void TearDown() {
                delete luSolver;
                delete matrixConverter;
            }
        };
    } // namespace LUSolverTests
} // namespace blitzdg
