// Copyright (C) 2017-2018  Derek Steinmoeller. 
// See COPYING and LICENSE files at project root for more details. 

#include <igloo/igloo_alt.h>
#include <blitz/array.h>
#include <LUSolver.hpp>
#include <SparseMatrixConverter.hpp>

using namespace igloo;
using namespace blitz;
using namespace std;

namespace LUSolverTests {
    const int N=5;
    const double eps=10*numeric_limits<double>::epsilon();

    Array<double, 1> b(N), x(N);
    Array<double, 2> A(N,N);

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
            luSolver = new LUSolver(&A, *matrixConverter);
        }

        It(Solves_Ax_equals_b)  {
            cout << "LUSolver" << endl;
            LUSolver & solver = *luSolver;
            Array <double, 1> soln(N);

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
}