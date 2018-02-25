#include <igloo/igloo_alt.h>
#include <blitz/array.h>
#include <EigenSolver.hpp>
#include <SparseMatrixConverter.hpp>

using namespace igloo;
using namespace std;
using namespace blitz;

namespace EigenSolverTests {
    const int N=5;
    const double eps=10*numeric_limits<double>::epsilon();
    const float epsf = 5.e-7;

    EigenSolver * eigenSolver = nullptr;
    SparseMatrixConverter * matrixConverter = nullptr;
    
    firstIndex ii;
    secondIndex jj;

    Array<double, 1> b(N), x(N);
    Array<double, 2> Adiag(N, N), Asymmetric(N, N);

    Describe(EigenSolver_Object) {
        void SetUp() {
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

            Adiag = 1,0,0,0,0,
                0,2,0,0,0,
                0,0,3,0,0,
                0,0,0,4,0,
                0,0,0,0,5;

            Asymmetric = 0,0.57735,0,0,0,
                        0.57735,-0,0.516398,0,0,
                        0,0.516398,-0,0.507093,0,
                        0,0,0.507093,-0,0.503953,
                        0,0,0,0.503953,-0;

            matrixConverter = new SparseMatrixConverter();
        }

        It(Should_Solve_Trivial_Symmetric_Eigenproblem) {
            eigenSolver = new EigenSolver(*matrixConverter);
            EigenSolver & solver = *eigenSolver;

            Array<double, 1> eigenvalues(5);
            eigenvalues = 0,0,0,0,0;
            Array<double, 2> eigenvectors(5,5);
            eigenvectors = 0,0,0,0,0,
                        0,0,0,0,0,
                        0,0,0,0,0,
                        0,0,0,0,0,
                        0,0,0,0,0;

            solver.solve(Adiag, eigenvalues, eigenvectors);
            
            Array<double, 2> expectedEvecs(5,5);
            
            expectedEvecs = 1,0,0,0,0,
                            0,1,0,0,0,
                            0,0,1,0,0,
                            0,0,0,1,0,
                            0,0,0,0,1;

            Assert::That(eigenvalues(0), Equals(1.));
            Assert::That(eigenvalues(1), Equals(2.));
            Assert::That(eigenvalues(2), Equals(3.));
            Assert::That(eigenvalues(3), Equals(4.));
            Assert::That(eigenvalues(4), Equals(5.));

            Array <double, 2> res(5,5);
            res = eigenvectors - expectedEvecs;
            Assert::That(sum(res(ii)*res(ii)), IsLessThan(eps));
        }

        It(Should_Solve_NonTrivial_Symmetric_Eigenproblem) {
            eigenSolver = new EigenSolver(*matrixConverter);
            EigenSolver & solver = *eigenSolver;

            Array<double, 1> eigenvalues(5);
            eigenvalues = 0,0,0,0,0;
            Array<double, 2> eigenvectors(5,5);
            eigenvectors = 0,0,0,0,0,
                        0,0,0,0,0,
                        0,0,0,0,0,
                        0,0,0,0,0,
                        0,0,0,0,0;

            solver.solve(Asymmetric, eigenvalues, eigenvectors);
            
            Array<double, 2> expectedEvecs(5,5);

            expectedEvecs = 0.344185,-0.540215,0.563165,-0.456254,0.253736,
                        -0.489198,0.456254,0.0711849,-0.540215,0.505587,
                            0.533334,2.06417e-16,-0.596285,-2.67144e-16,0.6,
                        -0.489198,-0.456254,0.0711849,0.540215,0.505587,
                            0.344185,0.540215,0.563165,0.456254,0.253736;

            Assert::That(eigenvalues(0) - -0.90618,    IsLessThan(epsf));
            Assert::That(eigenvalues(1) - -0.538469,   IsLessThan(epsf));
            Assert::That(eigenvalues(2) - 9.62592e-17, IsLessThan(epsf));
            Assert::That(eigenvalues(3) - 0.538469,    IsLessThan(epsf));
            Assert::That(eigenvalues(4) - 0.90618,     IsLessThan(epsf));

            Array<double, 2> res(5,5);
            res = eigenvectors - expectedEvecs;
            Assert::That(sum(res(ii)*res(ii)), IsLessThan(epsf));
        }
    };
}