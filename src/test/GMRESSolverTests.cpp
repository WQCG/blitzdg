#include "GMRESSolver.hpp"
#include <igloo/igloo_alt.h>
#include <iostream>

using namespace igloo;
using namespace blitzdg;

using std::cout;

namespace GMRESSolverTests {
    class Matvec {
        const matrix_type& mref_;
    public:
        Matvec(const matrix_type& mat)
            : mref_{ mat }
        {}

        bool operator()(const vector_type& in, vector_type& out) const {
            blitz::firstIndex i;
            blitz::secondIndex j;
            out = blitz::sum(mref_(i, j) * in(j), j);
            return true;
        }
    };

    // identity preconditioner
    class Precon {
    public:
        bool operator()(const vector_type& in, vector_type& out) const {
            out = in;
            return true;
        }
    };

    GMRESSolver* solver = nullptr;

    Describe(GMRESSolver_Object) {
        void SetUp() {
            solver = new GMRESSolver();
        }

        It(Solve_Ax_equals_b) {
            index_type N = 5;
            vector_type b(N), x(N);
            matrix_type A(N, N);

            A = 2.,3.,0.,0.,0.,
                3.,0.,4.,0.,6.,
                0.,-1.,-3.,2.,0.,
                0.,0.,1.,0.,0.,
                0.,4.,2.,0.,1.;

            b =  8.,
                45.,
                -3.,
                3.,
                19.;

            x = 1.,
                2.,
                3.,
                4.,
                5.;

            cout << "GMRESSolver: nonsingular system\n";
            GMRESSolver& gmres = *solver;
            vector_type soln(N);
            for (index_type i = 0; i < N; ++i)
                soln(i) = rand() / (real_type(1) + RAND_MAX);
            GMRESParams params;
            params.verbose = true;
            GMRESOut result = gmres.solve(Matvec(A), Precon(), b, soln, params);
            cout << result << "\n";
            cout << soln << "\n";
            Assert::That(result.flag, Equals(ConvFlag::success));
        }

        It(Solve_singular_linear_system) {
            index_type N = 5;
            vector_type b(N);
            matrix_type A(N, N);

            A = 2.,3.,0.,0.,2.,
                3.,0.,4.,0.,3.,
                0.,-1.,-3.,2.,0.,
                0.,0.,1.,0.,0.,
                0.,4.,2.,0.,0.;

            b =  8.,
                45.,
                -3.,
                3.,
                19.;

            cout << "GMRESSolver: singular system\n";
            GMRESSolver& gmres = *solver;
            vector_type soln(N);
            for (index_type i = 0; i < N; ++i)
                soln(i) = rand() / (real_type(1) + RAND_MAX);
            GMRESOut result = gmres.solve(Matvec(A), Precon(), b, soln);
            cout << result << "\n";
            cout << soln << "\n";
            Assert::That(result.flag, Equals(ConvFlag::maxits));
        }

        void TearDown() {
            delete solver;
        }
    };
}