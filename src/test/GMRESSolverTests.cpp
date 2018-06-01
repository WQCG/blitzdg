// Copyright (C) 2017-2018  Waterloo Quantitative Consulting Group, Inc.
// See COPYING and LICENSE files at project root for more details.

#include "GMRESSolver.hpp"
#include <igloo/igloo_alt.h>
#include <blitz/array.h>
#include <iostream>

using blitz::firstIndex;
using blitz::secondIndex;
using blitz::sum;
using std::cout;

namespace blitzdg {
    namespace GMRESSolverTests {
        using namespace igloo;
        class Matvec {
            real_matrix_type mat_;
        public:
            explicit Matvec(const real_matrix_type& mat)
                : mat_{ mat }
            {}

            bool operator()(const real_vector_type& in, real_vector_type& out) const {
                firstIndex i;
                secondIndex j;
                out = sum(mat_(i, j) * in(j), j);
                return true;
            }
        };

        // identity preconditioner
        class Precon {
        public:
            bool operator()(const real_vector_type& in, real_vector_type& out) const {
                out = in;
                return true;
            }
        };

        Describe(GMRESSolver_Object) {
            It(Solve_Ax_equals_b) {
                index_type N = 5;
                real_vector_type b(N), x(N);
                real_matrix_type A(N, N);
                
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
                GMRESSolver gmres;
                real_vector_type soln(N);
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
                real_vector_type b(N);
                real_matrix_type A(N, N);

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
                GMRESSolver gmres;
                real_vector_type soln(N);
                for (index_type i = 0; i < N; ++i)
                    soln(i) = rand() / (real_type(1) + RAND_MAX);
                GMRESOut result = gmres.solve(Matvec(A), Precon(), b, soln);
                cout << result << "\n";
                cout << soln << "\n";
                Assert::That(result.flag, Equals(ConvFlag::maxits));
            }
        };
    } // namespace GMRESSolverTests
} // namespace blitzdg