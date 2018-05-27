// Copyright (C) 2017-2018  Waterloo Quantitative Consulting Group, Inc.
// See COPYING and LICENSE files at project root for more details.

#include "EigenSolver.hpp"
#include "JacobiBuilders.hpp"
#include "Types.hpp"
#include <igloo/igloo_alt.h>

using std::cout;
using std::endl;
using std::numeric_limits;
using std::abs;

namespace blitzdg {

    namespace JacobiBuildersTests {
        using namespace igloo;

        const real_type eps = 10*numeric_limits<real_type>::epsilon();
        const float epsf = 5.7e-6;

        JacobiBuilders * Jacobi = nullptr;


        Describe(JacobiBuilders_Object) {
            void SetUp() {
                Jacobi = new JacobiBuilders();        
            }

            void TearDown() {
                delete Jacobi;
            }

            It(Should_Generate_0th_Order_Legendre_Polynomial) {
                cout << "Should_Generate_0th_Order_Legendre_Polynomial" << endl;
                real_vector_type x(3);
                x = -1.,0.,1.;
                real_vector_type  p(3);

                const JacobiBuilders & jacobi = *Jacobi;

                jacobi.computeJacobiPolynomial(x, 0.0, 0.0, 0, p);

                Assert::That(p(0), Equals(1./sqrt(2.)));
                Assert::That(p(1), Equals(1./sqrt(2.)));
                Assert::That(p(2), Equals(1./sqrt(2.)));
            }

            It(Should_Generate_1st_Order_Legendre_Polynomial) {
                cout << "Should_Generate_1st_Order_Legendre_Polynomial" << endl;
                real_vector_type x(3);
                x = -1.,0.,1.;
                real_vector_type  p(3);

                JacobiBuilders & jacobi = *Jacobi;

                jacobi.computeJacobiPolynomial(x, 0.0, 0.0, 1, p);

                Assert::That(p(0), Equals(-sqrt(3./2.)));
                Assert::That(p(1), Equals(0.0));
                Assert::That(p(2), Equals(sqrt(3./2.)));
            }


            It(Should_Generate_2nd_Order_Legendre_Polynomial) {
                cout << "Should_Generate_2nd_Order_Legendre_Polynomial" << endl;
                real_vector_type x(3);
                x = -1.,0.,1.;
                real_vector_type  p(3);

                JacobiBuilders & jacobi = *Jacobi;

                jacobi.computeJacobiPolynomial(x, 0.0, 0.0, 2, p);

                Assert::That(abs(p(0)-sqrt(5./2.)), IsLessThan(eps));
                Assert::That(abs(p(1)-(-sqrt(5./8.))), IsLessThan(eps));
                Assert::That(abs(p(2)-sqrt(5./2.)), IsLessThan(eps));
            }

            It(Should_Generate_0th_Order_Legendre_Polynomial_4pt_Grid) {
                cout << "Should_Generate_0th_Order_Legendre_Polynomial_4pt_Grid" << endl;
                real_vector_type x(4);
                x = -1,-0.447214,0.447214,1;
                real_vector_type p(4);

                JacobiBuilders & jacobi = *Jacobi;

                jacobi.computeJacobiPolynomial(x, 0.0, 0.0, 0, p);

                Assert::That(abs(p(0)-sqrt(1./2.)), IsLessThan(eps));
                Assert::That(abs(p(1)-sqrt(1./2.)), IsLessThan(eps));
                Assert::That(abs(p(2)-sqrt(1./2.)), IsLessThan(eps));
            } 

            It(Should_Generate_1st_Order_Legendre_Polynomial_4pt_Grid) {
                cout << "Should_Generate_1st_Order_Legendre_Polynomial_4pt_Grid" << endl;
                real_vector_type x(4);
                x = -1,-0.447214,0.447214,1;
                real_vector_type p(4);

                JacobiBuilders & jacobi = *Jacobi;

                jacobi.computeJacobiPolynomial(x, 0.0, 0.0, 1, p);

                Assert::That(abs(p(0)- -1.224744871391589), IsLessThan(epsf));
                Assert::That(abs(p(1)- -0.547722557505166), IsLessThan(epsf));
                Assert::That(abs(p(2)-  0.547722557505166), IsLessThan(epsf));
                Assert::That(abs(p(3)-  1.224744871391589), IsLessThan(epsf));
            } 

            It(Should_Generate_4th_Order_Quadrature_Points_and_Weights) {
                cout << "Should_Generate_4th_Order_Quadrature_Points_and_Weights" << endl;
                
                JacobiBuilders & jacobi = *Jacobi;
                
                real_vector_type x(5);
                real_vector_type w(5);

                jacobi.computeJacobiQuadWeights(0., 0., 4, x, w);

                Assert::That(abs(x(0) - -9.06179845938664e-01), IsLessThan(eps));
                Assert::That(abs(x(1) - -5.38469310105683e-01), IsLessThan(eps));
                Assert::That(abs(x(2) - -9.62591786604533e-17), IsLessThan(eps));
                Assert::That(abs(x(3) -  5.38469310105683e-01), IsLessThan(eps));
                Assert::That(abs(x(4) -  9.06179845938664e-01), IsLessThan(eps));

                Assert::That(abs(w(0) - 0.236926885056189), IsLessThan(eps));
                Assert::That(abs(w(1) - 0.478628670499366), IsLessThan(eps));
                Assert::That(abs(w(2) - 0.568888888888889), IsLessThan(eps));
                Assert::That(abs(w(3) - 0.478628670499367), IsLessThan(eps));
                Assert::That(abs(w(4) - 0.236926885056189), IsLessThan(eps));
            }

            It(Should_Generate_1st_Order_Quadrature_Points_and_Weights) {
                cout << "Should_Generate_1st_Order_Quadrature_Points_and_Weights" << endl;
                
                JacobiBuilders & jacobi = *Jacobi;
                
                real_vector_type x(2);
                real_vector_type w(2);

                jacobi.computeJacobiQuadWeights(0., 0., 1, x, w);

                Assert::That(abs(x(0) - -0.577350269189626), IsLessThan(eps));
                Assert::That(abs(x(1) -  0.577350269189626), IsLessThan(eps));

                Assert::That(abs(w(0) - 1.0), IsLessThan(eps));
                Assert::That(abs(w(1) - 1.0), IsLessThan(eps));
            }

            It(Should_Compute_2nd_Order_Quadrature_Points_and_Weights) {
                cout << "Should_Compute_2nd_Order_Quadrature_Points_and_Weights" << endl;

                JacobiBuilders & jacobi = *Jacobi;

                real_vector_type xx(3);
                real_vector_type ww(3);

                jacobi.computeJacobiQuadWeights(0., 0., 2, xx, ww);

                Assert::That(abs(xx(0) - -7.74596669241483e-01), IsLessThan(eps));
                Assert::That(abs(xx(1) -  0.0), IsLessThan(eps));
                Assert::That(abs(xx(2) -  7.74596669241483e-01), IsLessThan(eps));

                Assert::That(abs(ww(0) -  0.555555555555556), IsLessThan(eps));
                Assert::That(abs(ww(1) -  0.888888888888889), IsLessThan(eps));
                Assert::That(abs(ww(2) -  0.555555555555556), IsLessThan(eps));
            }

            It(Should_Generate_3rd_Order_Legendre_Gauss_Lobatto_Nodes) {
                cout << "Should_Generate_3rd_Order_Legendre_Gauss_Lobatto_Nodes" << endl;
                JacobiBuilders & jacobi = *Jacobi;
                
                real_vector_type x(4);
                jacobi.computeGaussLobottoPoints(0., 0., 3, x);

                Assert::That(abs(x(0) - -1), IsLessThan(eps));
                Assert::That(abs(x(1) - -0.447213595499958), IsLessThan(eps));
                Assert::That(abs(x(2) -  0.447213595499958), IsLessThan(eps));
                Assert::That(abs(x(3) -  1), IsLessThan(eps));
            }
        };
    }
}