// Copyright (C) 2017-2018  Derek Steinmoeller. 
// See COPYING and LICENSE files at project root for more details. 

#include <igloo/igloo_alt.h>
#include <blitz/array.h>
#include <Nodes1DProvisioner.hpp>
#include <LUSolver.hpp>
#include <EigenSolver.hpp>
#include <DirectSolver.hpp>

using namespace igloo;
using namespace blitz;
using namespace std;

namespace Nodes1DProvisionerTests {
    const int N=5;
    const double eps=10*numeric_limits<double>::epsilon();
    const float epsf = 5.7e-6;

    firstIndex ii;
    secondIndex jj;

    LUSolver * luSolver = nullptr;
    SparseMatrixConverter * matrixConverter = nullptr;
    Nodes1DProvisioner * nodes1DProvisioner = nullptr;
    EigenSolver * eigenSolver = nullptr;
    DirectSolver * directSolver = nullptr;

    Describe(Nodes1DProvisioner_Object) {
        void SetUp() {
            const int NOrder = 3;
            const int NumElements = 5;
            const double xmin = -1.0;
            const double xmax = 1.0;

            matrixConverter = new SparseMatrixConverter();
            eigenSolver = new EigenSolver(*matrixConverter);
            directSolver = new DirectSolver(*matrixConverter);
            nodes1DProvisioner = new Nodes1DProvisioner(NOrder, NumElements, xmin, xmax, *matrixConverter, *eigenSolver, *directSolver);
        }

        It(Should_Generate_0th_Order_Legendre_Polynomial) {
            cout << "Nodes1D" << endl;
            Array<double, 1> x(3);
            x = -1.,0.,1.;
            Array<double, 1>  p(3);

            Nodes1DProvisioner & nodes1D = *nodes1DProvisioner;

            nodes1D.computeJacobiPolynomial(x, 0.0, 0.0, 0, p);

            Assert::That(p(0), Equals(1./sqrt(2.)));
            Assert::That(p(1), Equals(1./sqrt(2.)));
            Assert::That(p(2), Equals(1./sqrt(2.)));
        }

        It(Should_Generate_1st_Order_Legendre_Polynomial) {
            cout << "Should_Generate_1st_Order_Legendre_Polynomial" << endl;
            Array<double, 1> x(3);
            x = -1.,0.,1.;
            Array<double, 1>  p(3);

            Nodes1DProvisioner & nodes1D = *nodes1DProvisioner;

            nodes1D.computeJacobiPolynomial(x, 0.0, 0.0, 1, p);

            Assert::That(p(0), Equals(-sqrt(3./2.)));
            Assert::That(p(1), Equals(0.0));
            Assert::That(p(2), Equals(sqrt(3./2.)));
        }


        It(Should_Generate_2nd_Order_Legendre_Polynomial) {
            cout << "Should_Generate_2nd_Order_Legendre_Polynomial" << endl;
            Array<double, 1> x(3);
            x = -1.,0.,1.;
            Array<double, 1>  p(3);

            Nodes1DProvisioner & nodes1D = *nodes1DProvisioner;

            nodes1D.computeJacobiPolynomial(x, 0.0, 0.0, 2, p);

            Assert::That(abs(p(0)-sqrt(5./2.)), IsLessThan(eps));
            Assert::That(abs(p(1)-(-sqrt(5./8.))), IsLessThan(eps));
            Assert::That(abs(p(2)-sqrt(5./2.)), IsLessThan(eps));
        }

        It(Should_Generate_0th_Order_Legendre_Polynomial_4pt_Grid) {
            cout << "Should_Generate_0th_Order_Legendre_Polynomial_4pt_Grid" << endl;
            Array<double, 1> x(4);
            x = -1,-0.447214,0.447214,1;
            Array<double, 1> p(4);

            Nodes1DProvisioner & nodes1D = *nodes1DProvisioner;

            nodes1D.computeJacobiPolynomial(x, 0.0, 0.0, 0, p);

            Assert::That(abs(p(0)-sqrt(1./2.)), IsLessThan(eps));
            Assert::That(abs(p(1)-sqrt(1./2.)), IsLessThan(eps));
            Assert::That(abs(p(2)-sqrt(1./2.)), IsLessThan(eps));
        } 

        It(Should_Generate_1st_Order_Legendre_Polynomial_4pt_Grid) {
            cout << "Should_Generate_1st_Order_Legendre_Polynomial_4pt_Grid" << endl;
            Array<double, 1> x(4);
            x = -1,-0.447214,0.447214,1;
            Array<double, 1> p(4);

            Nodes1DProvisioner & nodes1D = *nodes1DProvisioner;

            nodes1D.computeJacobiPolynomial(x, 0.0, 0.0, 1, p);

            Assert::That(abs(p(0)- -1.224744871391589), IsLessThan(epsf));
            Assert::That(abs(p(1)- -0.547722557505166), IsLessThan(epsf));
            Assert::That(abs(p(2)-  0.547722557505166), IsLessThan(epsf));
            Assert::That(abs(p(3)-  1.224744871391589), IsLessThan(epsf));
        } 

        It(Should_Generate_4th_Order_Quadrature_Points_and_Weights) {
            cout << "Should_Generate_4th_Order_Quadrature_Points_and_Weights" << endl;
            
            Nodes1DProvisioner & nodes1D = *nodes1DProvisioner;
            
            Array<double, 1> x(5);
            Array<double, 1> w(5);

            nodes1D.computeJacobiQuadWeights(0., 0., 4, x, w);

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
            
            Nodes1DProvisioner & nodes1D = *nodes1DProvisioner;
            
            Array<double, 1> x(2);
            Array<double, 1> w(2);

            nodes1D.computeJacobiQuadWeights(0., 0., 1, x, w);

            Assert::That(abs(x(0) - -0.577350269189626), IsLessThan(eps));
            Assert::That(abs(x(1) -  0.577350269189626), IsLessThan(eps));

            Assert::That(abs(w(0) - 1.0), IsLessThan(eps));
            Assert::That(abs(w(1) - 1.0), IsLessThan(eps));
        }

        It(Should_Compute_2nd_Order_Quadrature_Points_and_Weights) {
            cout << "Should_Compute_2nd_Order_Quadrature_Points_and_Weights" << endl;

            Nodes1DProvisioner & nodes1D = *nodes1DProvisioner;

            Array<double, 1> xx(3);
            Array<double, 1> ww(3);

            nodes1D.computeJacobiQuadWeights(0., 0., 2, xx, ww);

            Assert::That(abs(xx(0) - -7.74596669241483e-01), IsLessThan(eps));
            Assert::That(abs(xx(1) -  0.0), IsLessThan(eps));
            Assert::That(abs(xx(2) -  7.74596669241483e-01), IsLessThan(eps));

            Assert::That(abs(ww(0) -  0.555555555555556), IsLessThan(eps));
            Assert::That(abs(ww(1) -  0.888888888888889), IsLessThan(eps));
            Assert::That(abs(ww(2) -  0.555555555555556), IsLessThan(eps));
        }

        It(Should_Generate_3rd_Order_Legendre_Gauss_Lobatto_Nodes) {
            cout << "Should_Generate_3rd_Order_Legendre_Gauss_Lobatto_Nodes" << endl;
            Nodes1DProvisioner & nodes1D = *nodes1DProvisioner;
            
            Array<double, 1> x(4);
            nodes1D.computeGaussLobottoPoints(0., 0., 3, x);

            Assert::That(abs(x(0) - -1), IsLessThan(eps));
            Assert::That(abs(x(1) - -0.447213595499958), IsLessThan(eps));
            Assert::That(abs(x(2) -  0.447213595499958), IsLessThan(eps));
            Assert::That(abs(x(3) -  1), IsLessThan(eps));
        }

        It(Should_Build_3rd_Order_Vandermonde_Matrix) {
            cout << "Should_Build_3rd_Order_Vandermonde_Matrix" << endl;
            Nodes1DProvisioner & nodes1D = *nodes1DProvisioner;

            nodes1D.buildNodes();
            nodes1D.buildVandermondeMatrix();
            Array<double, 2> & V = nodes1D.get_V();

            Array<double, 2> expectedV(4,4);
            expectedV = 0.70711,-1.22474, 1.58114,-1.87083,
                        0.70711,-0.54772,-0.31623, 0.83666,
                        0.70711, 0.54772,-0.31623,-0.83666,
                        0.70711, 1.22474, 1.58114, 1.87083;

            Array<double, 2> res(4,4);
            res  = V - expectedV;
            Assert::That(sqrt(sum(res(ii)*res(ii))), IsLessThan(epsf));
        }

        It(Should_Build_4th_Order_GradVandermonde_Matrix) {
            cout << "Should_Build_4th_Order_GradVandermonde_Matrix" << endl;
            Nodes1DProvisioner & nodes1D = *nodes1DProvisioner;

            nodes1D.buildNodes();

            Array<double, 2> DVr(5,5);
            nodes1D.computeGradVandermonde(DVr);

            Array<double, 2> expectedDVr(5,5);
            expectedDVr = 0.00000,1.22474,-4.74342,11.22497,-21.21320,
                          0.00000,1.22474,-3.10530, 3.20713, -0.00000,
                          0.00000,1.22474,-0.00000,-2.80624,  0.00000,
                          0.00000,1.22474, 3.10530, 3.20713,  0.00000,
                          0.00000,1.22474 ,4.74342,11.22497, 21.21320;

            Array<double, 2> res(4,4);
            res = DVr - expectedDVr;
            Assert::That(sqrt(sum(res(ii)*res(ii))), IsLessThan(epsf));
        }

        It(Should_Build_3rd_Order_Differentiation_Matrix) {
            cout << "Should_Build_3rd_Order_Differentiation_Matrix" << endl;
            Nodes1DProvisioner & nodes1D = *nodes1DProvisioner;
            nodes1D.buildNodes();

            Array<double, 2> & Dr = nodes1D.get_Dr();

            Array<double, 2> expectedDr(4,4);
            expectedDr = -3.0000e+00, 4.0451e+00,-1.5451e+00, 5.0000e-01,
                         -8.0902e-01,-4.0540e-16, 1.1180e+00,-3.0902e-01,
                          3.0902e-01,-1.1180e+00, 6.2804e-16, 8.0902e-01,
                         -5.0000e-01, 1.5451e+00,-4.0451e+00, 3.0000e+00; 


            Array<double, 2> res(4,4);
            res = Dr - expectedDr;
            Assert::That(sqrt(sum(res(ii)*res(ii))), IsLessThan(epsf));
        }

        It(Should_Build_A_1D_X_Grid) {
            cout << "Should_Build_A_1D_X_Grid" << endl;

            Nodes1DProvisioner & nodes1D = *nodes1DProvisioner;

            nodes1D.buildNodes();

            Array<double, 2> & x = nodes1D.get_xGrid();

            Array<double,2> expectedx(4,5);
            expectedx = -1.000000,-0.600000,-0.200000,0.200000,0.600000,
                        -0.889443,-0.489443,-0.089443,0.310557,0.710557,
                        -0.710557,-0.310557, 0.089443,0.489443,0.889443,
                        -0.600000,-0.200000, 0.200000,0.600000,1.000000;

            Array<double, 2> res(4,5);
            res = x - expectedx;
            Assert::That(sqrt(sum(res(ii)*res(ii))), IsLessThan(epsf));
        }

        It(Should_Build_Element_To_Vertex_Connectivity) {
            cout << "Should_Build_Element_To_Vertex_Connectivity" << endl;
            Nodes1DProvisioner & nodes1D = *nodes1DProvisioner;

            nodes1D.buildNodes();

            Array<int, 2> EToV = nodes1D.get_EToV();

            Assert::That(EToV(0,0), Equals(0)); Assert::That(EToV(0,1), Equals(1));
            Assert::That(EToV(1,0), Equals(1)); Assert::That(EToV(1,1), Equals(2));
            Assert::That(EToV(2,0), Equals(2)); Assert::That(EToV(2,1), Equals(3));
            Assert::That(EToV(3,0), Equals(3)); Assert::That(EToV(3,1), Equals(4));
            Assert::That(EToV(4,0), Equals(4)); Assert::That(EToV(4,1), Equals(5));
        }

        It(Should_Compute_Jacobian) {
            cout << "Should_Compute_Jacobian" << endl;
            Nodes1DProvisioner & nodes1D = *nodes1DProvisioner;

            nodes1D.buildNodes();

            nodes1D.buildDr();
            nodes1D.computeJacobian();

            Array<double, 2> & J = nodes1D.get_J();
            Array<double, 2> & rx = nodes1D.get_rx();

            Array<double,2> expectedJ(4,5), expectedrx(4,5);
            expectedJ = 0.20000,0.20000,0.20000,0.20000,0.20000,
                        0.20000,0.20000,0.20000,0.20000,0.20000,
                        0.20000,0.20000,0.20000,0.20000,0.20000,
                        0.20000,0.20000,0.20000,0.20000,0.20000;

            expectedrx =5,5,5,5,5,
                        5,5,5,5,5,
                        5,5,5,5,5,
                        5,5,5,5,5;

            Array<double,2> resJ(4,5), resrx(4,5);
            resJ = J - expectedJ;
            resrx = rx - expectedrx;

            Assert::That(sqrt(sum(resJ(ii)*resJ(ii))), IsLessThan(epsf));
            Assert::That(sqrt(sum(resrx(ii)*resrx(ii))), IsLessThan(epsf));
        }

        It(Should_Build_1D_Lift_Operator) {
            cout << "Should_Build_1D_Lift_Operator" << endl;
            Nodes1DProvisioner & nodes1D = *nodes1DProvisioner;

            nodes1D.buildNodes();
            Array<double, 2> Lift = nodes1D.get_Lift();

            Array<double, 2> expectedLift(4,2);
            expectedLift =  8.00000,-2.00000,
                           -0.89443, 0.89443,
                            0.89443,-0.89443,
                           -2.00000, 8.00000;


            Array<double, 2> resLift(4,2);

            resLift = Lift - expectedLift;
            Assert::That(sqrt(sum(resLift*resLift)), IsLessThan(epsf));
        }

        It(Should_Build_1D_Connectivity_Matrices) {
            cout << "Should_Build_1D_Connectivity_Matrices" << endl;
            Nodes1DProvisioner & nodes1D = *nodes1DProvisioner;

            nodes1D.buildNodes();
            Array<int, 2> EToE = nodes1D.get_EToE();
            Array<int, 2> EToF = nodes1D.get_EToF();

            Array<int, 2> expectedEToE(5,2);
            expectedEToE = 0,1,
                           0,2,
                           1,3,
                           2,4,
                           3,4;

            Array<int, 2> expectedEToF(5,2);
            expectedEToF = 0,0,
                           1,0,
                           1,0,
                           1,0,
                           1,1;

            Array<int, 2> resEToE(5,2), resEToF(5,2);

            resEToE = EToE - expectedEToE;
            resEToF = EToF - expectedEToF;

            Assert::That(sqrt(sum(resEToE*resEToE)), Equals(0));
            Assert::That(sqrt(sum(resEToF*resEToF)), Equals(0));
        }
    };
}