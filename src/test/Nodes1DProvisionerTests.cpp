#include <igloo/igloo_alt.h>
#include <blitz/array.h>
#include <Nodes1DProvisioner.hpp>
#include <LUSolver.hpp>
#include <EigenSolver.hpp>
#include <MeshManager.hpp>

using namespace igloo;
using namespace blitz;
using namespace std;

namespace nodes1DProvisionerTests {
    const int N=5;
    const double eps=10*numeric_limits<double>::epsilon();
    const float epsf = 5.e-7;

    firstIndex ii;
    secondIndex jj;

    LUSolver * luSolver = nullptr;
    MeshManager * meshManager=nullptr;
    SparseMatrixConverter * matrixConverter = nullptr;
    Nodes1DProvisioner * nodes1DProvisioner = nullptr;
    EigenSolver * eigenSolver = nullptr;

    Describe(Nodes1DProvisioner_Object) {
        void SetUp() {
            const int NOrder = 3;
            const int NumElements = 5;
            const double xmin = -1.0;
            const double xmax = 1.0;

            matrixConverter = new SparseMatrixConverter();
            eigenSolver = new EigenSolver(*matrixConverter);

            nodes1DProvisioner = new Nodes1DProvisioner(NOrder, NumElements, xmin, xmax, *matrixConverter, *eigenSolver);
        }

        It(Should_Generate_0th_Order_Legendre_Polynomial) {
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
            
            Nodes1DProvisioner & nodes1D = *nodes1DProvisioner;
            
            Array<double, 1> x(2);
            Array<double, 1> w(2);

            nodes1D.computeJacobiQuadWeights(0., 0., 1, x, w);

            Assert::That(abs(x(0) - -0.577350269189626), IsLessThan(eps));
            Assert::That(abs(x(1) -  0.577350269189626), IsLessThan(eps));

            Assert::That(abs(w(0) - 1.0), IsLessThan(eps));
            Assert::That(abs(w(1) - 1.0), IsLessThan(eps));
        }

        It(Should_Generate_3rd_Order_Legendre_Gauss_Lobatto_Nodes) {
            Nodes1DProvisioner & nodes1D = *nodes1DProvisioner;
            
            Array<double, 1> x(4);
            nodes1D.computeGaussLobottoPoints(0., 0., 3, x);

            Assert::That(abs(x(0) - -1), IsLessThan(eps));
            Assert::That(abs(x(1) - -0.447213595499958), IsLessThan(eps));
            Assert::That(abs(x(2) -  0.447213595499958), IsLessThan(eps));
            Assert::That(abs(x(3) -  1), IsLessThan(eps));
        }

        It(Should_Build_3rd_Order_Vandermonde_Matrix) {
            Nodes1DProvisioner & nodes1D = *nodes1DProvisioner;

            nodes1D.buildNodes();
            nodes1D.buildVandermondeMatrix();
            Array<double, 2> & V = nodes1D.get_V();

            Array<double, 2> expectedV(4,4);
            expectedV = 0.70711,-1.22474, 1.58114,-1.87083,
                        0.70711,-0.54772,-0.31623,0.83666,
                        0.70711,0.54772,-0.31623,-0.83666,
                        0.70711,1.22474,1.58114,1.87083;

            Array<double, 2> res(4,4);
            res  = V - expectedV;
            Assert::That(sum(res(ii)*res(ii)), IsLessThan(epsf));
        }

        It(Should_Build_4th_Order_GradVandermonde_Matrix) {
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
            Assert::That(sum(res(ii)*res(ii)), IsLessThan(epsf));
        }
    };
}