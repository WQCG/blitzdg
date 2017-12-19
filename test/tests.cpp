#include <igloo/igloo_alt.h>
#include <blitz/array.h>
#include <LUSolver.hpp>

using namespace igloo;
using namespace blitz;

const int N=5;
const double eps=2.e-15;

Array<double,2> A(N,N), B(N,N), C(N,N), D(N,N);
Array<double,1> b(N), soln(N), d(N), e(N), x(N);

firstIndex ii;
secondIndex jj;

Describe(Simple_blitz_array_operations)
{
  void SetUp() {

    A = 2,3,0,0,0,
		    3,0,4,0,6,
		    0,-1,-3,2,0,
		    0,0,1,0,0,
		    0,4,2,0,1;

	  B = jj;

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

    d = 1,2,3,4;
    e = 2,3,4,5;
  }

  It(Properly_Multiplies_Pointwise)
  {
    Array<double, 1> result(N);
    result = d(ii)*e(ii);

    Assert::That(result(0), Equals(2.));
    Assert::That(result(1), Equals(6.));
    Assert::That(result(2), Equals(12.));
    Assert::That(result(3), Equals(20.));
  }

  It(Properly_Does_Dot_Product)
  {
    double result;
    result = sum(d * e);

    Assert::That(result, Equals(40.));
  }
  It(Properly_Does_Matrix_Vector_Product)
  {
    Array <double, 1> result(N);
    result = sum(A(ii,jj)*x(jj), jj);

    Assert::That(result(0), Equals(b(0)));
    Assert::That(result(1), Equals(b(1)));
    Assert::That(result(2), Equals(b(2)));
    Assert::That(result(3), Equals(b(3)));
    Assert::That(result(4), Equals(b(4)));
  }
};

Describe(LUSolver_Object)
{
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

  }
  It(Solves_Ax_equals_b) 
  {

    LUSolver luSolver(&A);
    Array <double, 1> soln(N);

	  // Compute LU factors.
	  luSolver.factorize();
    luSolver.solve(b, soln);
    
    Assert::That(abs(soln(0)-x(0)), IsLessThan(eps));
    Assert::That(abs(soln(1)-x(1)), IsLessThan(eps));
    Assert::That(abs(soln(2)-x(2)), IsLessThan(eps));
    Assert::That(abs(soln(3)-x(3)), IsLessThan(eps));
    Assert::That(abs(soln(4)-x(4)), IsLessThan(eps));
  }
};

int main(const int argc, const char *argv[])
{
  return TestRunner::RunAllTests(argc, argv);
  // exit code returns number of failed tests.
}
