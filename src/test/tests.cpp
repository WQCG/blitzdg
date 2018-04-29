// Copyright (C) 2017-2018  Derek Steinmoeller. 
// See COPYING and LICENSE files at project root for more details. 

#include <igloo/igloo_alt.h>
#include <blitz/array.h>
using namespace igloo;
using namespace blitz;
using namespace std;

const int N=5;
const double eps=10*numeric_limits<double>::epsilon();
const float epsf = 5.e-7;
firstIndex ii;
secondIndex jj;


Array<double, 1> b(N), x(N), d(N), e(N);
Array<double, 2> A(N,N), B(N,N);


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

int main(const int argc, const char *argv[])
{
  return TestRunner::RunAllTests(argc, argv);
  // exit code returns number of failed tests.
}
