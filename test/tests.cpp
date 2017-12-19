#include <igloo/igloo_alt.h>
#include <blitz/array.h>

using namespace igloo;
using namespace blitz;

const int N=5;

Array<double,2> A(N,N), B(N,N), C(N,N), D(N,N);
Array<double,1> b(N), soln(N), d(N), e(N);

firstIndex ii;
secondIndex jj;

Describe(Simple_blitz_array_operations)
{
  void SetUp()
  {

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

    d = 1,2,3,4;
    e = 2,3,4,5;
  }

  It(Properly_Multiplies_Pointwise)
  {
    Array<double, 1> result(N);
    result = d(ii)*e(ii);

    Assert::That(2., Equals(result(0)));
    Assert::That(6., Equals(result(1)));
    Assert::That(12., Equals(result(2)));
    Assert::That(20., Equals(result(3)));
  }

  It(Properly_Does_Dot_Product) 
  {
    double result;
    result = sum(d * e);

    Assert::That(40., Equals(result));
  }

  Describe(in_distorted_mode)
  {
    void SetUp()
    {
      //load other things
    }

    It(sounds_distorted)
    {
      Assert::That(1234, Equals(12345));
    }

    It(sounds_clean_when_I_switch_the_fuzzbox)
    {
      Assert::That("foo", Equals("bar"));
    }
  };
};

int main(const int argc, const char *argv[])
{
  return TestRunner::RunAllTests(argc, argv);
  // exit code returns number of failed tests.
}
