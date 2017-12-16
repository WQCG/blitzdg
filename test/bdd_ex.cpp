#include <igloo/igloo_alt.h>
using namespace igloo;

Describe(a_guitar_with_a_fuzzbox)
{
  void SetUp()
  {
    //load things
  }

  It(starts_in_clean_mode)
  {
    Assert::That(15, Equals(15));
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