// Copyright (C) 2017-2018  Derek Steinmoeller. 
// See COPYING and LICENSE files at project root for more details. 

#include <igloo/igloo_alt.h>

int main(const int argc, const char *argv[])
{
  using namespace igloo;
  return TestRunner::RunAllTests(argc, argv);
  // exit code returns number of failed tests.
}
