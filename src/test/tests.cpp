// Copyright (C) 2017-2018  Waterloo Quantitative Consulting Group, Inc.
// See COPYING and LICENSE files at project root for more details.

#include <igloo/igloo_alt.h>
#include <Warning.hpp>

using namespace igloo;
using namespace std;
using namespace blitzdg;

int main(const int argc, const char *argv[])
{
  printDisclaimer();
  return TestRunner::RunAllTests(argc, argv);
  // exit code returns number of failed tests.
}
