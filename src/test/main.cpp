// Copyright (C) 2017-2019  Waterloo Quantitative Consulting Group, Inc.
// See COPYING and LICENSE files at project root for more details.

#include <igloo/igloo_alt.h>
#include "Warning.hpp"

int main(const int argc, char **argv)
{
  using namespace igloo;
  blitzdg::printDisclaimer();
  return TestRunner::RunAllTests(argc, argv);
  // exit code returns number of failed tests.
}
