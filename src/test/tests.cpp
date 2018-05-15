// Copyright (C) 2017-2018  Waterloo Quantitative Consulting Group, Inc.
// See COPYING and LICENSE files at project root for more details.

#include <igloo/igloo_alt.h>
using namespace igloo;
using namespace std;

void printDisclaimer() {
	cout << "blitzdg, version 0.1.0a" << endl;
	cout << "Copyright (C) 2017-2018 Waterloo Quantitative Consulting Group, Inc." << endl;
	cout << "This is free software; see the source code for copying conditions." << endl;
	cout << "There is ABSOLUTELY NO WARRANTY; not even for MERCHANTABILITY or" << endl;
	cout << "FITNESS FOR A PARTICULAR PURPOSE." << endl << endl;
}

int main(const int argc, const char *argv[])
{
  printDisclaimer();
  return TestRunner::RunAllTests(argc, argv);
  // exit code returns number of failed tests.
}
