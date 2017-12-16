#include <iostream>
#include <blitz/array.h>
#include <LUFactorizer.hpp>

using namespace std;
using namespace blitz;

int main(int argc, char **argv) {
	cout << "hello, world\n" ;

	int N = 8;

    // Create two-dimensional arrays of float
    Array<double,2> A(N,N), B(N,N), C(N,N), D(N,N);
	Array<double, 2> * const & ptr = &A;

	firstIndex ii;
	secondIndex jj;

	A = ii;
	B = jj;

    C = A*B;
	D = A+B;

	cout << "C:" << C << endl ;
	cout << "D:" << D << endl ;

	cout << "ptr:" << ptr << endl;
	cout << "ptr (dereferenced):" << *ptr << endl;

	LUFactorizer factorizer(ptr);

	factorizer.factorize();

	B = factorizer.get_A();

	cout << B << endl;

    return 0;

}