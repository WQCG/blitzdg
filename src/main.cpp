#include <iostream>
#include <blitz/array.h>
#include <LUSolver.hpp>

using namespace std;
using namespace blitz;

int main(int argc, char **argv) {
	cout << "hello, world\n" ;

	int N = 5;

    // Create two-dimensional arrays of float
    Array<double,2> A(N,N), B(N,N), C(N,N), D(N,N);

	firstIndex ii;
	secondIndex jj;

	A = 2,3,0,0,0,
		3,0,4,0,6,
		0,-1,-3,2,0,
		0,0,1,0,0,
		0,4,2,0,1;
	B = jj;

    C = A*B;
	D = A+B;

	cout << "C:" << C << endl ;
	cout << "D:" << D << endl ;

	LUSolver luSolver(&A);

	// Compute LU factors.
	luSolver.factorize();

	// B gets a reference to A (not a copy)
	B = luSolver.get_A();

	//B and A will output the same values
	cout << B << endl;
	cout << A << endl;

	// But their memory addresses will be different.
	cout << &B << endl << &A << endl;

    return 0;
}