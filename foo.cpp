#include <iostream>
#include <blitz/array.h>

using namespace std;
using namespace blitz;

int main(int argc, char **argv) {
	cout << "hello, world\n" ;

	int N = 8;

    // Create two-dimensional arrays of float
    Array<float,2> A(N,N), B(N,N), C(N,N), D(N,N);


	firstIndex ii;
	secondIndex jj;

	//A = 10.0;
	//B = 20.0;

	A = ii;
	B = jj;

    C = A*B;
	D = A+B;

	cout << "C:" << C << endl ;
	cout << "D:" << D << endl ;

    return 0;

}