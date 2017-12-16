#include <iostream>
#include <blitz/array.h>
#include <rectangle.hpp>

using namespace std;
using namespace blitz;

int main(int argc, char **argv) {
	cout << "hello, world\n" ;

	int N = 8;

    // Create two-dimensional arrays of float
    Array<float,2> A(N,N), B(N,N), C(N,N), D(N,N);

	firstIndex ii;
	secondIndex jj;

	A = ii;
	B = jj;

    C = A*B;
	D = A+B;

	cout << "C:" << C << endl ;
	cout << "D:" << D << endl ;

	Rectangle rect;
  	rect.set_values (3,4);
  	cout << "area: " << rect.area();


    return 0;

}