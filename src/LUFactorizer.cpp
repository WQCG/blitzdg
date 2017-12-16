#include <suitesparse/umfpack.h>
#include <LUFactorizer.hpp>
#include <iostream>

using namespace std;

LUFactorizer::LUFactorizer(Array<double, 2> * const & Ain) {
    A = Ain;
}

void LUFactorizer::factorize() {
    cout << "Computing factorization!" << endl;
}

Array<double, 2> const & LUFactorizer::get_A() {
    return *A;
}


LUFactorizer::~LUFactorizer() {

}