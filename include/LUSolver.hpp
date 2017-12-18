#include <blitz/array.h>

using namespace blitz;

struct sparseTriplet {
    int nz;
    int *row;
    int *col;
    double *val;
};


class LUSolver {
    int N;
    Array<double, 2> * A;

    sparseTriplet triplet;

    void toSparseTriplet();

  public:
    LUSolver(Array<double, 2> * const &);
    
    Array<double, 2> & get_A();

    void factorize();

    Array<double, 2> const & solve(Array<double, 2> * const &, Array<double, 2> * &);

    ~LUSolver();
};

