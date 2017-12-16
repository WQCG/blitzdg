#include <blitz/array.h>

using namespace blitz;

class LUFactorizer {
    int N;
    const Array<double, 2> * A;
  public:
    LUFactorizer(Array<double, 2> * const &);
    
    Array<double, 2> const & get_A();

    void factorize();

    Array<double, 2> const & solve(const Array<double, 2>);

    ~LUFactorizer();
};

