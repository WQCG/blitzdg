#include <blitz/array.h>
#include <SparseTriplet.hpp>
#include <suitesparse/umfpack.h>

using namespace blitz;

class SparseMatrixConverter {

  public:
    SparseMatrixConverter();

    void fullToSparseTriplet(const Array<double, 2> & A, SparseTriplet & triplet);

    ~SparseMatrixConverter();
};
