#include <blitz/array.h>

using namespace std;
using namespace blitz;

class Nodes1DProvisioner {
    
    double Min_x;
    double Max_x;
    int NumElements;
    int NOrder;
    
    Array<double, 1> xGrid;
    Array<double, 1> rGrid;

    Array<double, 2> Dr;

    void computeJacobiPolynomial(Array<double,1> const & x, const double alpha, const double beta, const int N,  Array<double,1> & p);

  public:
    Nodes1DProvisioner(int NOrder, int NumElements, double xmin, double xmax);

    void buildNodes();

    void buildDr();

    Array<double, 2> & get_xGrid();
    Array<double, 1> & get_rGrid();
    Array<double, 2> & get_Dr();

    ~Nodes1DProvisioner();
};