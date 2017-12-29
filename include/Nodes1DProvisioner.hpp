#include <blitz/array.h>
#include <boost/algorithm/string.hpp>

using namespace std;
using namespace boost;

class Nodes1DProvisioner {
    
    double Min_x;
    double Max_x;
    int NumElements;
    int NOrder;
    
    Array<double, 1> xGrid;
    Array<double, 1> rGrid;

    Array<double, 2> Dr;

  public:
    Nodes1DProvisioner(int NOrder, int NumElements, double xmin, double xmax);

    void buildNodes();

    void buildDr();

    Array<double, 1> & get_xGrid();
    Array<double, 1> & get_rGrid();
    Array<double, 2> & get_Dr();

    ~Nodes1DProvisioner();
}