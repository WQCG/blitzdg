#include <blitz/array.h>
#include <boost/algorithm/string.hpp>

using namespace std;
using namespace boost;

class MeshManager {
    double * Vert;
    int * EToV;
    int Dim;
    int NumVerts;
    int ElementType;
    int NumElements;
    string CsvDelimeters;

    template<typename T>
    vector<int> readCsvFile(string csvFile, string delimiters, T * & result);

    template<typename T>
    void printArray(T * & arr, int numRows, int numCols);

  public:
    MeshManager();

    int get_Index(int row, int col, int numCols);

    // Read gmsh .msh file.
    void readMesh(string gmshInputFile);
    
    void readVertices(string vertFile);

    void readElements(string E2VFile);

    // Split up with metis.
    void partitionMesh(int numPartitions);

    double * & get_Vertices();

    int get_Dim();

    int get_NumVerts();

    int * & get_Elements();

    int get_NumElements();
	  int get_ElementType();

    void printVertices();
    void printElements();
    
    ~MeshManager();
};

