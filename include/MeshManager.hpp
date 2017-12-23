#include <blitz/array.h>

using namespace std;

class MeshManager {
    double * Vert;
    int Dim;
    int NumVerts;

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
    
    ~MeshManager();
};

