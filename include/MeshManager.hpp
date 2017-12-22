#include <blitz/array.h>

using namespace std;

class MeshManager {
    double * Vert;
    int Dim;
    int NumVerts;

  public:
    MeshManager();

    int GetIndex(int, int, int);

    // Read gmsh .msh file.
    void ReadMesh(string gmshInputFile);
    
    void ReadVertices(string vertFile);

    void ReadElementToVertex(string E2VFile);

    // Split up with metis.
    void PartitionMesh(int numPartitions);

    double * & GetVertices();

    int GetDim();

    int GetNumVerts();
    
    ~MeshManager();
};

