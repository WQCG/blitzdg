#include <blitz/array.h>

using namespace std;

class MeshManager {

  public:
    MeshManager();

    // Read gmsh .msh file.
    void ReadMesh(string gmshInputFile);
    
    void ReadMesh(string vertFile, string E2VFile);

    // Split up with metis.
    void PartitionMesh(int numPartitions);
    
    ~MeshManager();
};

