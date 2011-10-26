#include "GoTools/compositemodel/SurfaceModel.h"
#include "GoTools/compositemodel/CompositeModelFactory.h"
#include "GoTools/compositemodel/ftPointSet.h"
#include <fstream>
#include <stdlib.h> // For atof()

using namespace std;
using namespace Go;

int main( int argc, char* argv[] )
{
  if (argc != 4) {
    std::cout << "Input parameters : Input file, IGES or g2 (1/0), density"  << std::endl;
    exit(-1);
  }

  // Read input arguments
  std::ifstream file1(argv[1]);
  ALWAYS_ERROR_IF(file1.bad(), "Input file not found or file corrupt");

  double gap = 0.001;
  double neighbour = 0.01;
  double kink = 0.01;
  // Removed from corresponding use of fairingToolbox
  double approx = 0.001;
  int useIGES = atoi(argv[2]);
  double density = atof(argv[3]);

  CompositeModelFactory factory(approx, gap, neighbour, kink, 10.0*kink);

  CompositeModel *model;
  if (useIGES)
      model = factory.createFromIges(file1);
  else
      model = factory.createFromG2(file1);

  SurfaceModel *sfmodel = dynamic_cast<SurfaceModel*>(model);


  std::shared_ptr<ftPointSet> triang;
  triang = sfmodel->triangulate(density);

  vector<vector<int> > tri;
  vector<Vector3D> pos;

  // Get points and triangles
  triang->getPoints(pos);
  triang->getTriangles(tri);

  std::ofstream out("triangles.g2");
  for (size_t ki=0; ki<tri.size(); ki++)
    {
      out << "410 1 0 0" << std::endl;
      out << 3 << std::endl;
      Vector3D x1 = pos[tri[ki][0]];
      Vector3D x2 = pos[tri[ki][1]];
      Vector3D x3 = pos[tri[ki][2]];

      x1.write(out);
      out << " ";
      x2.write(out);
      out << std::endl;
      x2.write(out);
      out << " ";
      x3.write(out);
      out << std::endl;
      x3.write(out);
      out << " ";
      x1.write(out);
      out << std::endl;
   }

}

