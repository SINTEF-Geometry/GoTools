#include "GoTools/compositemodel/SurfaceModel.h"
#include "GoTools/compositemodel/ftSurface.h"
#include "GoTools/compositemodel/CompositeModelFactory.h"
#include "GoTools/compositemodel/RegularizeFaceSet.h"
#include <fstream>

//using namespace std;
using namespace Go;

int main( int argc, char* argv[] )
{
  if (argc != 3) {
    std::cout << "Input parameters : Input file(g2), output file" << std::endl;
    exit(-1);
  }

  // Read input arguments
  std::ifstream file1(argv[1]);
  ALWAYS_ERROR_IF(file1.bad(), "Input file not found or file corrupt");

  std::ofstream file2(argv[2]);

  double gap = 0.001; // 0.001;
  double neighbour = 0.01; // 0.01;
  double kink = 0.01;
  double approxtol = 0.01;

  CompositeModelFactory factory(approxtol, gap, neighbour, kink, 10.0*kink);

  CompositeModel *model = factory.createFromG2(file1);

  SurfaceModel *sfmodel = dynamic_cast<SurfaceModel*>(model);

   if (sfmodel)
  {
    std::vector<shared_ptr<ftSurface> > faces = sfmodel->allFaces();

    //RegularizeFaceSet reg(faces, gap, kink, true);
    RegularizeFaceSet reg(faces, gap, kink, false);
    std::vector<shared_ptr<ftSurface> > sub_faces = reg.getRegularFaces();

    for (size_t ki=0; ki<sub_faces.size(); ++ki)
      {
	shared_ptr<ParamSurface> surf = sub_faces[ki]->surface();
	surf->writeStandardHeader(file2);
	surf->write(file2);
      }
  }
}

