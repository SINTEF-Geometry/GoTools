#include "GoTools/compositemodel/SurfaceModel.h"
#include "GoTools/compositemodel/ftSurface.h"
#include "GoTools/compositemodel/CompositeModelFactory.h"
#include "GoTools/compositemodel/RegularizeFace.h"
#include <fstream>

using std::vector;
using std::shared_ptr;
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

    double gap = 0.001;
  double neighbour = 0.01;
  double kink = 0.01;
  double approxtol = 0.01;

  CompositeModelFactory factory(approxtol, gap, neighbour, kink, 10.0*kink);

  CompositeModel *model = factory.createFromG2(file1);

  SurfaceModel *sfmodel = dynamic_cast<SurfaceModel*>(model);

   if (sfmodel)
  {
    shared_ptr<ftSurface> face = sfmodel->getFace(0);

    RegularizeFace reg(face, gap, kink, neighbour);
    vector<shared_ptr<ftSurface> > sub_faces = reg.getRegularFaces();

    std::ofstream of("regularized_faces.g2");
    for (size_t ki=0; ki<sub_faces.size(); ++ki)
      {
	shared_ptr<ParamSurface> surf = sub_faces[ki]->surface();
	surf->writeStandardHeader(of);
	surf->write(of);
      }

    // Replace by spline surfaces
    shared_ptr<SurfaceModel> model2 =
      shared_ptr<SurfaceModel>(new SurfaceModel(approxtol, gap, neighbour,
						kink, 10.0*kink, sub_faces,
						true));

    model2->replaceRegularSurfaces();
    int nmb = model2->nmbEntities();
    for (int kr=0; kr<nmb; ++kr)
      {
	shared_ptr<ParamSurface> sf = model2->getSurface(kr);
	sf->writeStandardHeader(file2);
	sf->write(file2);
      }
  }
}


