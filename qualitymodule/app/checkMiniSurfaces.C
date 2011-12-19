#include "GoTools/compositemodel/SurfaceModel.h"
#include "GoTools/compositemodel/CompositeModelFactory.h"
#include "GoTools/qualitymodule/FaceSetQuality.h"
#include "GoTools/utils/Point.h"
#include "GoTools/compositemodel/ftPoint.h"
#include "GoTools/compositemodel/ftCurve.h"
#include <fstream>

using namespace Go;

int main( int argc, char* argv[] )
{
  if (argc != 3) {
    std::cout << "Input parameters : Input file on g2 format, mini edge size" << std::endl;
    exit(-1);
  }

  // Read input arguments
  std::ifstream file1(argv[1]);
  ALWAYS_ERROR_IF(file1.bad(), "Input file not found or file corrupt");

    double gap = 0.001;
  double neighbour = 0.01;
  double kink = 0.01;
  double approxtol = 0.01;

  CompositeModelFactory factory(approxtol, gap, neighbour, kink, 10.0*kink);

  shared_ptr<CompositeModel> model = shared_ptr<CompositeModel>(factory.createFromG2(file1));
  shared_ptr<SurfaceModel> sfmodel = dynamic_pointer_cast<SurfaceModel,CompositeModel>(model);

  FaceSetQuality quality(gap, kink, approxtol);
  quality.attach(sfmodel);

  double mini = atof(argv[2]);
  quality.setMiniElementSize(mini);

  vector<shared_ptr<ftSurface> > mini_surface;
  quality.miniSurfaces(mini_surface);

  std::cout << "Number of mini surfaces: " << mini_surface.size() << std::endl;

  std::ofstream out_file("mini_surf.g2");
  size_t ki;

  for (ki=0; ki<mini_surface.size(); ki++)
  {
      mini_surface[ki]->surface()->writeStandardHeader(out_file);
      mini_surface[ki]->surface()->write(out_file);
  }
  vector<shared_ptr<ftSurface> > mini_surface2;
  quality.miniSurfaces(mini_surface2);

  std::ofstream out_file2("mini_surf2.g2");

  for (ki=0; ki<mini_surface2.size(); ki++)
  {
      mini_surface2[ki]->surface()->writeStandardHeader(out_file2);
      mini_surface2[ki]->surface()->write(out_file2);
  }
  

}
