#include "GoTools/compositemodel/Body.h"
#include "GoTools/compositemodel/SurfaceModel.h"
#include "GoTools/compositemodel/CompositeModelFactory.h"
#include <fstream>

//using namespace std;
using namespace Go;

int main( int argc, char* argv[] )
{
  if (argc != 5) {
      std::cout << "Input parameters : centre (x,y,z), radius" << std::endl;
    exit(-1);
  }

  Point centre(atof(argv[1]), atof(argv[2]), atof(argv[3]));
  double radius = atof(argv[4]);
    double gap = 0.001;
  double neighbour = 0.01;
  double kink = 0.01;
  double approxtol = 0.01;

  CompositeModelFactory factory(approxtol, gap, neighbour, kink, 10.0*kink);

  shared_ptr<CompositeModel> model = (shared_ptr<CompositeModel>)(factory.createFromSphere(centre, radius));

  Body body(dynamic_pointer_cast<SurfaceModel,CompositeModel>(model));

  shared_ptr<SurfaceModel> sfmodel = body.getOuterShell();

  if (sfmodel.get())
  {
      std::ofstream out_file("sphere2.g2");
      int nmb = sfmodel->nmbEntities();
      for (int ki=0; ki<nmb; ki++)
      {
	  shared_ptr<ParamSurface> surf = sfmodel->getSurface(ki);

	  surf->writeStandardHeader(out_file);
	  surf->write(out_file);
      }
  }

}

