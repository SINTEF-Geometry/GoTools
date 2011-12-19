
#ifdef __BORLANDC__
#include <vcl.h>
#endif

#include "GoTools/compositemodel/SurfaceModel.h"
#include "GoTools/compositemodel/CompositeModelFactory.h"
#include <fstream>

//using namespace std;
using namespace Go;

#define TOPOLOGY_DEBUG

int main( int argc, char* argv[] )
{
  if (argc != 13) {
      std::cout << "Input parameters : corner (x,y,z), side 1 (x,y,z), other vec (x,y,z), l1, l2 l3" << std::endl;
    exit(-1);
  }

#ifdef __BORLANDC__
  using Go::Point;
#endif

  Point corner(atof(argv[1]), atof(argv[2]), atof(argv[3]));
  Point side1(atof(argv[4]), atof(argv[5]), atof(argv[6]));
  Point vec(atof(argv[7]), atof(argv[8]), atof(argv[9]));
  double l1 = atof(argv[10]);
  double l2 = atof(argv[11]);
  double l3 = atof(argv[12]);
		   
  double gap = 0.001;
  double neighbour = 0.01;
  double kink = 0.01;
  double approxtol = 0.01;

  CompositeModelFactory factory(approxtol, gap, neighbour, kink, 10.0*kink);

  CompositeModel *model = factory.createFromBox(corner, side1, vec, l1, l2, l3);

  SurfaceModel *sfmodel = dynamic_cast<SurfaceModel*>(model);

  if (sfmodel)
  {
      std::ofstream out_file("box.g2");
      int nmb = sfmodel->nmbEntities();
      for (int ki=0; ki<nmb; ki++)
      {
	  shared_ptr<ParamSurface> surf = sfmodel->getSurface(ki);

	  surf->writeStandardHeader(out_file);
	  surf->write(out_file);
      }
  }

  delete model;
}

	  

