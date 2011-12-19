#ifdef __BORLANDC__
#include <vcl.h>
#endif

#include "GoTools/compositemodel/CompositeCurve.h"
#include "GoTools/compositemodel/CompositeModelFactory.h"
#include <fstream>

//using namespace std;
using namespace Go;

int main( int argc, char* argv[] )
{
  if (argc != 7) {
      std::cout << "Input parameters : start point (x,y,z), end point (x,y,z)" << std::endl;
    exit(-1);
  }

#ifdef __BORLANDC__
  using Go::Point;
#endif

  Point startpt(atof(argv[1]), atof(argv[2]), atof(argv[3]));
  Point endpt(atof(argv[4]), atof(argv[5]), atof(argv[6]));
  double gap = 0.001;
  double neighbour = 0.01;
  double kink = 0.01;
  double approxtol = 0.01;

  CompositeModelFactory factory(approxtol, gap, neighbour, kink, 10.0*kink);

  CompositeModel *model = factory.createLineSegment(startpt, endpt);

  CompositeCurve *cvmodel = dynamic_cast<CompositeCurve*>(model);

  if (cvmodel)
  {
      std::ofstream out_file("line.g2");
      int nmb = cvmodel->nmbEntities();
      for (int ki=0; ki<nmb; ki++)
      {
	  shared_ptr<ParamCurve> crv = cvmodel->getCurve(ki);

	  crv->writeStandardHeader(out_file);
	  crv->write(out_file);
      }
  }

  delete model;
}
