#ifdef __BORLANDC__
#include <vcl.h>
#endif

#include "GoTools/compositemodel/SurfaceModel.h"
#include "GoTools/compositemodel/CompositeCurve.h"
#include "GoTools/compositemodel/CompositeModelFactory.h"
#include <fstream>

//using namespace std;
using std::shared_ptr;
using namespace Go;

int main( int argc, char* argv[] )
{
  if (argc != 2) {
    std::cout << "Input parameters : Input file on IGES format," << std::endl;
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

  CompositeModel *model = factory.createFromIges(file1);

  SurfaceModel *sfmodel = dynamic_cast<SurfaceModel*>(model);

  CompositeCurve *cvmodel = dynamic_cast<CompositeCurve*>(model);

  if (sfmodel)
  {
      std::ofstream out_file("iges.g2");
      int nmb = sfmodel->nmbEntities();
      for (int ki=0; ki<nmb; ki++)
      {
	  shared_ptr<ParamSurface> surf = sfmodel->getSurface(ki);

	  surf->writeStandardHeader(out_file);
	  surf->write(out_file);
      }
  }
  if (cvmodel)
    {
      std::ofstream out_file("iges.g2");
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
