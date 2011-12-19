#ifdef __BORLANDC__
#include <vcl.h>
#endif

#include "GoTools/compositemodel/SurfaceModel.h"
#include "GoTools/compositemodel/ftSurface.h"
#include "GoTools/compositemodel/CompositeModelFactory.h"
#include <fstream>

//using namespace std;
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

  CompositeModel *model = factory.createFromG2(file1);

  SurfaceModel *sfmodel = dynamic_cast<SurfaceModel*>(model);

  if (sfmodel)
  {
      int ki;
      std::ofstream out_file("fileg2.g2");
      int nmb = sfmodel->nmbEntities();
      vector<shared_ptr<ftSurface> > faces;
      for (ki=0; ki<nmb; ki++)
      {
	  shared_ptr<ParamSurface> surf = sfmodel->getSurface(ki);

	  surf->writeStandardHeader(out_file);
	  surf->write(out_file);

	  shared_ptr<ftSurface> curr_face = sfmodel->getFace(ki);
	  faces.push_back(curr_face);
      }

      shared_ptr<SurfaceModel> model2 = 
	shared_ptr<SurfaceModel>(new SurfaceModel(approxtol, gap,
						  neighbour, kink,
						  10.0*kink, faces, true));

      std::ofstream out_file2("bd.g2");
      int nmb_bd = sfmodel->nmbBoundaries();
      std::cout << "Number of boundaries: " << nmb_bd << std::endl;
      for (ki=0; ki<nmb_bd; ki++)
      {
	  ftCurve bdcv = sfmodel->getBoundary(ki);
	  bdcv.writeSpaceCurve(out_file2);
      }

  }

  delete model;
}
