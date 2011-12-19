#include "GoTools/compositemodel/SurfaceModel.h"
#include "GoTools/compositemodel/ftCurve.h"
#include "GoTools/compositemodel/CompositeModelFactory.h"
#include <fstream>

//using namespace std;
using namespace Go;

int main( int argc, char* argv[] )
{
  if (argc != 2) {
    std::cout << "Input parameters : Input file on g2 format" << std::endl;
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

  std::ofstream of1("boundary.g2");

  std::vector<shared_ptr<SurfaceModel> > models2 =
    sfmodel->getConnectedModels();

  std::cout << "Number of connected models: " << models2.size() << std::endl;

  
  // Outer boundaries
  int ki, kj;
  int nmb_bd = sfmodel->nmbBoundaries();
  std::cout << "Outer boundaries: " << nmb_bd << std::endl;
  for (ki=0; ki<nmb_bd; ++ki)
    {
      ftCurve bd = sfmodel->getBoundary(ki);
      bd.writeSpaceCurve(of1);
    }

  for (ki=0; ki<(int)models2.size(); ++ki)
    {
      std::ofstream of2("connected.g2");
      int nmb = models2[ki]->nmbEntities();
      for (kj=0; kj<nmb; kj++)
      {
	  shared_ptr<ParamSurface> surf = models2[ki]->getSurface(kj);

	  surf->writeStandardHeader(of2);
	  surf->write(of2);
      }
    }
}
