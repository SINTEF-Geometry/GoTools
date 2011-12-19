#include "GoTools/geometry/ParamSurface.h"
#include "GoTools/compositemodel/SurfaceModel.h"
#include "GoTools/compositemodel/CompositeModelFactory.h"
#include "GoTools/compositemodel/ftSurface.h"
#include "GoTools/intersections/Identity.h"
#include <fstream>

//using namespace std;
using namespace Go;

int main( int argc, char* argv[] )
{
  if (argc != 4) {
      std::cout << "Surface file1(g2), surface file2, tolerance" << std::endl;
    exit(-1);
  }

  std::ifstream file1(argv[1]);
  ALWAYS_ERROR_IF(file1.bad(), "Input file not found or file corrupt");

  std::ifstream file2(argv[2]);
  ALWAYS_ERROR_IF(file1.bad(), "Input file not found or file corrupt");

  double tol = atof(argv[3]);

  double gap = 0.001;
  double neighbour = 0.01;
  double kink = 0.01;
  double approxtol = 0.01;

  CompositeModelFactory factory(approxtol, gap, neighbour, kink, 10.0*kink);
  CompositeModel *model1 = factory.createFromG2(file1);
  CompositeModel *model2 = factory.createFromG2(file2);
  SurfaceModel *sfmodel1 = dynamic_cast<SurfaceModel*>(model1);
  SurfaceModel *sfmodel2 = dynamic_cast<SurfaceModel*>(model2);

  int nmb1 = sfmodel1->nmbEntities();
  int nmb2 = sfmodel2->nmbEntities();
  int ki, kj;

  Identity ident;
  for (ki=0; ki<nmb1; ++ki)
    {
      shared_ptr<ftSurface> face1 = sfmodel1->getFace(ki);
      for (kj=0; kj<nmb2; ++kj)
	{
	  shared_ptr<ftSurface> face2 = sfmodel2->getFace(kj);

	  // Check coincidence
	  
	  int coincide = ident.identicalSfs(face1->surface(),
					    face2->surface(),
					    tol);

	  if (coincide)
	    {
	      bool connected = face1->connectTwin(face2.get(), tol);
	      
	      std::cout << "Connected: " << connected << std::endl;
	    }
	}
    }
}
