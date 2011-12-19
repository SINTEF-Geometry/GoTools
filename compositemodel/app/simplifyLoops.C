#ifdef __BORLANDC__
#include <vcl.h>
#endif

#include "GoTools/compositemodel/SurfaceModel.h"
#include "GoTools/utils/Point.h"
#include "GoTools/compositemodel/CompositeModelFactory.h"

#include <iostream>
#include <fstream>
#include <stdlib.h> // For atof()
//using namespace std;
using namespace Go;




int main( int argc, char* argv[] )
{
#ifdef __BORLANDC__
  using Go::Point;
#endif

  // Test number of input arguments
  if (argc != 4)
    {
      std::cout << "Input arguments : File type (0=g2, 1=iges), Input file,";
      std::cout << " Output file" << std::endl;
      exit(-1);
    }


  double gap = 0.001;
  double neighbour = 0.01;
  double kink = 0.01;
  double bend = 0.1;
  double approx = 0.001;


  // Read input arguments
  int is_iges = atoi(argv[1]);
  std::ifstream file1(argv[2]);
  ALWAYS_ERROR_IF(file1.bad(), "Input file not found or file corrupt");

  std::ofstream out_file(argv[3]);

  CompositeModelFactory factory(approx, gap, neighbour, kink, bend);
  SurfaceModel* sm;
  if (is_iges)
    sm = (SurfaceModel*) factory.createFromIges(file1);
  else
    sm = (SurfaceModel*) factory.createFromG2(file1);

  if (!sm)
    exit(-1);

  double dist;
  bool modified = sm->simplifyTrimLoops(dist);
  std::cout << "Model modified: " << modified << std::endl;
  if (modified)
    {
      std::cout << "Error: " << dist << std::endl;

      int nmb = sm->nmbEntities();
      for (int ki=0; ki<nmb; ++ki)
	{
	  shared_ptr<ParamSurface> srf = sm->getSurface(ki);
	  srf->writeStandardHeader(out_file);
	  srf->write(out_file);
	}
    }
}

