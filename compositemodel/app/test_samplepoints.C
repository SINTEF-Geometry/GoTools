#ifdef __BORLANDC__
#include <vcl.h>
#endif

#include "GoTools/compositemodel/SurfaceModel.h"
#include "GoTools/compositemodel/CompositeModelFactory.h"
#include "GoTools/compositemodel/FaceUtilities.h"
#include "GoTools/utils/Point.h"
#include <fstream>
#include <stdlib.h> // For atof()

using namespace std;
using namespace Go;

int main( int argc, char* argv[] )
{
  if (argc != 4) {
    std::cout << "Input parameters : Input file, IGES or g2 (1/0), density"  << std::endl;
    exit(-1);
  }

  // Read input arguments
  std::ifstream file1(argv[1]);
  ALWAYS_ERROR_IF(file1.bad(), "Input file not found or file corrupt");

  double gap = 0.001;
  double neighbour = 0.01;
  double kink = 0.01;
  // Removed from corresponding use of fairingToolbox
  double approx = 0.001;
  int useIGES = atoi(argv[2]);
  double density = atof(argv[3]);

  CompositeModelFactory factory(approx, gap, neighbour, kink, 10.0*kink);

  CompositeModel *model;
  if (useIGES)
      model = factory.createFromIges(file1);
  else
      model = factory.createFromG2(file1);

  SurfaceModel *sfmodel = dynamic_cast<SurfaceModel*>(model);

  if (sfmodel)
    {
      std::vector<SamplePointData> points;
      sfmodel->fetchSamplePoints(density, points);

      size_t ki;
      int nmb_norm=0;
      std::ofstream out1("sample_points.g2");
      out1 << "400 1 0 4 255 0 0 255" << std::endl;
      out1 << points.size() << std::endl;
      for (ki=0; ki<points.size(); ++ki)
	{
	  out1 << points[ki].pos_ << std::endl;
	  if (points[ki].norm_[0] < MAXDOUBLE)
	    nmb_norm++;
	}
      out1 << "410 1 0 4 0  255 0 255" << std::endl;
      out1 << nmb_norm << std::endl;
      for (ki=0; ki<points.size(); ++ki)
	{
	  if (points[ki].norm_[0] < MAXDOUBLE)
	    {
	      out1 << points[ki].pos_ << " ";
	      out1 << points[ki].pos_+points[ki].norm_ << std::endl;
	    }
	}
       int break_point;
      break_point = 1;
    }

  delete model;
}

