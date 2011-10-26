#include "GoTools/compositemodel/SurfaceModel.h"
#include "GoTools/compositemodel/CompositeModelFactory.h"
#include "GoTools/qualitymodule/FaceSetQuality.h"
#include "GoTools/utils/Point.h"
#include "GoTools/compositemodel/ftPoint.h"
#include "GoTools/compositemodel/ftCurve.h"
#include <fstream>

using std::vector;
using std::shared_ptr;
using std::dynamic_pointer_cast;
using namespace Go;

int main( int argc, char* argv[] )
{
  if (argc != 2) {
    std::cout << "Input parameters : Input file on g2 format," << std::endl;
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

  vector<shared_ptr<ftPoint> > singular_points;
  vector<shared_ptr<ftCurve> > singular_curves;
  quality.vanishingSurfaceNormal(singular_points, singular_curves);

  std::cout << "Number of singular points: " << singular_points.size() << std::endl;
  std::cout << "Number of singular curves: " << singular_curves.size() << std::endl;

  std::ofstream out_file("surf_singularity.g2");
  size_t ki;
  if (singular_points.size() > 0)
  {
      out_file << "400 1 0 4 255 0 0 255" << std::endl;
      out_file << singular_points.size() << std::endl;
  }
  for (ki=0; ki<singular_points.size(); ki++)
  {
      const Point& pos = singular_points[ki]->position();
      out_file << pos[0] << "  ";
      out_file << pos[1] << "  ";
      out_file << pos[2] << "  " << std::endl;
      out_file << std::endl;
  }
  
  for (ki=0; ki<singular_curves.size(); ki++)
  {
      singular_curves[ki]->writeSpaceCurve(out_file);
  }

  vector<shared_ptr<ftPoint> > singular_points2;
  vector<shared_ptr<ftCurve> > singular_curves2;
  quality.vanishingSurfaceNormal(singular_points2, singular_curves2);

  std::ofstream out_file2("surf_singularity2.g2");
  if (singular_points2.size() > 0)
  {
      out_file2 << "400 1 0 4 255 0 0 255" << std::endl;
      out_file2 << singular_points2.size() << std::endl;
  }
  for (ki=0; ki<singular_points.size(); ki++)
  {
      const Point& pos = singular_points2[ki]->position();
      out_file2 << pos[0] << "  ";
      out_file2 << pos[1] << "  ";
      out_file2 << pos[2] << "  " << std::endl;
      out_file2 << std::endl;
  }
  
  for (ki=0; ki<singular_curves2.size(); ki++)
  {
      singular_curves2[ki]->writeSpaceCurve(out_file2);
  }

}
