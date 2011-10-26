#include "GoTools/compositemodel/SurfaceModel.h"
#include "GoTools/compositemodel/CompositeModelFactory.h"
#include "GoTools/qualitymodule/FaceSetQuality.h"
#include "GoTools/geometry/PointOnCurve.h"
#include <fstream>

using std::vector;
using std::shared_ptr;
using std::dynamic_pointer_cast;
using namespace Go;

int main( int argc, char* argv[] )
{
  if (argc != 3) {
    std::cout << "Input parameters : Input file on g2 format, minium curvature radius" << std::endl;
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

  double curvature_radius = atof(argv[2]);
  quality.setCurvatureRadius(curvature_radius);

  pair<shared_ptr<PointOnCurve>,double> min_curv_rad;
  vector<pair<shared_ptr<PointOnCurve>, double> > small_curv_rad;
  quality.cvCurvatureRadius(small_curv_rad, min_curv_rad);

  std::cout << "Number of small radiuses: " << small_curv_rad.size() << std::endl;
  std::cout << "Minimum curvature radius: " << min_curv_rad.second << std::endl;

  std::ofstream out_file("curvature_cvs.g2");
  for (size_t ki=0; ki<small_curv_rad.size(); ki++)
  {
      out_file << "400 1 0 4 255 0 0 255" << std::endl;
      out_file << "1 " << std::endl;
      Point pos = small_curv_rad[ki].first->getPos();
      out_file << pos[0] << "  ";
      out_file << pos[1] << "  ";
      out_file << pos[2] << "  " << std::endl;
      out_file << std::endl;
      small_curv_rad[ki].first->getCurve()->writeStandardHeader(out_file);
      small_curv_rad[ki].first->getCurve()->write(out_file);
      std::cout << "Curvature radius: " << small_curv_rad[ki].second << std::endl;
  }
  out_file << "400 1 0 4 0 255 0 255" << std::endl;
  out_file << "1 " << std::endl;
  Point pos2 = min_curv_rad.first->getPos();
  out_file << pos2[0] << "  ";
  out_file << pos2[1] << "  ";
  out_file << pos2[2] << "  " << std::endl;
  out_file << std::endl;
  min_curv_rad.first->getCurve()->writeStandardHeader(out_file);
  min_curv_rad.first->getCurve()->write(out_file);

  pair<shared_ptr<PointOnCurve>,double> min_curv_rad2;
  vector<pair<shared_ptr<PointOnCurve>, double> > small_curv_rad2;
  quality.cvCurvatureRadius(small_curv_rad2, min_curv_rad2);

  std::cout << "Number of small radiuses: " << small_curv_rad2.size() << std::endl;
  std::cout << "Minimum curvature radius: " << min_curv_rad2.second << std::endl;

  std::ofstream out_file2("curvature_cvs2.g2");
  for (size_t ki=0; ki<small_curv_rad2.size(); ki++)
  {
      out_file2 << "400 1 0 4 255 0 0 255" << std::endl;
      out_file2 << "1 " << std::endl;
      Point pos = small_curv_rad2[ki].first->getPos();
      out_file2 << pos[0] << "  ";
      out_file2 << pos[1] << "  ";
      out_file2 << pos[2] << "  " << std::endl;
      out_file2 << std::endl;
      small_curv_rad2[ki].first->getCurve()->writeStandardHeader(out_file2);
      small_curv_rad2[ki].first->getCurve()->write(out_file2);
      std::cout << "Curvature radius: " << small_curv_rad2[ki].second << std::endl;
  }
  out_file2 << "400 1 0 4 0 255 0 255" << std::endl;
  out_file2 << "1 " << std::endl;
  pos2 = min_curv_rad2.first->getPos();
  out_file2 << pos2[0] << "  ";
  out_file2 << pos2[1] << "  ";
  out_file2 << pos2[2] << "  " << std::endl;
  out_file2 << std::endl;
  min_curv_rad2.first->getCurve()->writeStandardHeader(out_file2);
  min_curv_rad2.first->getCurve()->write(out_file2);


}
