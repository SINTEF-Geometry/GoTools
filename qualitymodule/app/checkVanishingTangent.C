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

  FaceSetQuality quality(gap, kink, approxtol);
  quality.attach(sfmodel);

  vector<shared_ptr<PointOnCurve> > sing_points;
  vector<pair<shared_ptr<PointOnCurve>, shared_ptr<PointOnCurve> > > sing_curves;
  quality.vanishingCurveTangent(sing_points, sing_curves);

  std::cout << "Number of singular points: " << sing_points.size() << std::endl;
  std::cout << "Number of singular curves: " << sing_curves.size() << std::endl;

  std::ofstream out_file("crv_singularity.g2");
  size_t ki;
  if (sing_points.size() > 0)
  {
      out_file << "400 1 0 4 255 0 0 255" << std::endl;
      out_file << sing_points.size() << std::endl;
      for (ki=0; ki<sing_points.size(); ki++)
      {
	  Point pnt = sing_points[ki]->getPos();
	  out_file << pnt[0] << "  ";
	  out_file << pnt[1] << "  ";
	  out_file << pnt[2] << "  " << std::endl;
	  out_file << std::endl;
      }
  }
      
  for (ki=0; ki<sing_curves.size(); ki++)
  {
      shared_ptr<ParamCurve> crv = 
	  shared_ptr<ParamCurve>(sing_curves[ki].first->getCurve()->subCurve(sing_curves[ki].first->getPar(),
								 sing_curves[ki].second->getPar()));
      crv->writeStandardHeader(out_file);
      crv->write(out_file);
  }
      
  vector<shared_ptr<PointOnCurve> > sing_points2;
  vector<pair<shared_ptr<PointOnCurve>, shared_ptr<PointOnCurve> > > sing_curves2;
  quality.vanishingCurveTangent(sing_points2, sing_curves2);

  std::ofstream out_file2("crv_singularity2.g2");
  if (sing_points.size() > 0)
  {
      out_file2 << "400 1 0 4 255 0 0 255" << std::endl;
      out_file2 << sing_points2.size() << std::endl;
      for (ki=0; ki<sing_points2.size(); ki++)
      {
	  Point pnt = sing_points2[ki]->getPos();
	  out_file2 << pnt[0] << "  ";
	  out_file2 << pnt[1] << "  ";
	  out_file2 << pnt[2] << "  " << std::endl;
	  out_file2 << std::endl;
      }
  }
      
  for (ki=0; ki<sing_curves2.size(); ki++)
  {
      shared_ptr<ParamCurve> crv = 
	  shared_ptr<ParamCurve>(sing_curves2[ki].first->getCurve()->subCurve(sing_curves2[ki].first->getPar(),
								 sing_curves2[ki].second->getPar()));
      crv->writeStandardHeader(out_file2);
      crv->write(out_file2);
  }
      


}
