#include "GoTools/compositemodel/SurfaceModel.h"
#include "GoTools/compositemodel/CompositeModelFactory.h"
#include "GoTools/qualitymodule/FaceSetQuality.h"
#include "GoTools/compositemodel/PointOnEdge.h"
#include <fstream>

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

    double gap = 0.0001;
  double neighbour = 0.1;
  double kink = 0.001;
  double approxtol = 0.01;

  CompositeModelFactory factory(approxtol, gap, neighbour, kink, 10.0*kink);

  shared_ptr<CompositeModel> model = shared_ptr<CompositeModel>(factory.createFromG2(file1));
  shared_ptr<SurfaceModel> sfmodel = dynamic_pointer_cast<SurfaceModel,CompositeModel>(model);

  FaceSetQuality quality(sfmodel->getTolerances(), approxtol);
  quality.attach(sfmodel);

  vector<pair<shared_ptr<PointOnEdge>, shared_ptr<PointOnEdge> > > narrow;
  quality.narrowRegion(narrow);

  std::cout << "Number of narrow regions: " <<  narrow.size() << std::endl;
  std::ofstream out_file("narrow.g2");
  size_t ki;
  for (ki=0; ki<narrow.size(); ++ki)
  {
      shared_ptr<ParamCurve> crv1 = narrow[ki].first->edge()->geomCurve();
      crv1->writeStandardHeader(out_file);
      crv1->write(out_file);
      const Point& pos1 = narrow[ki].first->position();
      out_file << "400 1 0 4 255 0 0 255" << std::endl;
      out_file << "1 " << std::endl;
      out_file << pos1[0] << "  " << pos1[1] << "  " << pos1[2] << std::endl;

      crv1 = narrow[ki].second->edge()->geomCurve();
      crv1->writeStandardHeader(out_file);
      crv1->write(out_file);
      const Point& pos2 = narrow[ki].second->position();
      out_file << "400 1 0 4 0 255 0 255" << std::endl;
      out_file << "1 " << std::endl;
      out_file << pos2[0] << "  " << pos2[1] << "  " << pos2[2] << std::endl;

      shared_ptr<ParamSurface> srf = narrow[ki].first->edge()->face()->surface();
      srf->writeStandardHeader(out_file);
      srf->write(out_file);
  }

  vector<pair<shared_ptr<PointOnEdge>, shared_ptr<PointOnEdge> > > narrow2;
  quality.narrowRegion(narrow2);

  std::ofstream out_file2("narrow2.g2");
  for (ki=0; ki<narrow2.size(); ++ki)
  {
      shared_ptr<ParamCurve> crv1 = narrow2[ki].first->edge()->geomCurve();
      crv1->writeStandardHeader(out_file2);
      crv1->write(out_file2);
      const Point& pos1 = narrow2[ki].first->position();
      out_file2 << "400 1 0 4 255 0 0 255" << std::endl;
      out_file2 << "1 " << std::endl;
      out_file2 << pos1[0] << "  " << pos1[1] << "  " << pos1[2] << std::endl;

      crv1 = narrow2[ki].second->edge()->geomCurve();
      crv1->writeStandardHeader(out_file2);
      crv1->write(out_file2);
      const Point& pos2 = narrow2[ki].second->position();
      out_file2 << "400 1 0 4 0 255 0 255" << std::endl;
      out_file2 << "1 " << std::endl;
      out_file2 << pos2[0] << "  " << pos2[1] << "  " << pos2[2] << std::endl;

      shared_ptr<ParamSurface> srf = narrow2[ki].first->edge()->face()->surface();
      srf->writeStandardHeader(out_file2);
      srf->write(out_file2);
  }

}
