#include "GoTools/compositemodel/SurfaceModel.h"
#include "GoTools/compositemodel/CompositeModelFactory.h"
#include "GoTools/qualitymodule/FaceSetQuality.h"
#include "GoTools/compositemodel/PointOnEdge.h"
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

    double gap = 0.0001;
  double neighbour = 0.01;
  double kink = 0.001;
  double approxtol = 0.01;

  CompositeModelFactory factory(approxtol, gap, neighbour, kink, 10.0*kink);

  shared_ptr<CompositeModel> model = shared_ptr<CompositeModel>(factory.createFromG2(file1));
  shared_ptr<SurfaceModel> sfmodel = dynamic_pointer_cast<SurfaceModel,CompositeModel>(model);

  FaceSetQuality quality(gap, kink, approxtol);
  quality.attach(sfmodel);

  vector<pair<shared_ptr<PointOnEdge>, shared_ptr<PointOnEdge> > > intersections;
  quality.loopSelfIntersection(intersections);

  std::cout << "Number of self intersections: " <<  intersections.size() << std::endl;
  std::ofstream out_file("loop_self_ints.g2");
  size_t ki;
  for (ki=0; ki<intersections.size(); ++ki)
  {
      shared_ptr<ParamCurve> crv1 = intersections[ki].first->edge()->geomCurve();
      crv1->writeStandardHeader(out_file);
      crv1->write(out_file);
      const Point& pos1 = intersections[ki].first->position();
      out_file << "400 1 0 4 255 0 0 255" << std::endl;
      out_file << "1 " << std::endl;
      out_file << pos1[0] << "  " << pos1[1] << "  " << pos1[2] << std::endl;

      crv1 = intersections[ki].second->edge()->geomCurve();
      crv1->writeStandardHeader(out_file);
      crv1->write(out_file);
      const Point& pos2 = intersections[ki].second->position();
      out_file << "400 1 0 4 0 255 0 255" << std::endl;
      out_file << "1 " << std::endl;
      out_file << pos2[0] << "  " << pos2[1] << "  " << pos2[2] << std::endl;
  }

  vector<pair<shared_ptr<PointOnEdge>, shared_ptr<PointOnEdge> > > intersections2;
  quality.loopSelfIntersection(intersections2);

  std::ofstream out_file2("loop_self_ints2.g2");
  for (ki=0; ki<intersections2.size(); ++ki)
  {
      shared_ptr<ParamCurve> crv1 = intersections2[ki].first->edge()->geomCurve();
      crv1->writeStandardHeader(out_file2);
      crv1->write(out_file2);
      const Point& pos1 = intersections2[ki].first->position();
      out_file2 << "400 1 0 4 255 0 0 255" << std::endl;
      out_file2 << "1 " << std::endl;
      out_file2 << pos1[0] << "  " << pos1[1] << "  " << pos1[2] << std::endl;

      crv1 = intersections2[ki].second->edge()->geomCurve();
      crv1->writeStandardHeader(out_file2);
      crv1->write(out_file2);
      const Point& pos2 = intersections2[ki].second->position();
      out_file2 << "400 1 0 4 0 255 0 255" << std::endl;
      out_file2 << "1 " << std::endl;
      out_file2 << pos2[0] << "  " << pos2[1] << "  " << pos2[2] << std::endl;
  }

}
