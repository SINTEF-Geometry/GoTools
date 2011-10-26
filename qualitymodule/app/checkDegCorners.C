#include "GoTools/compositemodel/SurfaceModel.h"
#include "GoTools/compositemodel/CompositeModelFactory.h"
#include "GoTools/qualitymodule/FaceSetQuality.h"
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

  vector<shared_ptr<ftPoint> > deg_corners;
  quality.degenerateSfCorners(deg_corners);

  std::cout << "Number of degenerate corners: " << deg_corners.size() << std::endl;

  std::ofstream out_file("deg_corners.g2");
  for (size_t ki=0; ki<deg_corners.size(); ki++)
  {
      out_file << "400 1 0 4 255 0 0 255" << std::endl;
      out_file << "1 " << std::endl;
      const Point& pos = deg_corners[ki]->position();
      out_file << pos[0] << "  ";
      out_file << pos[1] << "  ";
      out_file << pos[2] << "  " << std::endl;
      out_file << std::endl;
      deg_corners[ki]->face()->surface()->writeStandardHeader(out_file);
      deg_corners[ki]->face()->surface()->write(out_file);
  }
  vector<shared_ptr<ftPoint> > deg_corners2;
  quality.degenerateSfCorners(deg_corners2);

  std::ofstream out_file2("deg_corners2.g2");
  for (size_t ki=0; ki<deg_corners2.size(); ki++)
  {
      out_file2 << "400 1 0 4 255 0 0 255" << std::endl;
      out_file2 << "1 " << std::endl;
      const Point& pos = deg_corners2[ki]->position();
      out_file2 << pos[0] << "  ";
      out_file2 << pos[1] << "  ";
      out_file2 << pos[2] << "  " << std::endl;
      out_file2 << std::endl;
      deg_corners2[ki]->face()->surface()->writeStandardHeader(out_file2);
      deg_corners2[ki]->face()->surface()->write(out_file2);
  }

}
