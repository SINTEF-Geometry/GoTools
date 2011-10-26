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
  vector<shared_ptr<ParamSurface> > deg_sfs;
  quality.degenSurfaces(deg_sfs);

  std::cout << "Number of degenerate surfaces: " << deg_sfs.size() << std::endl;

  std::ofstream out_file("deg_sfs.g2");
  for (size_t ki=0; ki<deg_sfs.size(); ki++)
  {
      deg_sfs[ki]->writeStandardHeader(out_file);
      deg_sfs[ki]->write(out_file);
  }

  vector<shared_ptr<ParamSurface> > deg_sfs2;
  quality.degenSurfaces(deg_sfs2);

  std::ofstream out_file2("deg_sfs2.g2");
  for (size_t ki=0; ki<deg_sfs2.size(); ki++)
  {
      deg_sfs2[ki]->writeStandardHeader(out_file2);
      deg_sfs2[ki]->write(out_file2);
  }

}
