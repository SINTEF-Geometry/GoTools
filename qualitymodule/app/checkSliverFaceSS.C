#include "GoTools/compositemodel/SurfaceModel.h"
#include "GoTools/compositemodel/CompositeModelFactory.h"
#include "GoTools/qualitymodule/FaceSetQuality.h"
#include <fstream>

using std::vector;
using namespace Go;

int main( int argc, char* argv[] )
{
  if (argc != 3 && argc != 4) {
    std::cout << "Input parameters : Input file on g2 format, thickness for thin side, [minimal proportion thick/thin side]" << std::endl;
    exit(-1);
  }

  // Read input arguments
  std::ifstream file1(argv[1]);
  ALWAYS_ERROR_IF(file1.bad(), "Input file not found or file corrupt");

  double thickness = atof(argv[2]);

    double gap = 0.01;
  double neighbour = 0.1;
  double kink = 0.01;
  double approxtol = 0.01;

  CompositeModelFactory factory(approxtol, gap, neighbour, kink, 10.0*kink);

  shared_ptr<CompositeModel> model = shared_ptr<CompositeModel>(factory.createFromG2(file1));
  shared_ptr<SurfaceModel> sfmodel = dynamic_pointer_cast<SurfaceModel,CompositeModel>(model);

  FaceSetQuality quality(gap, kink, approxtol);
  quality.attach(sfmodel);

  vector<shared_ptr<ParamSurface> > sliv_sfs;
  if (argc == 3)
    quality.sliverSurfaces(sliv_sfs, thickness);
  else
    {
      double factor = atof(argv[3]);
      quality.sliverSurfaces(sliv_sfs, thickness, factor);
    }
  std::cout << "Number of sliver faces: " << sliv_sfs.size() << std::endl;

  std::ofstream out_file("sliv_sfs.g2");
  for (size_t ki=0; ki<sliv_sfs.size(); ki++)
  {
      sliv_sfs[ki]->writeStandardHeader(out_file);
      sliv_sfs[ki]->write(out_file);
  }

  vector<shared_ptr<ParamSurface> > sliv_sfs2;
  if (argc == 3)
    quality.sliverSurfaces(sliv_sfs2, thickness);
  else
    {
      double factor = atof(argv[3]);
      quality.sliverSurfaces(sliv_sfs2, thickness, factor);
    }

  std::ofstream out_file2("sliv_sfs2.g2");
  for (size_t ki=0; ki<sliv_sfs2.size(); ki++)
  {
      sliv_sfs2[ki]->writeStandardHeader(out_file2);
      sliv_sfs2[ki]->write(out_file2);
  }

}
