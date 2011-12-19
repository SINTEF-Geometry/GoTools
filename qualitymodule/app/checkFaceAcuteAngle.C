#include "GoTools/compositemodel/SurfaceModel.h"
#include "GoTools/compositemodel/CompositeModelFactory.h"
#include "GoTools/qualitymodule/FaceSetQuality.h"
#include <fstream>

using std::vector;
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
  double kink = 0.1;
  double approxtol = 0.01;

  CompositeModelFactory factory(approxtol, gap, neighbour, kink, 10.0*kink);

  shared_ptr<CompositeModel> model = shared_ptr<CompositeModel>(factory.createFromG2(file1));
  shared_ptr<SurfaceModel> sfmodel = dynamic_pointer_cast<SurfaceModel,CompositeModel>(model);

  FaceSetQuality quality(gap, kink, approxtol);
  quality.attach(sfmodel);

  vector<pair<ftSurface*, ftSurface*> >  faces;
  quality.acuteFaceAngle(faces);

  std::cout << "Number of acute face angles: " << faces.size() << std::endl;

  std::ofstream out_file("acute_angle_sfs.g2");
  size_t ki;
  for (ki=0; ki<faces.size(); ki++)
  {
      shared_ptr<ParamSurface> surf1 = faces[ki].first->surface();
      surf1->writeStandardHeader(out_file);
      surf1->write(out_file);
      shared_ptr<ParamSurface> surf2 = faces[ki].second->surface();
      surf2->writeStandardHeader(out_file);
      surf2->write(out_file);
  }

  vector<pair<ftSurface*, ftSurface*> >  faces2;
  quality.acuteFaceAngle(faces2);

  std::ofstream out_file2("acute_angle_sfs2.g2");
  for (ki=0; ki<faces2.size(); ki++)
  {
      shared_ptr<ParamSurface> surf1 = faces2[ki].first->surface();
      surf1->writeStandardHeader(out_file2);
      surf1->write(out_file2);
      shared_ptr<ParamSurface> surf2 = faces2[ki].second->surface();
      surf2->writeStandardHeader(out_file2);
      surf2->write(out_file2);
  }

}
