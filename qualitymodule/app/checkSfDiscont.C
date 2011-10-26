#include "GoTools/compositemodel/SurfaceModel.h"
#include "GoTools/compositemodel/CompositeModelFactory.h"
#include "GoTools/qualitymodule/FaceSetQuality.h"
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

    double gap = 0.001;
  double neighbour = 0.01;
  double kink = 0.01;
  double approxtol = 0.01;

  CompositeModelFactory factory(approxtol, gap, neighbour, kink, 10.0*kink);

  shared_ptr<CompositeModel> model = shared_ptr<CompositeModel>(factory.createFromG2(file1));
  shared_ptr<SurfaceModel> sfmodel = dynamic_pointer_cast<SurfaceModel,CompositeModel>(model);

  FaceSetQuality quality(gap, kink, approxtol);
  quality.attach(sfmodel);

  vector<shared_ptr<ftSurface> >  g1discont_faces;
  vector<shared_ptr<ftSurface> >  c1discont_faces;
  quality.sfG1Discontinuity(g1discont_faces);
  quality.sfC1Discontinuity(c1discont_faces);

  std::cout << "Number of geometric discontinuous faces: " << g1discont_faces.size() << std::endl;
  std::cout << "Number of parametric discontinuous faces: " << c1discont_faces.size() << std::endl;

  std::ofstream out_file("g1discont_sfs.g2");
  size_t ki;
  for (ki=0; ki<g1discont_faces.size(); ki++)
  {
      shared_ptr<ParamSurface> surf1 = g1discont_faces[ki]->surface();
      surf1->writeStandardHeader(out_file);
      surf1->write(out_file);
  }

  std::ofstream out_file2("c1discont_sfs.g2");
  for (ki=0; ki<c1discont_faces.size(); ki++)
  {
      shared_ptr<ParamSurface> surf1 = c1discont_faces[ki]->surface();
      surf1->writeStandardHeader(out_file2);
      surf1->write(out_file2);
  }

  vector<shared_ptr<ftSurface> >  g1discont_faces2;
  vector<shared_ptr<ftSurface> >  c1discont_faces2;
  quality.sfG1Discontinuity(g1discont_faces2);
  quality.sfC1Discontinuity(c1discont_faces2);

  std::ofstream out_filen("g1discont_sfs2.g2");
  for (ki=0; ki<g1discont_faces2.size(); ki++)
  {
      shared_ptr<ParamSurface> surf1 = g1discont_faces2[ki]->surface();
      surf1->writeStandardHeader(out_filen);
      surf1->write(out_filen);
  }

  std::ofstream out_filen2("c1discont_sfs2.g2");
  for (ki=0; ki<c1discont_faces2.size(); ki++)
  {
      shared_ptr<ParamSurface> surf1 = c1discont_faces2[ki]->surface();
      surf1->writeStandardHeader(out_filen2);
      surf1->write(out_filen2);
  }
}
