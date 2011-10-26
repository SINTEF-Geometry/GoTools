#include "GoTools/compositemodel/SurfaceModel.h"
#include "GoTools/compositemodel/CompositeModelFactory.h"
#include "GoTools/qualitymodule/FaceSetQuality.h"
#include "GoTools/qualitymodule/FaceSetRepair.h"
#include <fstream>

using std::vector;
using std::pair;
using std::shared_ptr;
using std::dynamic_pointer_cast;
using namespace Go;

int main( int argc, char* argv[] )
{
  if (argc != 3) {
    std::cout << "Input parameters : Input file on g2 format, repair?" << std::endl;
    exit(-1);
  }

  // Read input arguments
  std::ifstream file1(argv[1]);
  ALWAYS_ERROR_IF(file1.bad(), "Input file not found or file corrupt");

  int do_repair = atoi(argv[2]);

    double gap = 0.001;
  double neighbour = 0.01;
  double kink = 0.01;
  double approxtol = 0.01;

  CompositeModelFactory factory(approxtol, gap, neighbour, kink, 10.0*kink);

  shared_ptr<CompositeModel> model = shared_ptr<CompositeModel>(factory.createFromG2(file1));
  shared_ptr<SurfaceModel> sfmodel = dynamic_pointer_cast<SurfaceModel,CompositeModel>(model);

  shared_ptr<FaceSetQuality> quality = 
    shared_ptr<FaceSetQuality>(new FaceSetQuality(gap, kink, approxtol));
  quality->attach(sfmodel);

  vector<pair<shared_ptr<ftSurface>, shared_ptr<ftSurface> > >  coinc_faces;
  vector<pair<shared_ptr<ftSurface>, shared_ptr<ftSurface> > >  embedded_faces;
  quality->identicalOrEmbeddedFaces(coinc_faces, embedded_faces);

  std::cout << "Number of pairs of identical faces: " << coinc_faces.size() << std::endl;
  std::cout << "Number of pairs of embedded faces: " << embedded_faces.size() << std::endl;

  std::ofstream out_file("identical_sfs.g2");
  size_t ki;
  for (ki=0; ki<coinc_faces.size(); ki++)
  {
      shared_ptr<ParamSurface> surf1 = coinc_faces[ki].first->surface();
      surf1->writeStandardHeader(out_file);
      surf1->write(out_file);
      shared_ptr<ParamSurface> surf2 = coinc_faces[ki].second->surface();
      surf2->writeStandardHeader(out_file);
      surf2->write(out_file);
  }
  std::ofstream out_file2("embedded_sfs.g2");
  for (ki=0; ki<embedded_faces.size(); ki++)
  {
      shared_ptr<ParamSurface> surf1 = embedded_faces[ki].first->surface();
      surf1->writeStandardHeader(out_file2);
      surf1->write(out_file2);
      shared_ptr<ParamSurface> surf2 = embedded_faces[ki].second->surface();
      surf2->writeStandardHeader(out_file2);
      surf2->write(out_file2);
  }
  if (do_repair)
    {
      shared_ptr<FaceSetRepair> repair = 
	shared_ptr<FaceSetRepair>(new FaceSetRepair(quality));
 
      repair->identicalAndEmbeddedFaces();

    }

  vector<pair<shared_ptr<ftSurface>, shared_ptr<ftSurface> > >  coinc_faces2;
  vector<pair<shared_ptr<ftSurface>, shared_ptr<ftSurface> > >  embedded_faces2;
  quality->identicalOrEmbeddedFaces(coinc_faces2, embedded_faces2);


  std::ofstream out_filen("identical_sfs2.g2");
  for (ki=0; ki<coinc_faces2.size(); ki++)
  {
      shared_ptr<ParamSurface> surf1 = coinc_faces2[ki].first->surface();
      surf1->writeStandardHeader(out_filen);
      surf1->write(out_filen);
      shared_ptr<ParamSurface> surf2 = coinc_faces2[ki].second->surface();
      surf2->writeStandardHeader(out_filen);
      surf2->write(out_filen);
  }
  std::ofstream out_filen2("embedded_sfs2.g2");
  for (ki=0; ki<embedded_faces2.size(); ki++)
  {
      shared_ptr<ParamSurface> surf1 = embedded_faces2[ki].first->surface();
      surf1->writeStandardHeader(out_filen2);
      surf1->write(out_filen2);
      shared_ptr<ParamSurface> surf2 = embedded_faces2[ki].second->surface();
      surf2->writeStandardHeader(out_filen2);
      surf2->write(out_filen2);
  }

}
