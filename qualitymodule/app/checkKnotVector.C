#include "GoTools/compositemodel/SurfaceModel.h"
#include "GoTools/compositemodel/CompositeModelFactory.h"
#include "GoTools/qualitymodule/FaceSetQuality.h"
#include <fstream>

using std::vector;
using namespace Go;

int main( int argc, char* argv[] )
{
  if (argc != 3) {
    std::cout << "Input parameters : Input file on g2 format, parameter tolerance" << std::endl;
    exit(-1);
  }

  // Read input arguments
  std::ifstream file1(argv[1]);
  ALWAYS_ERROR_IF(file1.bad(), "Input file not found or file corrupt");

    double gap = 0.001;
  double neighbour = 0.01;
  double kink = 0.01;
  double approxtol = 0.01;

  double tol = atof(argv[2]);

  CompositeModelFactory factory(approxtol, gap, neighbour, kink, 10.0*kink);

  shared_ptr<CompositeModel> model = shared_ptr<CompositeModel>(factory.createFromG2(file1));
  shared_ptr<SurfaceModel> sfmodel = dynamic_pointer_cast<SurfaceModel,CompositeModel>(model);

  FaceSetQuality quality(gap, kink, approxtol);
  quality.attach(sfmodel);

  vector<shared_ptr<ParamCurve> > cv_knots;
  vector<shared_ptr<ParamSurface> > sf_knots;
  quality.indistinctKnots(cv_knots, sf_knots, tol);

  std::cout << "Trimming curves with indistinct knots: " << cv_knots.size() << std::endl;
  std::cout << "Surfaces with indistinct knots: " << sf_knots.size() << std::endl;

  std::ofstream out_file("indistinct_knots.g2");
  size_t ki;
  for (ki=0; ki<cv_knots.size(); ki++)
  {
      cv_knots[ki]->writeStandardHeader(out_file);
      cv_knots[ki]->write(out_file);
  }

  for (ki=0; ki<sf_knots.size(); ki++)
  {
      sf_knots[ki]->writeStandardHeader(out_file);
      sf_knots[ki]->write(out_file);
  }

  vector<shared_ptr<ParamCurve> > cv_knots2;
  vector<shared_ptr<ParamSurface> > sf_knots2;
  quality.indistinctKnots(cv_knots2, sf_knots2, tol);

  std::ofstream out_file2("indistinct_knots2.g2");
  for (ki=0; ki<cv_knots2.size(); ki++)
  {
      cv_knots2[ki]->writeStandardHeader(out_file2);
      cv_knots2[ki]->write(out_file2);
  }

  for (ki=0; ki<sf_knots2.size(); ki++)
  {
      sf_knots2[ki]->writeStandardHeader(out_file2);
      sf_knots2[ki]->write(out_file2);
  }

}
