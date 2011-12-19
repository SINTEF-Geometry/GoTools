#include "GoTools/compositemodel/SurfaceModel.h"
#include "GoTools/compositemodel/CompositeModelFactory.h"
#include "GoTools/qualitymodule/FaceSetQuality.h"
#include <fstream>

using std::vector;
using namespace Go;

int main( int argc, char* argv[] )
{
  if (argc != 2) {
    std::cout << "Input parameters : Input file on g2 format, " << std::endl;
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

  vector<shared_ptr<ParamCurve> >  g1_crvs;
  vector<shared_ptr<ParamCurve> >  c1_crvs;
  quality.cvC1G1Discontinuity(c1_crvs, g1_crvs);

  std::cout << "Number c1 discontinuity crvs: " << c1_crvs.size() << std::endl;
  std::cout << "Number g1 discontinuity crvs: " << g1_crvs.size() << std::endl;

  std::ofstream out_file("c1dis_cvs.g2");
  size_t ki;
  for (ki=0; ki<c1_crvs.size(); ki++)
  {
      c1_crvs[ki]->writeStandardHeader(out_file);
      c1_crvs[ki]->write(out_file);
   }

  std::ofstream out_file2("g1dis_cvs.g2");
  for (ki=0; ki<g1_crvs.size(); ki++)
  {
      g1_crvs[ki]->writeStandardHeader(out_file2);
      g1_crvs[ki]->write(out_file2);
   }

  vector<shared_ptr<ParamCurve> >  g1_crvs2;
  vector<shared_ptr<ParamCurve> >  c1_crvs2;
  quality.cvC1G1Discontinuity(c1_crvs2, g1_crvs2);


  std::ofstream out_filen("c1dis_cvs2.g2");
  for (ki=0; ki<c1_crvs2.size(); ki++)
  {
      c1_crvs2[ki]->writeStandardHeader(out_filen);
      c1_crvs2[ki]->write(out_filen);
   }

  std::ofstream out_filen2("g1dis_cvs2.g2");
  for (ki=0; ki<g1_crvs2.size(); ki++)
  {
      g1_crvs2[ki]->writeStandardHeader(out_filen2);
      g1_crvs2[ki]->write(out_filen2);
   }

}
