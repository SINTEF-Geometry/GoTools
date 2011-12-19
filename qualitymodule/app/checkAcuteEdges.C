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

//   // Alternative set of tolerances related to
//   // StepReader/isoplugin. @jbt
//   double gap = 0.0001;
//   double neighbour = 0.001;
//   double kink = 0.01;
//   double approxtol = 0.0001;
  
  CompositeModelFactory factory(approxtol, gap, neighbour, kink, 10.0*kink);

  shared_ptr<CompositeModel> model = shared_ptr<CompositeModel>(factory.createFromG2(file1));
  shared_ptr<SurfaceModel> sfmodel = dynamic_pointer_cast<SurfaceModel,CompositeModel>(model);

  FaceSetQuality quality(gap, kink, approxtol);
  quality.attach(sfmodel);

  vector<pair<ftEdge*, ftEdge*> > acute;
  quality.acuteEdgeAngle(acute);  

  std::cout << "Number of acute edges: " << acute.size() << std::endl;

  std::ofstream out_file("acute_edge.g2");
  size_t ki;
  for (ki=0; ki<acute.size(); ++ki)
  {
      shared_ptr<ParamCurve> crv1 = acute[ki].first->geomCurve();
      crv1->writeStandardHeader(out_file);
      crv1->write(out_file);
      crv1 = acute[ki].second->geomCurve();
      crv1->writeStandardHeader(out_file);
      crv1->write(out_file);
  }

  vector<pair<ftEdge*, ftEdge*> > acute2;
  quality.acuteEdgeAngle(acute2);  

  std::ofstream out_file2("acute_edge2.g2");
  for (ki=0; ki<acute2.size(); ++ki)
  {
      shared_ptr<ParamCurve> crv1 = acute2[ki].first->geomCurve();
      shared_ptr<SplineCurve> tmp = 
	shared_ptr<SplineCurve>(crv1->geometryCurve());
      tmp->writeStandardHeader(out_file2);
      tmp->write(out_file2);
      crv1 = acute2[ki].second->geomCurve();
      tmp = shared_ptr<SplineCurve>(crv1->geometryCurve());
      tmp->writeStandardHeader(out_file2);
      tmp->write(out_file2);
  }


}

