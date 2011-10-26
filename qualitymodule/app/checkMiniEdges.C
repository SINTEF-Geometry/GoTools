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
  if (argc != 3) {
    std::cout << "Input parameters : Input file on g2 format, mini edge size" << std::endl;
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

  double mini = atof(argv[2]);
  quality.setMiniElementSize(mini);

  vector<shared_ptr<ftEdge> >  mini_edges;
  quality.miniEdges(mini_edges);

  std::cout << "Number mini edges: " << mini_edges.size() << std::endl;

  std::ofstream out_file("mini_cvs.g2");
  size_t ki;
  for (ki=0; ki<mini_edges.size(); ki++)
  {
      shared_ptr<ParamCurve> crv1 = mini_edges[ki]->geomCurve();
      crv1->writeStandardHeader(out_file);
      crv1->write(out_file);
   }
  vector<shared_ptr<ftEdge> >  mini_edges2;
  quality.miniEdges(mini_edges2);


  std::ofstream out_file2("mini_cvs2.g2");
  for (ki=0; ki<mini_edges2.size(); ki++)
  {
      shared_ptr<ParamCurve> crv1 = mini_edges2[ki]->geomCurve();
      crv1->writeStandardHeader(out_file2);
      crv1->write(out_file2);
   }

}
