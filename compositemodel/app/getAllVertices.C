#include "GoTools/compositemodel/SurfaceModel.h"
#include "GoTools/compositemodel/CompositeModelFactory.h"
#include <fstream>

//using namespace std;
using std::shared_ptr;
using std::dynamic_pointer_cast;
using namespace Go;

int main( int argc, char* argv[] )
{
  if (argc != 2) {
    std::cout << "Input parameters : Input file on g2 format" << std::endl;
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

  vector<shared_ptr<Vertex> > vertices;
  sfmodel->getAllVertices(vertices);

  std::ofstream out_file("vertices.g2");
  for (size_t ki=0; ki<vertices.size(); ki++)
  {
      Point pnt1 = vertices[ki]->getVertexPoint();

      out_file << "400 1 0 4 255 0 0 255 \n";
      out_file << "1 \n";
      out_file << pnt1[0]  << " " << pnt1[1] << " " << pnt1[2] << "\n";
  }
}

