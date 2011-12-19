#include "GoTools/compositemodel/SurfaceModel.h"
#include "GoTools/compositemodel/CompositeModelFactory.h"
#include "GoTools/qualitymodule/FaceSetQuality.h"
#include "GoTools/qualitymodule/FaceSetRepair.h"
#include <fstream>

using std::vector;
using namespace Go;

int main( int argc, char* argv[] )
{
  if (argc != 3) {
    std::cout << "Input parameters : Input file on g2 format, repair" << std::endl;
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


  vector<pair<shared_ptr<Vertex>, shared_ptr<Vertex> > > coinc_vertex;
  quality->identicalVertices(coinc_vertex);

  std::cout << "Number of pairs of identical vertices: " << coinc_vertex.size() << std::endl;

  std::ofstream out_file("identical_vertices.g2");
  for (size_t ki=0; ki<coinc_vertex.size(); ki++)
  {
      Point pnt1 = coinc_vertex[ki].first->getVertexPoint();
      Point pnt2 = coinc_vertex[ki].second->getVertexPoint();

      out_file << "400 1 0 4 0 255 0 255 \n";
      out_file << "1 \n";
      out_file << pnt1[0]  << " " << pnt1[1] << " " << pnt1[2] << "\n";
      out_file << "400 1 0 4 255 0 0 255 \n";
      out_file << "1 \n";
      out_file << pnt2[0]  << " " << pnt2[1] << " " << pnt2[2] << "\n";
  }

  if (do_repair)
    {
      shared_ptr<FaceSetRepair> repair = 
	shared_ptr<FaceSetRepair>(new FaceSetRepair(quality));
 
      repair->identicalVertices();

    }

  vector<pair<shared_ptr<Vertex>, shared_ptr<Vertex> > > coinc_vertex2;
  quality->identicalVertices(coinc_vertex2);

  std::ofstream out_file2("identical_vertices2.g2");
  for (size_t ki=0; ki<coinc_vertex2.size(); ki++)
  {
      Point pnt1 = coinc_vertex2[ki].first->getVertexPoint();
      Point pnt2 = coinc_vertex2[ki].second->getVertexPoint();

      out_file2 << "400 1 0 4 0 255 0 255 \n";
      out_file2 << "1 \n";
      out_file2 << pnt1[0]  << " " << pnt1[1] << " " << pnt1[2] << "\n";
      out_file2 << "400 1 0 4 255 0 0 255 \n";
      out_file2 << "1 \n";
      out_file2 << pnt2[0]  << " " << pnt2[1] << " " << pnt2[2] << "\n";
  }
}

