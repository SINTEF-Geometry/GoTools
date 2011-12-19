#ifdef __BORLANDC__
#include <vcl.h>
#endif

#include "GoTools/compositemodel/CompositeCurve.h"
#include "GoTools/compositemodel/CompositeModelFactory.h"
#include "GoTools/geometry/LineCloud.h"
#include <fstream>
#include <stdlib.h> // For atof()

using namespace std;
using namespace Go;

int main( int argc, char* argv[] )
{
  if (argc != 4) {
    std::cout << "Input parameters : Input file, n, density"  << std::endl;
    exit(-1);
  }

  // Read input arguments
  std::ifstream file1(argv[1]);
  ALWAYS_ERROR_IF(file1.bad(), "Input file not found or file corrupt");

  double gap = 0.001;
  double neighbour = 0.01;
  double kink = 0.01;
  // Removed from corresponding use of fairingToolbox
  double approx = 0.001;
  int n = atoi(argv[2]);
  double density = atof(argv[3]);

  CompositeModelFactory factory(approx, gap, neighbour, kink, 10.0*kink);

  CompositeModel *model;
  model = factory.createFromIges(file1);
  
  CompositeCurve *cvmodel = dynamic_cast<CompositeCurve*>(model);


  std::vector<shared_ptr<GeneralMesh> > meshes;
  cvmodel->tesselate(&n, meshes);

  size_t ki;
  std::ofstream out1("lines1.g2");
  for (ki=0; ki<meshes.size(); ++ki)
    {
      double *nodes = meshes[ki]->vertexArray();
      int num_vx = meshes[ki]->numVertices();
      vector<double> seg;
      Point pt1(&nodes[0],&nodes[3]);
      Point pt2;
      for (int kj=1; kj<num_vx; kj++)
	{
	  pt2 = Point(&nodes[3*kj], &nodes[3*(kj+1)]);
	  seg.insert(seg.end(), pt1.begin(), pt1.end());
	  seg.insert(seg.end(), pt2.begin(), pt2.end());
	  pt1 = pt2;
	}
      LineCloud line_seg(&seg[0], (int)seg.size()/6);
      line_seg.writeStandardHeader(out1);
      line_seg.write(out1);
    }
  
  std::vector<shared_ptr<GeneralMesh> > meshes2;
  cvmodel->tesselate(density, meshes2);

  std::ofstream out2("lines2.g2");
  for (ki=0; ki<meshes2.size(); ++ki)
    {
      double *nodes = meshes2[ki]->vertexArray();
      int num_vx = meshes2[ki]->numVertices();
      vector<double> seg;
      Point pt1(&nodes[0],&nodes[3]);
      Point pt2;
      for (int kj=1; kj<num_vx; kj++)
	{
	  pt2 = Point(&nodes[3*kj], &nodes[3*(kj+1)]);
	  seg.insert(seg.end(), pt1.begin(), pt1.end());
	  seg.insert(seg.end(), pt2.begin(), pt2.end());
	  pt1 = pt2;
	}
      LineCloud line_seg(&seg[0], (int)seg.size()/6);
      line_seg.writeStandardHeader(out2);
      line_seg.write(out2);
     }
  
  int break_point;
  break_point = 1;
}

