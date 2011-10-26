#include "GoTools/compositemodel/SurfaceModel.h"
#include "GoTools/compositemodel/CompositeModelFactory.h"
#include "GoTools/tesselator/LineStrip.h"
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

  if (sfmodel.get())
    {
      int nmb_bd = sfmodel->nmbBoundaries();
      int ki;
      std::ofstream out1("boundaries.g2");
      for (ki=0; ki<nmb_bd; ++ki)
	{
	  ftCurve bd = sfmodel->getBoundary(ki);
	  std::vector<shared_ptr<LineStrip> > bd_mesh;
	  bd.tesselate(20, bd_mesh);
	  for (size_t kj=0; kj<bd_mesh.size(); ++kj)
	    {
	      LineCloud lines;
	      std::vector<double> tmp_lines;
	      int nmb_vert = bd_mesh[kj]->numVertices();
	      double *vertices = bd_mesh[kj]->vertexArray();
	      for (int kr=0; kr<nmb_vert-1; ++kr)
		tmp_lines.insert(tmp_lines.end(), vertices+kr*3, 
				 vertices+(kr+2)*3);

	      lines.setCloud(&tmp_lines[0], nmb_vert-1);
	      lines.writeStandardHeader(out1);
	      lines.write(out1);
	    }
	}

      std::ofstream out2("g1Disc.g2");
      ftCurve g1disc = sfmodel->getG1Disconts();
      std::vector<shared_ptr<LineStrip> > g1d_mesh;
      g1disc.tesselate(20, g1d_mesh);
      for (size_t kj=0; kj<g1d_mesh.size(); ++kj)
	{
	  LineCloud lines;
	  vector<double> tmp_lines;
	  int nmb_vert = g1d_mesh[kj]->numVertices();
	  double *vertices = g1d_mesh[kj]->vertexArray();
	  for (int kr=0; kr<nmb_vert-1; ++kr)
	    tmp_lines.insert(tmp_lines.end(), vertices+kr*3, 
			     vertices+(kr+2)*3);

	  lines.setCloud(&tmp_lines[0], nmb_vert-1);
	  lines.writeStandardHeader(out2);
	  lines.write(out2);
	}
    }
}

