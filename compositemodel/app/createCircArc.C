#ifdef __BORLANDC__
#include <vcl.h>
#endif

#include "GoTools/compositemodel/CompositeCurve.h"
#include "GoTools/compositemodel/CompositeModelFactory.h"
#include <fstream>

//using namespace std;
using std::shared_ptr;
using namespace Go;

int main( int argc, char* argv[] )
{
  if (argc != 11) {
      std::cout << "Input parameters : start point (x,y,z), mid point (x,y,z), angle (in degrees), axis (x,y,z)" << std::endl;
    exit(-1);
  }

#ifdef __BORLANDC__
  using Go::Point;
#endif

  Point startpt(atof(argv[1]), atof(argv[2]), atof(argv[3]));
  Point midpt(atof(argv[4]), atof(argv[5]), atof(argv[6]));
  double angle = atof(argv[7]);
  Point axis(atof(argv[8]), atof(argv[9]), atof(argv[10]));
  double gap = 0.001;
  double neighbour = 0.01;
  double kink = 0.01;
  double approxtol = 0.01;

  angle = (angle*M_PI)/180.0;

  CompositeModelFactory factory(approxtol, gap, neighbour, kink, 10.0*kink);

  CompositeModel *model = factory.createCircularArc(midpt, startpt, angle, axis);

  CompositeCurve *cvmodel = dynamic_cast<CompositeCurve*>(model);

  if (cvmodel)
  {
      std::ofstream out_file("circ.g2");
      int nmb = cvmodel->nmbEntities();
      for (int ki=0; ki<nmb; ki++)
      {
	  shared_ptr<ParamCurve> crv = cvmodel->getCurve(ki);

	  crv->writeStandardHeader(out_file);
	  crv->write(out_file);
      }

      out_file << "400 1 0 4 255 0 0 255" << std::endl;
      out_file << "2" << std::endl;
      out_file << midpt << std::endl;
      out_file << startpt << std::endl;


      std::ofstream out1("circ_mesh.g2");
      std::vector<std::shared_ptr<GeneralMesh> > meshes;
      int res[2];
      res[0] = res[1] = 100;
      cvmodel->tesselate(res, meshes);
      for (size_t kj=0; kj<meshes.size(); ++kj)
	{
	  LineCloud lines;
	  vector<double> tmp_lines;
	  int nmb_vert = meshes[kj]->numVertices();
	  double *vertices = meshes[kj]->vertexArray();
	  for (int kj=0; kj<nmb_vert-1; ++kj)
	    tmp_lines.insert(tmp_lines.end(), vertices+kj*3, vertices+(kj+2)*3);

	  lines.setCloud(&tmp_lines[0], nmb_vert-1);
	  lines.writeStandardHeader(out1);
	  lines.write(out1);
	}
  }

  delete model;
}
