#include "GoTools/compositemodel/SurfaceModel.h"
#include "GoTools/compositemodel/CompositeModelFactory.h"
#include <fstream>

//using namespace std;
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

  std::ofstream out_file("bd_cvs.g2");
  vector<shared_ptr<ParamCurve> > bd_cvs;
  int nmb_sfs = sfmodel->nmbEntities();
  for (int ki=0; ki<nmb_sfs; ++ki)
    {
      shared_ptr<ftSurface> face = sfmodel->getFace(ki);
      vector<shared_ptr<ftEdgeBase> > edges = face->createInitialEdges();
      for (size_t kj=0; kj<edges.size(); ++kj)
	{
	  ftEdge *curr_edge = edges[kj]->geomEdge();
	  shared_ptr<ParamCurve> curr_crv = curr_edge->geomCurve();
	  shared_ptr<ParamCurve> curr_crv2 = 
	    shared_ptr<ParamCurve>(curr_crv->subCurve(curr_edge->tMin(),
						      curr_edge->tMax()));
	  shared_ptr<CurveOnSurface> sf_cv = 
	    dynamic_pointer_cast<CurveOnSurface, ParamCurve>(curr_crv2);
	  if (sf_cv.get())
	    {
	      sf_cv->spaceCurve()->writeStandardHeader(out_file);
	      sf_cv->spaceCurve()->write(out_file);
	    }
	  else
	    {
	      curr_crv2->writeStandardHeader(out_file);
	      curr_crv2->write(out_file);
	    }
	}
    }
}
