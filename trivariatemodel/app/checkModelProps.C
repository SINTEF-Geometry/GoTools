#include <fstream>
#include "GoTools/trivariate/SplineVolume.h"
#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/trivariatemodel/ftVolume.h"
#include "GoTools/trivariatemodel/VolumeModel.h"
#include "GoTools/compositemodel/ftSurface.h"
#include "GoTools/compositemodel/EdgeVertex.h"
#include "GoTools/geometry/Utils.h"
#include "GoTools/trivariate/SurfaceOnVolume.h"
#include "GoTools/geometry/CurveOnSurface.h"

using namespace Go;
using std::cout;
using std::endl;
using std::ifstream;
using std::vector;

int main(int argc, char* argv[] )
{
  if (argc != 2 && argc != 3)
      cout << "Usage: " << "infile2, (Insert knots)" << endl;

  ifstream is2(argv[1]);
  ALWAYS_ERROR_IF(is2.bad(), "Bad or no input filename");
  int insert = 0;
  if (argc == 3)
    insert = atoi(argv[2]);

  vector<shared_ptr<ftVolume> > volumes;
  
  //int ki;
  while (!is2.eof())
    {
      // Read volume from file
      ObjectHeader head;
      is2 >> head;
      shared_ptr<SplineVolume> vol2;
      vol2 = shared_ptr<SplineVolume>(new SplineVolume());
      vol2->read(is2);

      shared_ptr<ParamVolume> pvol
          = dynamic_pointer_cast<ParamVolume, SplineVolume>(vol2);
      volumes.push_back(shared_ptr<ftVolume>(new ftVolume(pvol)));
      // volumes.push_back(shared_ptr<ftVolume>(new ftVolume(vol2)));

      eatwhite(is2);
    }

  double gap = 0.0001; //1.0e-5; //0.0001;
  double neighbour = 0.01; //1.0e-3; //0.01;
  double kink = 0.01;
  double bend = 0.1; // 0.05;
  shared_ptr<VolumeModel> model = 
    shared_ptr<VolumeModel>(new VolumeModel(volumes, gap, neighbour, 
					    kink, bend));

  int nmb1 = model->nmbEntities();
  for (int ki=0; ki<nmb1; ++ki)
    {
      shared_ptr<ftVolume> vol = model->getBody(ki);
      shared_ptr<SplineVolume> svol= model->getSplineVolume(ki);

      shared_ptr<SurfaceModel> sfmodel = vol->getShell(0);
      int nmb2 = sfmodel->nmbEntities();
      vector<shared_ptr<ftSurface> > faces = sfmodel->allFaces();
      for (int kj=0; kj<nmb2; ++kj)
	{
	  vector<shared_ptr<ftEdge> > edges = faces[kj]->getAllEdges();
	  shared_ptr<SurfaceOnVolume> surf =
	    dynamic_pointer_cast<SurfaceOnVolume,ParamSurface>(faces[kj]->surface());
	  for (size_t kr=0; kr<edges.size(); ++kr)
	    {
	      shared_ptr<CurveOnSurface> cv =
		dynamic_pointer_cast<CurveOnSurface,ParamCurve>(edges[kr]->geomCurve());
	      int break_debug3;
	      break_debug3 = 0;
	    }
	  int break_debug2;
	  break_debug2 = 0;
	}
      int break_debug;
      break_debug = 0;
    }

  vector<shared_ptr<EdgeVertex> > radedg;
  model->getRadialEdges(radedg);
  std::cout << "Number of radial edges: " << radedg.size() << std::endl;

  if (radedg.size() > 0)
    {
      vector<ftSurface*> adj_face = radedg[0]->getAdjacentFaces();
      vector<Body*> adj_bod = radedg[0]->getAdjacentBodies();
      std::cout << "Nmb faces: " << adj_face.size();
      std::cout <<", nmb bodies: " << adj_bod.size() << std::endl;

      if (adj_bod.size() > 0)
	{
	  shared_ptr<SurfaceModel> sfmodel = adj_bod[0]->getOuterShell();
	  vector<shared_ptr<ftEdge> > bd_edges = sfmodel->getBoundaryEdges();
	  std::cout << "Shell, nmb boundary edges: " << bd_edges.size() << std::endl;
	  vector<shared_ptr<ftEdge> > inner_edges = 
	    sfmodel->getUniqueInnerEdges();
	  std::cout << "Shell, nmb inner edges: " << inner_edges.size() << std::endl;
	  if (inner_edges.size() > 0)
	    {
	      vector<ftSurface*> adj_face = inner_edges[0]->getAdjacentFaces();
	      std::cout << "Nmb faces: " << adj_face.size() << std::endl;
	    }
	}
    }

  vector<shared_ptr<ftEdge> > edg;
  model->uniqueNonRadialEdges(edg);
  std::cout << "Number of non-radial edges: " << edg.size() << std::endl;
  if (edg.size() > 0)
    {
      vector<ftSurface*> adj_face = edg[0]->getAdjacentFaces();
      std::cout << "Nmb faces: " << adj_face.size() << std::endl;
    }

  vector<shared_ptr<ftSurface> > bd_faces = model->getBoundaryFaces();
  std::cout << "Number of boundary faces: " << bd_faces.size() << std::endl;
  if (bd_faces.size() > 0)
    {
      vector<Body*> adj_bod = bd_faces[0]->getAdjacentBodies();
      std::cout <<"Nmb bodies: " << adj_bod.size() << std::endl;
    }
      

  vector<shared_ptr<ftSurface> > inner_faces = model->getUniqueInnerFaces();
  std::cout << "Number of inner faces: " << inner_faces.size() << std::endl;
  if (inner_faces.size() > 0)
    {
      vector<Body*> adj_bod = inner_faces[0]->getAdjacentBodies();
      std::cout <<"Nmb bodies: " << adj_bod.size() << std::endl;
    }

  vector<shared_ptr<Vertex> > vx;
  model->getAllVertices(vx);
  std::cout << "Number of vertices: " << vx.size() << std::endl;

}
