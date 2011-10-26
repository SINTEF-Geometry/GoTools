#include "GoTools/compositemodel/SurfaceModel.h"
#include "GoTools/compositemodel/CompositeModelFactory.h"
#include "GoTools/qualitymodule/FaceSetQuality.h"
#include <fstream>

using std::vector;
using std::endl;
using std::shared_ptr;
using namespace Go;

int main( int argc, char* argv[] )
{
  if (argc != 5 && argc != 6) {
    std::cout << "Input parameters : Input file on g2 format, get edge-vertex, get face-vertex, get face-edge [refinement tolerance factor]" << std::endl;
    exit(-1);
  }

  // Read input arguments
  std::ifstream file1(argv[1]);
  ALWAYS_ERROR_IF(file1.bad(), "Input file not found or file corrupt");

  bool get_EV = atoi(argv[2]) != 0;
  bool get_FV = atoi(argv[3]) != 0;
  bool get_FE = atoi(argv[4]) != 0;

  double gap = 0.001;
  double neighbour = 0.01;
  double kink = 0.01;
  double approxtol = 0.01;

  CompositeModelFactory factory(approxtol, gap, neighbour, kink, 10.0*kink);

  //shared_ptr<CompositeModel> model = shared_ptr<CompositeModel>(factory.createFromG2(file1));
  //shared_ptr<SurfaceModel> sfmodel = dynamic_pointer_cast<SurfaceModel,CompositeModel>(model);


  shared_ptr<SurfaceModel> sfmodel(dynamic_cast<SurfaceModel*>(factory.createFromG2(file1)));

  FaceSetQuality quality(gap, kink, approxtol);
  quality.attach(sfmodel);

  //  if (argc == 6)
  //  quality.manipulategap(atof(argv[5]));

  if (get_EV)
    {
      vector<pair<ftEdge*, shared_ptr<Vertex> > > edge_vertices;
      quality.edgeVertexDistance(edge_vertices);
      std::ofstream out_file("distance_EV.g2");
      for (size_t i = 0; i < edge_vertices.size(); i++)
	{
	  SplineCurve* sc = edge_vertices[i].first->geomCurve()->geometryCurve();
	  sc->writeStandardHeader(out_file);
	  sc->write(out_file);
	  Point pt = edge_vertices[i].second->getVertexPoint();
	  out_file << "400 1 0 4 255 0 0 255" << endl;
	  out_file << pt << endl;
	}
      std::cout << "Found " << edge_vertices.size() << " fair distance edge/vertex pairs" << endl;

      vector<pair<ftEdge*, shared_ptr<Vertex> > > edge_vertices2;
      quality.edgeVertexDistance(edge_vertices2);
      std::ofstream out_file2("distance_EV2.g2");
      for (size_t i = 0; i < edge_vertices2.size(); i++)
	{
	  SplineCurve* sc = edge_vertices2[i].first->geomCurve()->geometryCurve();
	  sc->writeStandardHeader(out_file2);
	  sc->write(out_file2);
	  Point pt = edge_vertices2[i].second->getVertexPoint();
	  out_file2 << "400 1 0 4 255 0 0 255" << endl;
	  out_file2 << pt << endl;
	}
    }

  if (get_FV)
    { 
      vector<pair<ftSurface*, shared_ptr<Vertex> > > face_vertices;
      quality.faceVertexDistance(face_vertices);
      std::ofstream out_file("distance_FV.g2");
      for (size_t i = 0; i < face_vertices.size(); i++)
	{
	  shared_ptr<ParamSurface> ps = face_vertices[i].first->asFtSurface()->surface();
	  ps->writeStandardHeader(out_file);
	  ps->write(out_file);
	  Point pt = face_vertices[i].second->getVertexPoint();
	  out_file << "400 1 0 4 255 0 0 255" << endl;
	  out_file << pt << endl;
	}
      std::cout << "Found " << face_vertices.size() << " fair distance face/vertex pairs" << endl;

      vector<pair<ftSurface*, shared_ptr<Vertex> > > face_vertices2;
      quality.faceVertexDistance(face_vertices2);
      std::ofstream out_file2("distance_FV2.g2");
      for (size_t i = 0; i < face_vertices2.size(); i++)
	{
	  shared_ptr<ParamSurface> ps = face_vertices2[i].first->asFtSurface()->surface();
	  ps->writeStandardHeader(out_file2);
	  ps->write(out_file2);
	  Point pt = face_vertices2[i].second->getVertexPoint();
	  out_file2 << "400 1 0 4 255 0 0 255" << endl;
	  out_file2 << pt << endl;
	}
    }

  if (get_FE)
    {
      vector<pair<ftSurface*, ftEdge*> > face_edges;
      quality.faceEdgeDistance(face_edges);
      std::ofstream out_file("distance_FE.g2");
      for (size_t i = 0; i < face_edges.size(); i++)
	{
	  shared_ptr<ParamSurface> ps = face_edges[i].first->asFtSurface()->surface();
	  ps->writeStandardHeader(out_file);
	  ps->write(out_file);
	  SplineCurve* sc = face_edges[i].second->geomCurve()->geometryCurve();
	  sc->writeStandardHeader(out_file);
	  sc->write(out_file);
	}
      std::cout << "Found " << face_edges.size() << " fair distance face/edge pairs" << endl;

      vector<pair<ftSurface*, ftEdge*> > face_edges2;
      quality.faceEdgeDistance(face_edges2);
      std::ofstream out_file2("distance_FE2.g2");
      for (size_t i = 0; i < face_edges2.size(); i++)
	{
	  shared_ptr<ParamSurface> ps = face_edges2[i].first->asFtSurface()->surface();
	  ps->writeStandardHeader(out_file2);
	  ps->write(out_file2);
	  SplineCurve* sc = face_edges2[i].second->geomCurve()->geometryCurve();
	  sc->writeStandardHeader(out_file2);
	  sc->write(out_file2);
	}
    }

}
