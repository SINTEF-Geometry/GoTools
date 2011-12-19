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
  double kink = 0.01;
  double approxtol = 0.01;

  CompositeModelFactory factory(approxtol, gap, neighbour, kink, 10.0*kink);

  shared_ptr<CompositeModel> model = shared_ptr<CompositeModel>(factory.createFromG2(file1));
  shared_ptr<SurfaceModel> sfmodel = dynamic_pointer_cast<SurfaceModel,CompositeModel>(model);

  FaceSetQuality quality(gap, kink, approxtol);
  quality.attach(sfmodel);

  vector<pair<shared_ptr<ftEdge>, shared_ptr<ftEdge> > >  coinc_edges;
  vector<pair<shared_ptr<ftEdge>, shared_ptr<ftEdge> > >  embedded_edges;
  quality.identicalOrEmbeddedEdges(coinc_edges, embedded_edges);

  std::cout << "Number of pairs of identical edges: " << coinc_edges.size() << std::endl;
  std::cout << "Number of pairs of embedded edges: " << embedded_edges.size() << std::endl;

  std::ofstream out_file("identical_cvs.g2");
  size_t ki, kj;
  vector<shared_ptr<ParamSurface> > edge_id_sfs;
  for (ki=0; ki<coinc_edges.size(); ki++)
  {
      shared_ptr<ParamCurve> crv1 = coinc_edges[ki].first->geomCurve();
      crv1->writeStandardHeader(out_file);
      crv1->write(out_file);
      shared_ptr<ParamCurve> crv2 = coinc_edges[ki].second->geomCurve();
      crv2->writeStandardHeader(out_file);
      crv2->write(out_file);
      shared_ptr<ParamSurface> surf = coinc_edges[ki].first->face()->surface();
      for (kj=0; kj<edge_id_sfs.size(); kj++)
	  if (edge_id_sfs[kj].get() == surf.get())
	      break;
      if (kj == edge_id_sfs.size())
	  edge_id_sfs.push_back(surf);
      surf = coinc_edges[ki].second->face()->surface();
      for (kj=0; kj<edge_id_sfs.size(); kj++)
	  if (edge_id_sfs[kj].get() == surf.get())
	      break;
      if (kj == edge_id_sfs.size())
	  edge_id_sfs.push_back(surf);
   }
  std::ofstream out_file2("embedded_cvs.g2");
  for (ki=0; ki<embedded_edges.size(); ki++)
  {
      shared_ptr<ParamCurve> crv1 = embedded_edges[ki].first->geomCurve();
      crv1->writeStandardHeader(out_file2);
      crv1->write(out_file2);
      shared_ptr<ParamCurve> crv2 = embedded_edges[ki].second->geomCurve();
      crv2->writeStandardHeader(out_file2);
      crv2->write(out_file2);
      shared_ptr<ParamSurface> surf = embedded_edges[ki].first->face()->surface();
      for (kj=0; kj<edge_id_sfs.size(); kj++)
	  if (edge_id_sfs[kj].get() == surf.get())
	      break;
      if (kj == edge_id_sfs.size())
	  edge_id_sfs.push_back(surf);
      surf = embedded_edges[ki].second->face()->surface();
      for (kj=0; kj<edge_id_sfs.size(); kj++)
	  if (edge_id_sfs[kj].get() == surf.get())
	      break;
      if (kj == edge_id_sfs.size())
	  edge_id_sfs.push_back(surf);
  }
  vector<pair<shared_ptr<ftEdge>, shared_ptr<ftEdge> > >  coinc_edges2;
  vector<pair<shared_ptr<ftEdge>, shared_ptr<ftEdge> > >  embedded_edges2;
  quality.identicalOrEmbeddedEdges(coinc_edges2, embedded_edges2);

  std::ofstream out_filen("identical_cvs2.g2");
  for (ki=0; ki<coinc_edges2.size(); ki++)
  {
      shared_ptr<ParamCurve> crv1 = coinc_edges2[ki].first->geomCurve();
      crv1->writeStandardHeader(out_filen);
      crv1->write(out_filen);
      shared_ptr<ParamCurve> crv2 = coinc_edges2[ki].second->geomCurve();
      crv2->writeStandardHeader(out_filen);
      crv2->write(out_filen);
   }
  std::ofstream out_filen2("embedded_cvs2.g2");
  for (ki=0; ki<embedded_edges2.size(); ki++)
  {
      shared_ptr<ParamCurve> crv1 = embedded_edges2[ki].first->geomCurve();
      crv1->writeStandardHeader(out_filen2);
      crv1->write(out_filen2);
      shared_ptr<ParamCurve> crv2 = embedded_edges2[ki].second->geomCurve();
      crv2->writeStandardHeader(out_filen2);
      crv2->write(out_filen2);
  }

  std::ofstream out_file3("id_edge_sfs.g2");
  for (ki=0; ki<edge_id_sfs.size(); ki++)
  {
      edge_id_sfs[ki]->writeStandardHeader(out_file3);
      edge_id_sfs[ki]->write(out_file3);
  }

}
