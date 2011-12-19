//===========================================================================
//                                                                           
// File: CurveModel.C                                                    
//                                                                           
// Created: June 2009
//                                                                           
// Author: Vibeke Skytt
//                                                                           
// Revision: 
//                                                                           
// Description:
//                                                                           
//===========================================================================

#include "GoTools/compositemodel/CurveModel.h"
#include "GoTools/compositemodel/ftEdge.h"
#include "GoTools/compositemodel/Vertex.h"
#include "GoTools/compositemodel/CompositeCurve.h"
#include "GoTools/utils/CurvatureUtils.h"
#include "GoTools/tesselator/CurveTesselator.h"
#include "GoTools/tesselator/TesselatorUtils.h"

using std::vector;

namespace Go
{
  //===========================================================================
  CurveModel::CurveModel(double gap,   // Gap between adjacent curves
			 double neighbour,  // Threshold for whether curves are adjacent
			 double kink,  // Kink between adjacent curves
			 double bend, // Intended G1 discontinuity between adjacent curves
			 vector<shared_ptr<ParamCurve> >& curves)
  //===========================================================================
    : CompositeModel(gap, neighbour, kink, bend)
  {
    if (curves.empty())
      return;

    // Make edges based on the input curves
    edges_.reserve(curves.size());
    for (size_t ki=0; ki<curves.size(); ++ki)
      edges_.push_back(shared_ptr<ftEdge>(new ftEdge(curves[ki],
						     curves[ki]->startparam(),
						     curves[ki]->endparam())));
    // Make topology
    buildTopology();

  }

  //===========================================================================
  CurveModel::~CurveModel()
  //===========================================================================
  {
    // Empty destructor
  }

  //===========================================================================
  CurveModel* CurveModel::clone() const
  //===========================================================================
  {
    vector<shared_ptr<ParamCurve> > curves(edges_.size());
    for (size_t ki=0; ki<edges_.size(); ++ki)
      curves[ki] = shared_ptr<ParamCurve>(edges_[ki]->geomCurve()->clone());

    return new CurveModel(toptol_.gap, toptol_.neighbour,
			  toptol_.kink, toptol_.bend, curves);
  }

  //===========================================================================
  int CurveModel::nmbEntities() const
  //===========================================================================
  {
      return (int)edges_.size();
  }

  //===========================================================================
  shared_ptr<ParamCurve> CurveModel::getCurve(int idx) const
  //===========================================================================
  {
    return edges_[idx]->geomCurve();
  }

  //===========================================================================
  int CurveModel::getIndex(ParamCurve* curve) const
  //===========================================================================
  {
    for (size_t ki=0; ki<edges_.size(); ++ki)
      if (edges_[ki]->geomCurve().get() == curve)
	return (int)ki;

    return -1;
  }

  //===========================================================================
  void 
  CurveModel::evaluate(int idx,      // Index
		       double par[], // Parameter value
		       Point& pnt) const
  //===========================================================================
  {
    pnt = edges_[idx]->point(par[0]);
  }

  //===========================================================================
  void 
  CurveModel::evaluate(int idx,      // Index
		       double par[], // Parameter value
		       int nder,     // Number of derivatives to compute, 0=only position
		       std::vector<Point>& der) const
  //===========================================================================
  {
    edges_[idx]->geomCurve()->point(der, par[0], nder);
  }

//===========================================================================
void 
CurveModel::closestPoint(Point& pnt,     // Input point
			 Point& clo_pnt, // Found closest point
			 int& idx,           // Index of curve where the closest point is found
			 double clo_par[],   // Parameter value corresponding to the closest point
			 double& dist)       // Distance between input point and found closest point
//===========================================================================
  {
  }

//===========================================================================
shared_ptr<IntResultsModel> CurveModel::intersect(const ftLine& line)
//===========================================================================
{
	return shared_ptr<IntResultsModel>();
}

//===========================================================================
shared_ptr<IntResultsModel> CurveModel::intersect_plane(const ftPlane& plane)
//===========================================================================
{
	return shared_ptr<IntResultsModel>();
}
//===========================================================================
void 
CurveModel::extremalPoint(Point& dir,     // Direction
			  Point& ext_pnt, // Found extremal point
			  int& idx,       // Index of curve where the extremal point is found
			  double ext_par[]) 
//===========================================================================
{
	return;
}

//===========================================================================
BoundingBox CurveModel::boundingBox() 
//===========================================================================
{
  BoundingBox box;
  if (edges_.size() == 0)
    return box;

  box = edges_[0]->geomCurve()->boundingBox();
  for (size_t ki=0; ki<edges_.size(); ++ki)
    {
      BoundingBox bb = edges_[ki]->geomCurve()->boundingBox();
      box.addUnionWith(bb);
    }

  return box;
}

//===========================================================================
BoundingBox CurveModel::boundingBox(int idx) const
//===========================================================================
{
  return edges_[idx]->geomCurve()->boundingBox();
}

//===========================================================================
bool CurveModel::isDegenerate(int idx) const
//===========================================================================
{
  if (idx < 0 || idx >= (int)edges_.size())
    return true;  // Illegal index
  return edges_[idx]->geomCurve()->isDegenerate(toptol_.gap);
}

//===========================================================================
double CurveModel::curvature(int idx, // Index of curve
			     double *par) const
//===========================================================================
{
  vector<Point> der(3), unitder(3);
  evaluate(idx, par, 2, der);
  (void)curvatureRadius(der, unitder);
  double curvature = unitder[2].length();
  return curvature;
}

//===========================================================================
void CurveModel::turn(int idx)
//===========================================================================
{
}

//===========================================================================
void CurveModel::turn()
//===========================================================================
{
}

//===========================================================================
void CurveModel::tesselate(vector<shared_ptr<GeneralMesh> >& meshes) const
//===========================================================================
{
  int res = 100;
  tesselate(&res, meshes);
}

//===========================================================================
  void CurveModel::tesselate(int resolution[],
			     vector<shared_ptr<GeneralMesh> >& meshes) const
//===========================================================================
{
  meshes.clear();
  for (size_t ki=0; ki<edges_.size(); ++ki)
    {
      CurveTesselator tesselator(*edges_[ki]->geomCurve().get());
      tesselator.changeRes(resolution[0]);
      shared_ptr<GeneralMesh> mesh = tesselator.getMesh();
      meshes.push_back(mesh);
    }
}

//===========================================================================
  void CurveModel::tesselate(double density,
			     vector<shared_ptr<GeneralMesh> >& meshes) const
//===========================================================================
{
  int min_nmb = 5;
  int max_nmb = (int)(1000000.0/(int)edges_.size());

  for (size_t ki=0; ki<edges_.size(); ++ki)
    {
      shared_ptr<ParamCurve> curve = edges_[ki]->geomCurve();
      double len = curve->estimatedCurveLength();
      int res = (int)(len/density);
      res = std::max(min_nmb, std::min(res, max_nmb));

      CurveTesselator tesselator(*curve.get());
      tesselator.changeRes(res);
      shared_ptr<GeneralMesh> mesh = tesselator.getMesh();
      meshes.push_back(mesh);
    }

}

  //===========================================================================
  void CurveModel::tesselatedCtrPolygon(vector<shared_ptr<LineCloud> >& ctr_pol) const
  //===========================================================================
  {
    for (size_t ki=0; ki<edges_.size(); ++ki)
      {
	shared_ptr<LineCloud> curr_pol = TesselatorUtils::getCtrPol(edges_[ki]->geomCurve().get());
	ctr_pol.push_back(curr_pol);
      }
  }

  bool nmb_vertex(Vertex* v1, Vertex* v2)
  {
    if (v1->nmbUniqueEdges() == 2)
      return false;
    else if (v2->nmbUniqueEdges() == 2)
      return true;
    else if (v1->nmbUniqueEdges() < v2->nmbUniqueEdges())
      return true;
    else
      return false;
  }

 //===========================================================================
  vector<shared_ptr<CompositeCurve> > CurveModel::fetchCompositeCurves() const
  //===========================================================================
  {
    // Fetch all vertices bounding the edges
    std::set<Vertex*> all_vertices;  // All vertices in the model 
                                                 // represented once
    size_t ki;
    for (ki=0; ki<edges_.size(); ++ki)
      {
	vector<Vertex*>  curr_vertices(2);
	curr_vertices[0] = edges_[ki]->getVertex(true).get();
	curr_vertices[1] = edges_[ki]->getVertex(false).get();
	all_vertices.insert(curr_vertices.begin(), curr_vertices.end());
      }
    
    // Remove all vertices, not being an endpoint of a composite curve
    vector<Vertex*> vertices(all_vertices.begin(), all_vertices.end());
    std::sort(vertices.begin(), vertices.end(), nmb_vertex);

    vector<shared_ptr<CompositeCurve> > composite_curves;
    vector<ftEdge*> collected_edges;
    Vertex *vx;
    int kj;
    for (ki=0; ki<vertices.size(); ++ki)
      {
	int nmb_crvs = vertices[ki]->nmbUniqueEdges();
	for (kj=0; kj<nmb_crvs; ++kj)
	  {
	    // Check if the edge is used already
	    // Remember this is a curve model. Only one half edge exists
	    // for each unique edge
	    vx = vertices[ki];
	    ftEdge *first = vx->getEdge(kj);
	    ftEdge *curr = first;
	    vector<ftEdge*>::iterator res = 
	      std::find(collected_edges.begin(), collected_edges.end(), curr);
	    if (res != collected_edges.end())
	      continue;   // Path already collected

	    collected_edges.push_back(curr); // Save start edge of path

	    vector<shared_ptr<ParamCurve> > path;
	    path.push_back(curr->geomCurve());
	    while (true)
	      {
		shared_ptr<Vertex> v1, v2;
		curr->getVertices(v1, v2);
		Vertex *next = (v1.get() == vx) ? v2.get() : v1.get();
		int nmb = next->nmbUniqueEdges();
		if (nmb != 2)
		  break;  // We have reached the end of the path

		vector<ftEdge*> next_edges = next->uniqueEdges();
		curr = (next_edges[0] == curr) ? next_edges[1] : next_edges[0];
		if  (curr == first)
		  break;  // We have reached the start edge of a loop

		collected_edges.push_back(curr);
		path.push_back(curr->geomCurve());
		vx = next;
	      }

	    CompositeCurve *comp_crv = 
	      new CompositeCurve(toptol_.gap, toptol_.neighbour, 
				 toptol_.kink, toptol_.bend, path);
	    composite_curves.push_back(shared_ptr<CompositeCurve>(comp_crv));
	  }
      }
    return composite_curves;
  }

  //===========================================================================
  void CurveModel::buildTopology()
  //===========================================================================
  {
    // Fetch all vertices bounding the edges
    std::set<shared_ptr<Vertex> > all_vertices;  // All vertices in the model 
                                       // represented once
    size_t ki, kj;
    for (ki=0; ki<edges_.size(); ++ki)
      {
	vector<shared_ptr<Vertex> > curr_vertices(2);
	curr_vertices[0] = edges_[ki]->getVertex(true);
	curr_vertices[1] = edges_[ki]->getVertex(false);
	all_vertices.insert(curr_vertices.begin(), curr_vertices.end());
      }
    
    // Check distance bewteen all pairs of vertices
    vector<shared_ptr<Vertex> > vertices(all_vertices.begin(), 
					 all_vertices.end());
    for (ki=0; ki<vertices.size(); ++ki)
      for (kj=ki+1; kj<vertices.size(); ++kj)
	{
	  double dist = 
	    vertices[ki]->getVertexPoint().dist(vertices[kj]->getVertexPoint());
	  if (dist < toptol_.neighbour)
	    {
	      vector<ftEdge*> connected_edges = vertices[kj]->allEdges();
	      vertices[ki]->joinVertex(vertices[kj]);
	      for (size_t kr=0; kr<connected_edges.size(); ++kr)
		connected_edges[kr]->replaceVertex(vertices[kj], vertices[ki]);
	      vertices.erase(vertices.begin() + kj);
	      break;
	    }
	}
    
  }
} // namespace Go
