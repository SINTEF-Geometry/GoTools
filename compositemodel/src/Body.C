//===========================================================================
//                                                                           
// File: Body.C                                                       
//                                                                           
// Created: September 2008
//                                                                           
// Author: Vibeke Skytt
//                                                                           
// Revision: 
//                                                                           
// Description:
//                                                                           
//===========================================================================

// VSK. September 2008
// Note that body requires a set of test functionality
// Check that each shell is closed
// Check that the shells don't intersect
// Check that the inner void shells is inside the outer shell
// These tests are not implemented yet
// If the contained shells have different topology tolerances, the body gets the
// largest tolerances. Should topological results as neighbourhood, gap, kink etc. 
// be updated with respect to the largest tolerances?
// The functionality is currently very lean, and will be extended when needed

#include "GoTools/compositemodel/Body.h"
#include "GoTools/geometry/SplineCurve.h"

using std::vector;
using std::shared_ptr;

namespace Go
{

//---------------------------------------------------------------------------
Body::Body()
    : toptol_(0.0, 0.0, 0.0, 0.0)
//---------------------------------------------------------------------------
{
}

//---------------------------------------------------------------------------
Body::Body(const CoordinateSystem<3> xyz)
    : coordinate_(xyz), toptol_(0.0, 0.0, 0.0, 0.0)
//---------------------------------------------------------------------------
{
}

//---------------------------------------------------------------------------
Body::Body(const CoordinateSystem<3> xyz, vector<shared_ptr<SurfaceModel> >& shells)
    : coordinate_(xyz), shells_(shells), toptol_(0.0, 0.0, 0.0, 0.0)
//---------------------------------------------------------------------------
{
    double gap=0.0, neighbour=0.0, kink=0.0, bend=0.0;
    for (size_t ki=0; ki<shells.size(); ki++)
    {
	tpTolerances toptol = shells[ki]->getTolerances();
	gap = std::max(toptol.gap, gap);
	neighbour = std::max(toptol.neighbour, neighbour);
	kink = std::max(toptol.kink, kink);
	bend = std::max(toptol.bend, bend);
    }
    toptol_ = tpTolerances(gap, neighbour, kink, bend);

    addBodyPointers();
}

//---------------------------------------------------------------------------
Body::Body(const CoordinateSystem<3> xyz, shared_ptr<SurfaceModel>  shell)
    : coordinate_(xyz), toptol_(shell->getTolerances())
//---------------------------------------------------------------------------
{
    shells_.push_back(shell);
    addBodyPointers();
}

Body::Body(vector<shared_ptr<SurfaceModel> >& shells)
    : shells_(shells), toptol_(0.0, 0.0, 0.0, 0.0)
//---------------------------------------------------------------------------
{
    double gap=0.0, neighbour=0.0, kink=0.0, bend=0.0;
    for (size_t ki=0; ki<shells.size(); ki++)
    {
	tpTolerances toptol = shells[ki]->getTolerances();
	gap = std::max(toptol.gap, gap);
	neighbour = std::max(toptol.neighbour, neighbour);
	kink = std::max(toptol.kink, kink);
	bend = std::max(toptol.bend, bend);
    }
    toptol_ = tpTolerances(gap, neighbour, kink, bend);

    addBodyPointers();
}

//---------------------------------------------------------------------------
Body::Body(shared_ptr<SurfaceModel>  shell)
    : toptol_(shell->getTolerances())
//---------------------------------------------------------------------------
{
    shells_.push_back(shell);
    addBodyPointers();
}

//---------------------------------------------------------------------------
Body::~Body()
//---------------------------------------------------------------------------
{
}

//---------------------------------------------------------------------------
void Body::addshell(shared_ptr<SurfaceModel> shell)
//---------------------------------------------------------------------------
{
    shells_.push_back(shell);
    if (shells_.size() == 1)
    {
	toptol_ = shell->getTolerances();
    }
    else
    {
	double gap=0.0, neighbour=0.0, kink=0.0, bend=0.0;
	tpTolerances toptol = shell->getTolerances();
	gap = std::max(toptol_.gap, toptol.gap);
	neighbour = std::max(toptol_.neighbour, toptol.neighbour);
	kink = std::max(toptol_.kink, toptol.kink);
	bend = std::max(toptol_.bend, toptol.bend);

	toptol_ = tpTolerances(gap, neighbour, kink, bend);
    }
    addBodyPointers();
}

//---------------------------------------------------------------------------
BoundingBox Body::boundingBox() const
//---------------------------------------------------------------------------
{
  if (shells_.size() > 0)
    return shells_[0]->boundingBox();
  else
    {
      BoundingBox dummy;
      return dummy;
    }
}

//---------------------------------------------------------------------------
vector<shared_ptr<Vertex> > Body:: vertices() const
//---------------------------------------------------------------------------
{
  // Collect all vertices
 vector<shared_ptr<Vertex> > vertices;
 std::set<shared_ptr<Vertex> > all_vertices;  // All vertices in the model represented once

  for (size_t ki=0; ki<shells_.size(); ++ki)
    {
      vector<shared_ptr<Vertex> > curr_vertices;
      shells_[ki]->getAllVertices(curr_vertices);
      all_vertices.insert(curr_vertices.begin(), curr_vertices.end());
    }

  vertices.insert(vertices.end(), all_vertices.begin(), all_vertices.end());
  return vertices;
}


//---------------------------------------------------------------------------
bool Body::areNeighbours(Body *other, shared_ptr<ftSurface>& bd_face1,
			 shared_ptr<ftSurface>& bd_face2, int adj_idx) const
//---------------------------------------------------------------------------
{
  for (size_t ki=0; ki<shells_.size(); ++ki)
    {
      int nmb_faces = shells_[ki]->nmbEntities();
      for (int kj=0; kj<nmb_faces; ++kj)
	{
	  shared_ptr<ftSurface> curr = shells_[ki]->getFace(kj);
	  if (!curr->twin())
	    continue;

	  for (size_t kr=0; kr<other->shells_.size(); ++kr)
	    {
	      int nmb_faces2 = other->shells_[kr]->nmbEntities();
	      for (int kh=0; kh<nmb_faces2; ++kh)
		{
		  shared_ptr<ftSurface> curr2 = other->shells_[kr]->getFace(kh);
		  if (!curr2->twin())
		    continue;

		  if (curr->twin() == curr2.get() && curr2->twin() == curr.get())
		    {
		      // Decrement local adjecency index, return when desired index is found
		      if (adj_idx-- > 0)
			continue;
		      bd_face1 = curr;
		      bd_face2 = curr2;
		      return true;
		    }
		}
	    }
	}
    }
  return false;
}

//---------------------------------------------------------------------------
  bool Body::isInside(const Point& pnt)
//---------------------------------------------------------------------------
{
  double tol = 1.0e-10;

  // Fetch the midpoint of the bounding box
  BoundingBox box = boundingBox();
  Point mid = 0.5*(box.low() + box.high());
  if (mid.dist(pnt) <= toptol_.neighbour)
    mid = box.low();

  // Make a curve through the input point and this midpoint
  Point vec = pnt - mid;
  vec.normalize();
  Point vec2 = box.high() - box.low();
  double len = vec2.length();

  Point start = pnt - len*vec;
  Point end = pnt + len*vec;
  shared_ptr<SplineCurve> crv = 
    shared_ptr<SplineCurve>(new SplineCurve(start, -len, end, len));

  // Intersect all boundary shells with the curve
  vector<bool> segment;
  vector<pair<ftPoint, double> > int_pts;
  size_t ki, kj;
  for (ki=0; ki<shells_.size(); ++ki)
    {
      vector<bool> seg0;
      vector<pair<ftPoint, double> > int_pts0 = 
	shells_[ki]->intersect(crv, seg0);
      int_pts.insert(int_pts.end(), int_pts0.begin(), int_pts0.end());
      segment.insert(segment.end(), seg0.begin(), seg0.end());
    }

  // Remove touch points
  for (ki=0; ki<int_pts.size();)
    {
      if (segment[ki])
	{
	  int_pts.erase(int_pts.begin()+ki);
	  segment.erase(segment.begin()+ki);
	}
      else
	ki++;
    }

  // Remove duplicates
  for (ki=0; ki<int_pts.size(); ++ki)
    for (kj=ki+1; kj<int_pts.size(); )
      {
	if (fabs(int_pts[ki].second - int_pts[kj].second) < tol)
	  int_pts.erase(int_pts.begin()+kj);
	else
	  kj++;
      }


  // Count number of intersection points on either side of the
  // input point
  int nmb1, nmb2;
  for (nmb1=0, nmb2=0, ki=0; ki<int_pts.size(); ++ki)
    {
      if (int_pts[ki].second < 0.0)
	nmb1++;
      if (int_pts[ki].second > 0.0)
	nmb2++;
    }

  if (nmb1+nmb2 < (int)int_pts.size())
    return true;  // On boundary

  if (nmb1 % 2 == 0 || nmb2 % 2 == 0)
    return false;  // Outside
  else
    return true;
}

//---------------------------------------------------------------------------
  void Body::addBodyPointers()
//---------------------------------------------------------------------------
{
  for (size_t ki=0; ki<shells_.size(); ++ki)
    {
      int nmb = shells_[ki]->nmbEntities();
      for (int kj=0; kj<nmb; ++kj)
	{
	  shared_ptr<ftSurface> face = shells_[ki]->getFace(kj);
	  face->setBody(this);
	}
    }
}

//---------------------------------------------------------------------------
  shared_ptr<SurfaceModel> Body::getShell(ftSurface* face) const
//---------------------------------------------------------------------------
{
  shared_ptr<SurfaceModel> model;
  for (size_t ki=0; ki<shells_.size(); ++ki)
    {
      int nmb = shells_[ki]->nmbEntities();
      int kj;
      for (kj=0; kj<nmb; ++kj)
	{
	  shared_ptr<ftSurface> curr = shells_[ki]->getFace(kj);
	  if (curr.get() == face)
	    {
	      model = shells_[ki];
	      break;
	    }
	}
      if (kj < nmb)
	break;
    }
  return model;
}

} // namespace Go
