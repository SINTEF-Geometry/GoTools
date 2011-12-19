//#define DEBUG_REG
#include "GoTools/compositemodel/RegularizeUtils.h"
#include "GoTools/geometry/BoundedUtils.h"
#include <fstream>

using namespace Go;
using std::vector;

//==========================================================================
vector<shared_ptr<ftSurface> > 
RegularizeUtils::divideVertex(shared_ptr<ftSurface> face,
			      shared_ptr<Vertex> vx, 
			      vector<shared_ptr<Vertex> > cand_vx,
			      ftEdge* cand_edge,
			      double epsge, double tol2, double angtol,
			      vector<shared_ptr<Vertex> > non_corner)
//==========================================================================
{
  shared_ptr<ParamSurface> surf = face->surface();
  RectDomain dom = surf->containingDomain();

  // Get the plane with which to divide the current face to get subdivision
  // information
  Point pnt;
  Point normal;
  getDivisionPlane(face, vx, epsge, pnt, normal);

  // Fetch boundary curve information
  size_t kr, kh;
  vector<ftEdge*> vx_edg = vx->getFaceEdges(face.get());
  vector<shared_ptr<ParamCurve> > vx_cvs;
  for (kr=0; kr<vx_edg.size(); ++kr)
    {
      shared_ptr<ParamCurve> tmp = vx_edg[kr]->geomCurve();
      for (kh=0; kh<vx_cvs.size(); ++kh)
	if (vx_cvs[kh].get() == tmp.get())
	  break;
      if (kh == vx_cvs.size())
	vx_cvs.push_back(tmp);
    }

  int close_idx;
  double close_dist;
  Point close_par;
  getClosestBoundaryPar(face, vx, vx_cvs, pnt, epsge, 
			close_idx, close_dist, close_par);

  // Traverse all candidate vertices and check if any is feasible for 
  // division
  // Compute distance between vertices and found plane
  // Select vertex with minimum distance
  int min_idx = -1;
  double min_frac = 1.0e8;
  double max_frac = 0.0;
  Point vx_par = vx->getFacePar(face.get());
  Point edge_par;
  double level_ang = M_PI/3.0; // M_PI/4.0; //M_PI/6.0;
  double min_ang = 1.0e8;
  double fac = 0.2;
  double fac2 = 2.0;
  double min_dist = 1.0e8;
  double tol = 1.0e-4;

  Point curr_vx_par1, curr_vx_par;
  for (int ki=0; ki<(int)cand_vx.size(); ++ki)
    {
      Point curr_vx_par2 = cand_vx[ki]->getFacePar(face.get());
      Point vec = cand_vx[ki]->getVertexPoint() - vx->getVertexPoint();
      double dist = vec.length();
      double dist1 = fabs(vx_par[0]-curr_vx_par2[0]);
      double dist2 = fabs(vx_par[1]-curr_vx_par2[1]);
      dist1 /= (dom.umax() - dom.umin());
      dist2 /= (dom.vmax() - dom.vmin());
      double frac = std::min(dist1, dist2);
      double ang = vec.angle(normal);
      ang = fabs(0.5*M_PI - ang);

     // Check if the vertex is associated the same underlying curve
      // as the initial vertex. In that case, it is not a candidate
      // for split
      vector<ftEdge*> vx_edg2 = cand_vx[ki]->getFaceEdges(face.get());
      for (kr=0; kr<vx_edg2.size(); ++kr)
	{
	  shared_ptr<ParamCurve> cv = vx_edg2[kr]->geomCurve();
	  for (kh=0; kh<vx_cvs.size(); ++kh)
	    if (cv.get() == vx_cvs[kh].get())
	      break;
	  if (kh < vx_cvs.size())
	    break;
	}
      if (kr < vx_edg2.size())
	{
	  curr_vx_par = curr_vx_par1 = curr_vx_par2;
	  continue;  // Vertex not allowed for split
	}

      // Check for corners between the current and the candidate vertex
      if (!cornerInShortestPath(vx, cand_vx[ki], face, angtol))
	  continue;

      // Compute the angle between the vector from the split vertex to the
      // previous choice of destination vertex and the corresponding vector
      // for the current on in the parameter domain. If these vectors are almost
      // parallel, the closest candidate will be choosen
      Point vec1 = (ki == 0) ? Point(0.0, 0.0) : curr_vx_par - vx_par;
      Point vec2 = curr_vx_par2 - vx_par;
      double par_ang = vec1.angle(vec2);
      double par_limit = 0.05*M_PI;
      
      if (fabs(frac-min_frac) < tol && fabs(ang-min_ang) < tol &&
	  dist < min_dist)
	{
	  min_ang = ang;
	  min_dist = dist;
	  min_frac = frac;
	  min_idx = ki;
	  curr_vx_par = curr_vx_par2;
	}
      else if ((frac < min_frac && ang < level_ang && 
		(ang < fac2*min_ang ||  dist < fac2*min_dist)) ||
	       dist < fac*min_dist)
	{
	  min_ang = ang;
	  min_dist = dist;
	  min_frac = frac;
	  min_idx = ki;
	  curr_vx_par = curr_vx_par2;
	}
      else if (par_ang < par_limit && dist < min_dist)
	{
	  min_ang = ang;
	  min_dist = dist;
	  min_frac = frac;
	  min_idx = ki;
	  curr_vx_par = curr_vx_par2;
	}
      max_frac = std::max(max_frac, std::max(dist1,dist2));

      
      // Check edge between vertices
      ftEdge *curr_edge = NULL;
      if (ki > 0)
	cand_vx[ki-1]->getCommonEdgeInFace(cand_vx[ki].get(),
					       face.get());

      if (curr_edge)
	{
	  // Perform closest point to vertex
	  double seed = 0.5*(curr_edge->tMin() + curr_edge->tMax());
	  double clo_par, clo_dist;
	  Point clo_pt;
	  curr_edge->closestPoint(pnt, clo_par, clo_pt, clo_dist, &seed);
	  Point face_seed = 0.5*(curr_vx_par1 + curr_vx_par2);
	  Point curr_e_par = curr_edge->faceParameter(clo_par, 
						      face_seed.begin());

	  dist1 = fabs(vx_par[0]-curr_e_par[0]);
	  dist2 = fabs(vx_par[1]-curr_e_par[1]); 
	  dist1 /= (dom.umax() - dom.umin());
	  dist2 /= (dom.vmax() - dom.vmin());
	  frac = std::min(dist1, dist2);
	  double mfac = 10.0;
	  if (frac < mfac*min_frac)
	    {
	      min_frac = frac;
	      min_idx = -1;
	      edge_par = curr_e_par;
	    }
	  max_frac = std::max(max_frac, std::max(dist1,dist2));
	}
      curr_vx_par = curr_vx_par1 = curr_vx_par2;
    }

  vector<shared_ptr<CurveOnSurface> > trim_segments;
  shared_ptr<BoundedSurface> bd_sf;
  double ang = 0.0;
  if (min_idx >= 0)
    {
      Point vec = cand_vx[min_idx]->getVertexPoint() - vx->getVertexPoint();
      ang = vec.angle(normal);
    }
  if (min_idx >= 0 && (/*min_frac < fac*max_frac && */
		       fabs(0.5*M_PI - ang) < level_ang))
    {
      // Find division curve between vertices
      Point parval2 = cand_vx[min_idx]->getFacePar(face.get());
      trim_segments = BoundedUtils::getTrimCrvsParam(surf, vx_par,
						     parval2, epsge,
						     bd_sf);

      // Check output
      Point par1, par2;
      for (kr=0; kr<trim_segments.size(); ++kr)
	{
	  par1 = 
	    trim_segments[kr]->parameterCurve()->point(trim_segments[kr]->startparam());
	  par2 = 
	    trim_segments[kr]->parameterCurve()->point(trim_segments[kr]->endparam());
	  if (par1.dist(parval2) < epsge || par2.dist(parval2) < epsge)
	    break;
	}
      if (kr == trim_segments.size() ||
	  (trim_segments.size()>1 && 
	   par1.dist(vx_par)>epsge && par2.dist(vx_par)>epsge))
	trim_segments.clear();
	      
      // Remove intersections not connected with the initial point
      for (kr=0; kr<trim_segments.size(); )
	{
	  Point pos1 = trim_segments[kr]->ParamCurve::point(trim_segments[kr]->startparam());
	  Point pos2 = trim_segments[kr]->ParamCurve::point(trim_segments[kr]->endparam());
	  if (pos1.dist(pnt) > epsge && pos2.dist(pnt) > epsge)
	    trim_segments.erase(trim_segments.begin()+kr);
	  else
	    kr++;
	}
    }
  if (false /*trim_segments.size() == 0 && cand_edge*/)
    {
      // Let the division curve end at a point of the given edge
      double tmid = 0.5*(cand_edge->tMin() + cand_edge->tMax());
      double clo_par, clo_dist;
      Point clo_pt;
      Point parval2;
      cand_edge->closestPoint(pnt, clo_par, clo_pt, clo_dist, &tmid);
      double p_len = cand_edge->tMax() - cand_edge->tMin();
      double lenfac = 0.1;
      if (clo_par - cand_edge->tMin() < lenfac*p_len ||
	  cand_edge->tMax() - clo_par < lenfac*p_len)
	parval2 = cand_edge->faceParameter(tmid);
      else
	{
	  Point face_seed = cand_edge->faceParameter(tmid);
	  parval2 = cand_edge->faceParameter(clo_par, face_seed.begin()); 
	}
      trim_segments = BoundedUtils::getTrimCrvsParam(surf, vx_par,
						     parval2, epsge,
						     bd_sf);
      // Check output
      Point par1, par2;
      for (kr=0; kr<trim_segments.size(); ++kr)
	{
	  par1 = 
	    trim_segments[kr]->parameterCurve()->point(trim_segments[kr]->startparam());
	  par2 = 
	    trim_segments[kr]->parameterCurve()->point(trim_segments[kr]->endparam());
	  if (par1.dist(parval2) < epsge || par2.dist(parval2) < epsge)
	    break;
	}
      if (kr == trim_segments.size() ||
	  (trim_segments.size()>1 && 
	   par1.dist(vx_par)>epsge && par2.dist(vx_par)>epsge))
	trim_segments.clear();
	      
    }
  if (trim_segments.size() == 0 && min_frac < fac*max_frac &&
	   edge_par.dimension() == 2)
    {
      // Let the division curve end at an edge closest point
      trim_segments = BoundedUtils::getTrimCrvsParam(surf, vx_par,
						     edge_par, epsge,
						     bd_sf);
      // Check output
      Point par1, par2;
      for (kr=0; kr<trim_segments.size(); ++kr)
	{
	  par1 = 
	    trim_segments[kr]->parameterCurve()->point(trim_segments[kr]->startparam());
	   par2 = 
	    trim_segments[kr]->parameterCurve()->point(trim_segments[kr]->endparam());
	  if (par1.dist(edge_par) < epsge || par2.dist(edge_par) < epsge)
	    break;
	}
      if (kr == trim_segments.size() ||
	  (trim_segments.size()>1 && 
	   par1.dist(vx_par)>epsge && par2.dist(vx_par)>epsge))
	trim_segments.clear();
	      
    }
  if (trim_segments.size() == 0)
    {
      // Check if a constant parameter curve is a feasible 
      // split curve. First evaluate constant parameter tangents
      vector<Point> pts(3);
      surf->point(pts, vx_par[0], vx_par[1], 1);
      double ang1 = pts[1].angle(normal);
      double ang2 = pts[2].angle(normal);
      ang1 = fabs(0.5*M_PI - ang1);
      ang2 = fabs(0.5*M_PI - ang2);
      double d1 = 0.0, d2 = 0.0;
      if (close_idx >= 0)
	{
	  d1 = fabs(vx_par[1] - close_par[1]);
	  d2 = fabs(vx_par[0] - close_par[0]);
	}
      if ((ang1 < 0.25*ang2 ||
	   (close_idx>=0 && (d1 < 0.01*d2 || d1 < epsge))) 
	  && ang1 < level_ang)
	{
	  Point parval1(2), parval2(2);
	  parval1[1] = parval2[1] = vx_par[1];
	  parval1[0] = dom.umin();
	  parval2[0] = dom.umax();
	  trim_segments = BoundedUtils::getTrimCrvsParam(surf, parval1,
							 parval2, epsge,
							 bd_sf);
 	}
      else if ((ang2 < 0.25*ang1 ||
		(close_idx>=0 && (d2 < 0.01*d1 || d2 < epsge))) 
	       && ang2 < level_ang)
	{
	  Point parval1(2), parval2(2);
	  parval1[0] = parval2[0] = vx_par[0];
	  parval1[1] = dom.vmin();
	  parval2[1] = dom.vmax();
	  trim_segments = BoundedUtils::getTrimCrvsParam(surf, parval1,
							 parval2, epsge,
							 bd_sf);
 	}
      else
	{
	  // Find intersections between the face and this plane
	  trim_segments = BoundedUtils::getPlaneIntersections(surf, pnt,
							      normal, epsge,
							      bd_sf);
	}

      // Remove intersections not connected with the initial point
      for (kr=0; kr<trim_segments.size(); )
	{
	  Point pos1 = trim_segments[kr]->ParamCurve::point(trim_segments[kr]->startparam());
	  Point pos2 = trim_segments[kr]->ParamCurve::point(trim_segments[kr]->endparam());
	  if (pos1.dist(pnt) > epsge && pos2.dist(pnt) > epsge)
	    trim_segments.erase(trim_segments.begin()+kr);
	  else
	    kr++;
	}

    }

  // Check feasability of intersections
  //MESSAGE("Some code is missing to check feasability of intersection");

#ifdef DEBUG_REG
  std::ofstream out_file("split_segments.g2");
  for (size_t kj=0; kj<trim_segments.size(); ++kj)
    {
      shared_ptr<ParamCurve> cv = trim_segments[kj]->spaceCurve();
      cv->writeStandardHeader(out_file);
      cv->write(out_file);
    }
#endif

  // Define faces
  vector<shared_ptr<BoundedSurface> > sub_sfs =
    BoundedUtils::splitWithTrimSegments(bd_sf, trim_segments, epsge);

#ifdef DEBUG_REG
  std::ofstream of("split_surf.g2");
  for (kr=0; kr<sub_sfs.size(); ++kr)
    {
      sub_sfs[kr]->writeStandardHeader(of);
      sub_sfs[kr]->write(of);
    }
#endif

  // Create faces
  vector<shared_ptr<ftSurface> > faces = createFaces(sub_sfs, face,
						     epsge, tol2, angtol,
						     non_corner);
  return faces;
}


//==========================================================================
vector<shared_ptr<ftSurface> > 
RegularizeUtils::createFaces(vector<shared_ptr<BoundedSurface> >& sub_sfs,
			     shared_ptr<ftSurface> face,
			     double epsge, double tol2, double angtol,
			     vector<shared_ptr<Vertex> > non_corner)
//==========================================================================
{
  // Create faces
  vector<shared_ptr<ftSurface> > faces;
  faces.reserve(sub_sfs.size());
  for (size_t kj=0; kj<sub_sfs.size(); ++kj)
    {
      shared_ptr<ftSurface> curr =
	shared_ptr<ftSurface>(new ftSurface(sub_sfs[kj], -1));
      curr->setBody(face->getBody());
      (void)curr->createInitialEdges(epsge, angtol/*, true*/);

      // Check if any non-corner vertices belongs to the new face
      double frac = 0.001;
      vector<shared_ptr<ftEdge> > edges = curr->getAllEdges();
      for (size_t kr=0; kr<non_corner.size(); kr++)
	{
	  Point vx_pt = non_corner[kr]->getVertexPoint();
	  for (size_t kh=0; kh<edges.size(); ++kh)
	    {
	      double t1 = edges[kh]->tMin();
	      double t2 = edges[kh]->tMax();
	      double par, dist;
	      Point close;
	      edges[kh]->closestPoint(vx_pt, par, close, dist);
	      if (dist < tol2 && 
		  (par-t1 > frac*(t2-t1) && t2-par > frac*(t2-t1)))
		{
		  shared_ptr<ftEdgeBase> new_edge = edges[kh]->split2(par);
		  shared_ptr<ftEdge> curr_edge = 
		    dynamic_pointer_cast<ftEdge,ftEdgeBase>(new_edge);
		  if (curr_edge.get())
		    {
		      shared_ptr<Vertex> vx = curr_edge->getVertex(true);
// 		      if (vx.get())
// 			vx->joinVertex(non_corner[kr]);
		      edges.push_back(curr_edge);
		    }
		  //curr->updateBoundaryLoops(new_edge);
		  break;
		}
	    }
	}
	      
      faces.push_back(curr);
    }

  return faces;
}

//==========================================================================
  void RegularizeUtils::getDivisionPlane(shared_ptr<ftSurface> face,
					 shared_ptr<Vertex> vx,
					 double epsge,
					 Point& pnt,
					 Point& normal)
//==========================================================================
{
  // Get parameter corresponding to the vertex in the given face
  Point param = vx->getFacePar(face.get());

  // Get the edges corresponding to this face meeting in the vertex
  vector<ftEdge*> edges =  vx->getFaceEdges(face.get());

  // Theoretically a face can have an equal number of edges meeting in 
  // a given vertex. Assume two, and use the two first edges
  ASSERT(edges.size() >= 2);

  // Define plane
  pnt = face->point(param[0], param[1]);
  Point norm = face->normal(param[0], param[1]);
  norm.normalize();
  double t1 = edges[0]->parAtVertex(vx.get());
  double t2 = edges[1]->parAtVertex(vx.get());
  Point tan1 = edges[0]->tangent(t1);
  Point tan2 = edges[1]->tangent(t2);
  tan1.normalize();
  tan2.normalize();
  Point vec = 0.5*(tan1 + tan2);
  if (vec.length() < epsge)
    vec = tan1;

  // Project the vector into the tangent plane of the face
  // normal = vec - (vec*norm)*norm;
  // if (normal.length() < epsge)
    // normal = vec%norm;

  normal = vec;
  normal.normalize();
}

//==========================================================================
void 
RegularizeUtils::getClosestBoundaryPar(shared_ptr<ftSurface> face,
				       shared_ptr<Vertex> vx,
				       vector<shared_ptr<ParamCurve> >& vx_cvs,
				       const Point& pnt,
				       double epsge,
				       int& close_idx, double& close_dist,
				       Point& close_par)
//==========================================================================
{
  close_idx = -1;
  close_dist = 1.0e8;

  // Fetch information about boundary curves
  size_t kr, kh;
  vector<shared_ptr<ftEdge> > all_edg = face->getAllEdges();
  vector<shared_ptr<ParamCurve> > cvs;
  for (kr=0; kr<all_edg.size(); ++kr)
    {
      shared_ptr<ParamCurve> tmp = all_edg[kr]->geomCurve();
      for (kh=0; kh<cvs.size(); ++kh)
	if (cvs[kh].get() == tmp.get())
	  break;
      if (kh == cvs.size())
	cvs.push_back(tmp);
    }
 
  
  int idx1 = -1, idx2 = -1;
  for (kr=0; kr<cvs.size(); ++kr)
    {
      for (kh=0; kh<vx_cvs.size(); ++kh)
	if (vx_cvs[kh].get() == cvs[kr].get())
	  {
	    if (idx1 < 0)
		idx1 = (int)kr;
	    else
		idx2 = (int)kr;
	  }
    }
  if (idx2 < 0)
    {
      if (idx1 == 0)
	{
	  cvs.erase(cvs.begin(), 
		    cvs.begin()+std::min((int)(cvs.size()),2));
	  if (cvs.size() > 0)
	    cvs.erase(cvs.end()-1);
	}
      else if (idx1 >= (int)cvs.size()-1)
	{
	  cvs.erase(cvs.begin()+idx1-1, cvs.begin()+idx1+1);
	  if (cvs.size() > 0)
	    cvs.erase(cvs.begin());
	}
      else
	cvs.erase(cvs.begin()+idx1-1, cvs.begin()+idx1+2);
    }
  else
    {
      cvs.erase(cvs.begin()+idx2);
      cvs.erase(cvs.begin()+idx1);
    }

  Point close_pnt;
  double close_t;
  shared_ptr<ParamSurface> surf = face->surface();
  for (kr=0; kr<cvs.size(); ++kr)
    {
      double dist, upar, vpar;
      Point close_pnt2;
      cvs[kr]->closestPoint(pnt, cvs[kr]->startparam(), cvs[kr]->endparam(),
			    close_t, close_pnt, dist);
      if (dist < close_dist)
	{
	  close_dist = dist;
	  close_idx = (int)kr;

	  surf->closestPoint(close_pnt, upar, vpar, close_pnt2, dist,
			     epsge);
	  close_par = Point(upar,vpar);
	}
    }

}


//==========================================================================
bool
  RegularizeUtils::cornerInShortestPath(shared_ptr<Vertex> vx1,
					shared_ptr<Vertex> vx2,
					shared_ptr<ftSurface> face,
					double angtol)
//==========================================================================
{
  // Find shortest path
  vector<ftEdge*> edg1 = vx1->getFaceEdges(face.get());
  
  vector<ftEdge*> shortest_path;
  size_t ki;
  for (ki=0; ki<edg1.size(); ++ki)
    {
      shared_ptr<Vertex> other_vx = edg1[ki]->getOtherVertex(vx1.get());
      vector<ftEdge*> path;
      bool found = getPath(edg1[ki], other_vx, vx2, face, path);
      if (found && (shortest_path.size() == 0 || path.size() < shortest_path.size()))
	shortest_path = path;
    }

  // Check for corners
  for (ki=1; ki<shortest_path.size(); ++ki)
    {
      shared_ptr<Vertex> common_vx = 
	shortest_path[ki-1]->getCommonVertex(shortest_path[ki]);
      double t1 = shortest_path[ki-1]->parAtVertex(common_vx.get());
      double t2 = shortest_path[ki]->parAtVertex(common_vx.get());
      Point tan1 = shortest_path[ki-1]->tangent(t1);
      Point tan2 = shortest_path[ki]->tangent(t2);
      double ang = tan1.angle(tan2);
      if (ang > angtol)
	return true;
    }

  return false;
}


//==========================================================================
bool
  RegularizeUtils::getPath(ftEdge* edg, shared_ptr<Vertex> vx,
			   shared_ptr<Vertex> last, shared_ptr<ftSurface> face,
			   vector<ftEdge*>& path)
//==========================================================================
{
  path.push_back(edg);

  shared_ptr<Vertex> other = edg->getOtherVertex(vx.get());
  if (other.get() == last.get())
    return true;

  vector<ftEdge*> edges = other->getFaceEdges(face.get());
  vector<ftEdge*> shortest_path;
  for (size_t ki=0; ki<edges.size(); ++ki)
    {
      if (edges[ki] == edg)
	continue;
      size_t kj;
      for (kj=0; kj<path.size(); ++kj)
	if (path[kj] == edges[ki])
	  break;
      if (kj < path.size())
	continue;

      vector<ftEdge*> curr_path = path;  
      bool found = getPath(edges[ki], other, last, face, curr_path);
      if (found && (shortest_path.size() == 0 || curr_path.size() < shortest_path.size()))
	shortest_path = curr_path;
    }
  
  if (shortest_path.size() > 0)
    {
      path = shortest_path;
      return true;
    }
  else
    return false;
}
