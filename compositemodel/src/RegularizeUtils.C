#define DEBUG_REG
#include "GoTools/compositemodel/RegularizeUtils.h"
#include "GoTools/geometry/BoundedUtils.h"
#include <fstream>

using namespace Go;
using std::vector;
using std::pair;
using std::make_pair;

//==========================================================================
vector<shared_ptr<ftSurface> > 
RegularizeUtils::divideVertex(shared_ptr<ftSurface> face,
			      shared_ptr<Vertex> vx, 
			      vector<shared_ptr<Vertex> > cand_vx,
			      ftEdge* cand_edge,
			      double epsge, double tol2, double angtol,
			      vector<shared_ptr<Vertex> > non_corner,
			      const Point& centre, const Point& axis)
//==========================================================================
{
#ifdef DEBUG_REG
  if (cand_vx.size() > 0)
{
  std::ofstream ofvx("cand_vx2.g2");
  ofvx << "400 1 0 4 155 100 0 255" << std::endl;
  ofvx << cand_vx.size() << std::endl;
  for (size_t kj=0; kj<cand_vx.size(); ++kj)
    ofvx << cand_vx[kj]->getVertexPoint() << std::endl;
}
#endif
  
  shared_ptr<ParamSurface> surf = face->surface();
  RectDomain dom = surf->containingDomain();
  Point vx_point = vx->getVertexPoint();

  // Statistics
  double parfrac = getMaxParFrac(face);

  // Get the plane with which to divide the current face to get subdivision
  // information
  Point pnt;
  Point normal;
  if (centre.dimension() > 0)
    {
      pnt = centre;
      normal = (vx_point - centre).cross(axis);
    }
  else
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

  // Fetch adjacent vertices
  vector<shared_ptr<Vertex> > next_vx(vx_edg.size());
    for (kr=0; kr<vx_edg.size(); ++kr)
      next_vx[kr] = vx_edg[kr]->getOtherVertex(vx.get());

  int close_idx;
  double close_dist;
  Point close_par;
  getClosestBoundaryPar(face, vx, vx_cvs, vx_point, epsge, 
			close_idx, close_dist, close_par);
#ifdef DEBUG_REG
  std::ofstream ofvx("closest_vx2.g2");
  ofvx << "400 1 0 4 0 100 155 255" << std::endl;
  ofvx << "1" << std::endl;
  ofvx << face->point(close_par[0],close_par[1]) << std::endl;
#endif

  // Traverse all candidate vertices and check if any is feasible for 
  // division
  // Compute distance between vertices and found plane
  // Select vertex with minimum distance
  int min_idx = -1;
  double min_frac = MAXDOUBLE;
  double max_frac = 0.0;
  double level_frac = std::max(std::min(0.5, 100.0*parfrac), 0.1);
  Point vx_par = vx->getFacePar(face.get());
  Point edge_par;
  double level_ang = M_PI/3; // M_PI/2.0; // M_PI/4.0; //M_PI/6.0;
  double min_ang = 1.0e8;
  double fac = 0.5; // 0.2;
  double fac2 = 2.0;
  double fac3 = 0.05;
  double min_dist = MAXDOUBLE;
  double tol = 1.0e-4;
  double curr_rad_dist = MAXDOUBLE;

  Point curr_vx_par1, curr_vx_par;
  double d1, d2, rad_dist;
  vector<shared_ptr<CurveOnSurface> > trim_segments;
  shared_ptr<BoundedSurface> bd_sf;
  if (centre.dimension() > 0)
    d1 = vx_point.dist(centre);
  while (true) 
    {
      min_idx = -1;
      min_frac = MAXDOUBLE;
      max_frac = 0.0;
      min_ang = 1.0e8;
      min_dist = MAXDOUBLE;
      curr_rad_dist = MAXDOUBLE;
      for (int ki=0; ki<(int)cand_vx.size(); ++ki)
	{
	  Point curr_vx_par2 = cand_vx[ki]->getFacePar(face.get());
	  Point vec = cand_vx[ki]->getVertexPoint() - vx_point;
	  double dist = vec.length();
	  double dist1 = fabs(vx_par[0]-curr_vx_par2[0]);
	  double dist2 = fabs(vx_par[1]-curr_vx_par2[1]);
	  dist1 /= (dom.umax() - dom.umin());
	  dist2 /= (dom.vmax() - dom.vmin());
	  double frac = std::min(dist1, dist2);
	  double ang = vec.angle(normal);
	  ang = fabs(0.5*M_PI - ang);
	  if (centre.dimension() > 0)
	    {
	      d2 = cand_vx[ki]->getVertexPoint().dist(centre);
	      rad_dist = fabs(d1 - d2);
	    }
	  else 
	    rad_dist = MAXDOUBLE;

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
	  Point vec1 = (curr_vx_par.dimension() == 0) ? Point(0.0, 0.0) : curr_vx_par - vx_par;
	  Point vec2 = curr_vx_par2 - vx_par;
	  double par_ang = vec1.angle(vec2);
	  double par_limit = 0.05*M_PI;
      
	  if (curr_rad_dist < epsge)
	    {
	      if (rad_dist < curr_rad_dist)
		{
		  curr_rad_dist = rad_dist;
		  min_ang = ang;
		  min_dist = dist;
		  min_frac = frac;
		  min_idx = ki;
		  curr_vx_par = curr_vx_par2;
		}
	    }
	  else
	    {
	      if (fabs(frac-min_frac) < tol && fabs(ang-min_ang) < tol &&
		  dist < min_dist)
		{
		  curr_rad_dist = rad_dist;
		  min_ang = ang;
		  min_dist = dist;
		  min_frac = frac;
		  min_idx = ki;
		  curr_vx_par = curr_vx_par2;
		}
	      else if ((frac < 0.9*min_frac && ang < level_ang && 
			dist < fac2*min_dist) || dist < fac*min_dist)
		{
		  curr_rad_dist = rad_dist;
		  min_ang = ang;
		  min_dist = dist;
		  min_frac = frac;
		  min_idx = ki;
		  curr_vx_par = curr_vx_par2;
		}
	      else if (par_ang < par_limit && dist < min_dist)
		{
		  curr_rad_dist = rad_dist;
		  min_ang = ang;
		  min_dist = dist;
		  min_frac = frac;
		  min_idx = ki;
		  curr_vx_par = curr_vx_par2;
		}
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
	      curr_edge->closestPoint(vx_point, clo_par, clo_pt, clo_dist, &seed);
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

      if (min_idx < 0)
	break;  // No vertex is choosen

      double ang = 0.0;
      if (min_idx >= 0)
	{
#ifdef DEBUG_REG
	  std::ofstream ofcurr("curr_cand_vx.g2");
	  ofcurr << "400 1 0 4 155 0 100 255" << std::endl;
	  ofcurr << 1 << std::endl;
	  ofcurr << cand_vx[min_idx]->getVertexPoint() << std::endl;
#endif
	  
	  Point vec = cand_vx[min_idx]->getVertexPoint() - vx_point;
	  ang = vec.angle(normal);
	}

      if (min_idx >= 0 && curr_rad_dist < epsge)
	{
	  // Perform cylinder intersection
	  double cyl_rad = 0.5*(d1 + d2);
	  trim_segments = BoundedUtils::getCylinderIntersections(surf, centre, 
								 axis, cyl_rad,
								 epsge, bd_sf);
      
	  // Remove intersections not connected with the initial point
	  for (kr=0; kr<trim_segments.size(); )
	    {
	      Point pos1 = trim_segments[kr]->ParamCurve::point(trim_segments[kr]->startparam());
	      Point pos2 = trim_segments[kr]->ParamCurve::point(trim_segments[kr]->endparam());
	      if (pos1.dist(vx_point) > epsge && pos2.dist(vx_point) > epsge)
		trim_segments.erase(trim_segments.begin()+kr);
	      else
		kr++;
	    }
	}

      // Avoid vertices which is inline with the neighbouring vertices
      // to the split vertex
      Point vec1 = cand_vx[min_idx]->getVertexPoint() - vx_point;
     for (kr=0; kr<next_vx.size(); ++kr)
	{
	  Point vec2 = next_vx[kr]->getVertexPoint() - vx_point;
	  double vx_ang = vec1.angle(vec2);
	  if (vx_ang < angtol)
	    break;
	}
     if (kr < next_vx.size())
       {
	 // This candidate vertex is not a good coice for splitting
	 ;
       }
     else if (trim_segments.size() == 0 && min_idx >= 0 && 
	  (fac*min_dist < close_dist || min_frac < fac3*max_frac) &&
	  (min_frac < level_frac*max_frac || 
	   fabs(0.5*M_PI - ang) < level_frac*level_ang))
	// if (min_idx >= 0 && (min_frac < fac*max_frac || 
	// 		       fabs(0.5*M_PI - ang) < level_ang))
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
	      if (pos1.dist(vx_point) > epsge && pos2.dist(vx_point) > epsge)
		trim_segments.erase(trim_segments.begin()+kr);
	      else
		kr++;
	    }
	}
      if (trim_segments.size() == 0)
	{
	  // The choosen vertex did not work. Remove it from the pool and try again
	  cand_vx.erase(cand_vx.begin()+min_idx);
	}
      else 
	break;
    }
  // if (false /*trim_segments.size() == 0 && cand_edge*/)
  //   {
  //     // Let the division curve end at a point of the given edge
  //     double tmid = 0.5*(cand_edge->tMin() + cand_edge->tMax());
  //     double clo_par, clo_dist;
  //     Point clo_pt;
  //     Point parval2;
  //     cand_edge->closestPoint(pnt, clo_par, clo_pt, clo_dist, &tmid);
  //     double p_len = cand_edge->tMax() - cand_edge->tMin();
  //     double lenfac = 0.1;
  //     if (clo_par - cand_edge->tMin() < lenfac*p_len ||
  // 	  cand_edge->tMax() - clo_par < lenfac*p_len)
  // 	parval2 = cand_edge->faceParameter(tmid);
  //     else
  // 	{
  // 	  Point face_seed = cand_edge->faceParameter(tmid);
  // 	  parval2 = cand_edge->faceParameter(clo_par, face_seed.begin()); 
  // 	}
  //     trim_segments = BoundedUtils::getTrimCrvsParam(surf, vx_par,
  // 						     parval2, epsge,
  // 						     bd_sf);
  //     // Check output
  //     Point par1, par2;
  //     for (kr=0; kr<trim_segments.size(); ++kr)
  // 	{
  // 	  par1 = 
  // 	    trim_segments[kr]->parameterCurve()->point(trim_segments[kr]->startparam());
  // 	  par2 = 
  // 	    trim_segments[kr]->parameterCurve()->point(trim_segments[kr]->endparam());
  // 	  if (par1.dist(parval2) < epsge || par2.dist(parval2) < epsge)
  // 	    break;
  // 	}
  //     if (kr == trim_segments.size() ||
  // 	  (trim_segments.size()>1 && 
  // 	   par1.dist(vx_par)>epsge && par2.dist(vx_par)>epsge))
  // 	trim_segments.clear();
	      
  //   }
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
	  trim_segments = BoundedUtils::getPlaneIntersections(surf, vx_point,
							      normal, epsge,
							      bd_sf);
	}

      // Remove intersections not connected with the initial point
      for (kr=0; kr<trim_segments.size(); )
	{
	  Point pos1 = trim_segments[kr]->ParamCurve::point(trim_segments[kr]->startparam());
	  Point pos2 = trim_segments[kr]->ParamCurve::point(trim_segments[kr]->endparam());
	  if (pos1.dist(vx_point) > epsge && pos2.dist(vx_point) > epsge)
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
  // Sort surfaces according to then number of loops
  for (size_t ki=0; ki<sub_sfs.size(); ++ki)
    {
      int nmb1 = sub_sfs[ki]->numberOfLoops();
      for (size_t kj=ki+1; kj<sub_sfs.size(); ++kj)
	{
	  int nmb2 = sub_sfs[kj]->numberOfLoops();
	  if (nmb2 > nmb1)
	    {
	      std::swap(sub_sfs[ki], sub_sfs[kj]);
	      std::swap(nmb1, nmb2);
	    }
	}
    }

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
  int nmb_corners = 0;
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
	nmb_corners++;
    }

  if (nmb_corners >= 2)
    return true;
  else
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

  if (vx.get() == last.get())
    return true;

  vector<ftEdge*> edges = vx->getFaceEdges(face.get());
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

      shared_ptr<Vertex> other = edges[ki]->getOtherVertex(vx.get());
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

//==========================================================================
bool
  RegularizeUtils::noExtension(shared_ptr<Vertex> vx, ftSurface* face,
			       shared_ptr<Vertex>& vx2, pair<Point, Point>& co_par1, 
			       pair<Point, Point>& co_par2, int& dir1, int& dir2,
			       double& val1, double& val2, double angtol, 
			       bool check_constant_curve)
//==========================================================================
{
  double ptol = 1.0e-8;
  size_t ki, kj, kr;

 // Fetch faces and emove faces from other bodies
  vector<ftSurface*> vx_faces = vx->faces();
  Body *bd0 = face->getBody();
  for (kj=0; kj<vx_faces.size();)
    {
      if (vx_faces[kj]->getBody() != bd0)
	vx_faces.erase(vx_faces.begin()+kj);
      else
	kj++;
    }
  // Check number of faces
  if (vx_faces.size() != 3)
    return false;

  // Get edges
  vector<ftEdge*> edges = vx->uniqueEdges();

  // Look for an edge connecting two T-vertices where 3 edges meet
  for (ki=0; ki<edges.size(); ++ki)
    {
      vx2 = edges[ki]->getOtherVertex(vx.get());
      vector<ftSurface*> vx_faces2 = vx2->faces();

      // Remove faces from other bodies
      for (kj=0; kj<vx_faces2.size();)
	{
	  if (vx_faces2[kj]->getBody() != bd0)
	    vx_faces2.erase(vx_faces2.begin()+kj);
	  else
	    kj++;
	}

      if (vx_faces2.size() != 3)
	continue;

      // Fetch the face which is not adjacent to the initial vertex
      for (kj=0; kj<vx_faces2.size(); ++kj)
	{
	  for (kr=0; kr<vx_faces.size(); ++kr)
	    if (vx_faces[kr] == vx_faces2[kj])
	      break;
	  if (kr == vx_faces.size())
	    break;  // Face found
	}
      if (kj == vx_faces2.size())
	continue;  // No such face is found

      // Check if vertex is a non-corner in this face
      vector<ftEdge*> edges2 = vx2->getFaceEdges(vx_faces2[kj]);
      if (edges2.size() != 2)
	continue;

      // Check angle
      double t1 = edges2[0]->parAtVertex(vx2.get());
      double t2 = edges2[1]->parAtVertex(vx2.get());
      Point tan1 = edges2[0]->tangent(t1);
      Point tan2 = edges2[1]->tangent(t2);
      if (tan1.angle(tan2) > angtol)
	continue;  // The faces meet in a corner
      // Check the continuity between the faces meeting in the
      // identified edge
      if (!edges[ki]->twin())
	continue;
      if (!edges[ki]->hasConnectivityInfo())
	continue;  // No continuity info
      shared_ptr<FaceConnectivity<ftEdgeBase> > info = 
	edges[ki]->getConnectivityInfo();

      int status = info->WorstStatus();
      if (status > 1)
	continue;  // Not G1 or almost G1

      if (check_constant_curve)
	{
	  // For the time being, check if the boundary between the
	  // faces are constant parameter curves in both parameter directions
	  shared_ptr<ParamCurve> cv1 = edges[ki]->geomCurve();
	  shared_ptr<ParamCurve> cv2 = edges[ki]->twin()->geomEdge()->geomCurve();
	  shared_ptr<CurveOnSurface> sf_cv1 = 
	    dynamic_pointer_cast<CurveOnSurface,ParamCurve>(cv1);
	  shared_ptr<CurveOnSurface> sf_cv2 = 
	    dynamic_pointer_cast<CurveOnSurface,ParamCurve>(cv2);
	  if (!sf_cv1.get() || 
	      !sf_cv1->isConstantCurve(ptol, dir1, val1))
	    continue;
	  if (!sf_cv2.get() || 
	      !sf_cv2->isConstantCurve(ptol, dir2, val2))
	    continue;

	  // Check if the next or previous edge follow the same constant boundary
	  int dir;
	  double val;
	  shared_ptr<ParamCurve> tmp_crv = edges[ki]->next()->geomEdge()->geomCurve();
	  shared_ptr<CurveOnSurface> sf_crv = 
	    dynamic_pointer_cast<CurveOnSurface,ParamCurve>(tmp_crv);
	  if (sf_crv.get())
	    sf_crv->isConstantCurve(ptol, dir, val);
	  else
	    dir = -1;
	  if (dir == dir1)
	    continue;
	  tmp_crv = edges[ki]->prev()->geomEdge()->geomCurve();
	  sf_crv = dynamic_pointer_cast<CurveOnSurface,ParamCurve>(tmp_crv);
	  if (sf_crv.get())
	    sf_crv->isConstantCurve(ptol, dir, val);
	  else
	    dir = -1;
	  if (dir == dir1)
	    continue;
	  tmp_crv = edges[ki]->twin()->next()->geomEdge()->geomCurve();
	  sf_crv = dynamic_pointer_cast<CurveOnSurface,ParamCurve>(tmp_crv);
	  if (sf_crv.get())
	    sf_crv->isConstantCurve(ptol, dir, val);
	  else
	    dir = -1;
	  if (dir == dir2)
	    continue;
	  tmp_crv = edges[ki]->twin()->prev()->geomEdge()->geomCurve();
	  sf_crv = dynamic_pointer_cast<CurveOnSurface,ParamCurve>(tmp_crv);
	  if (sf_crv.get())
	    sf_crv->isConstantCurve(ptol, dir, val);
	  else
	    dir = -1;
	  if (dir == dir2)
	    continue;

	  Point pt1 = cv1->point(edges[ki]->tMin());
	  Point pt2 = cv2->point(edges[ki]->twin()->tMin());
	  Point pt3 = cv2->point(edges[ki]->twin()->tMax());
	  if (pt1.dist(pt2) < pt1.dist(pt3))
	    {
	      co_par1 = make_pair(sf_cv1->parameterCurve()->point(edges[ki]->tMin()),
				  sf_cv2->parameterCurve()->point(edges[ki]->twin()->tMin()));
	      co_par2 = make_pair(sf_cv1->parameterCurve()->point(edges[ki]->tMax()),
				  sf_cv2->parameterCurve()->point(edges[ki]->twin()->tMax()));
	    }
	  else
	    {
	      co_par1 = make_pair(sf_cv1->parameterCurve()->point(edges[ki]->tMin()),
				  sf_cv2->parameterCurve()->point(edges[ki]->twin()->tMax()));
	      co_par2 = make_pair(sf_cv1->parameterCurve()->point(edges[ki]->tMax()),
				  sf_cv2->parameterCurve()->point(edges[ki]->twin()->tMin()));
	    }
	  dir1--;
	  dir2--;

	  return true;   // A possible insignificant splitting curve is found
	}
      else 
	return true;
    }

  return false;  // The vertex has an extension beyon the first adjacent face
}

//==========================================================================
double
RegularizeUtils::getMaxParFrac(shared_ptr<ftSurface> face)
//==========================================================================
{
  // Fetch all edges
  vector<shared_ptr<ftEdge> > edges = face->getAllEdges();

  double maxparfrac = 0.0;
  for (size_t ki=0; ki<edges.size(); ++ki)
    {
      Point par1 = edges[ki]->faceParameter(edges[ki]->tMin());
      Point par2 = edges[ki]->faceParameter(edges[ki]->tMax());
      double parfrac = std::min(fabs(par1[0]-par2[0]),fabs(par1[1]-par2[1]))/
	std::max(fabs(par1[0]-par2[0]),fabs(par1[1]-par2[1]));
      maxparfrac = std::max(maxparfrac, parfrac);
    }
  return maxparfrac;
}

 
