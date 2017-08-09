/*
 * Copyright (C) 1998, 2000-2007, 2010, 2011, 2012, 2013 SINTEF ICT,
 * Applied Mathematics, Norway.
 *
 * Contact information: E-mail: tor.dokken@sintef.no                      
 * SINTEF ICT, Department of Applied Mathematics,                         
 * P.O. Box 124 Blindern,                                                 
 * 0314 Oslo, Norway.                                                     
 *
 * This file is part of GoTools.
 *
 * GoTools is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Affero General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version. 
 *
 * GoTools is distributed in the hope that it will be useful,        
 * but WITHOUT ANY WARRANTY; without even the implied warranty of         
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the          
 * GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public
 * License along with GoTools. If not, see
 * <http://www.gnu.org/licenses/>.
 *
 * In accordance with Section 7(b) of the GNU Affero General Public
 * License, a covered work must retain the producer line in every data
 * file that is created or manipulated using GoTools.
 *
 * Other Usage
 * You can be released from the requirements of the license by purchasing
 * a commercial license. Buying such a license is mandatory as soon as you
 * develop commercial activities involving the GoTools library without
 * disclosing the source code of your own applications.
 *
 * This file may be used in accordance with the terms contained in a
 * written agreement between you and SINTEF ICT. 
 */

//#define DEBUG_REG

#include "GoTools/compositemodel/RegularizeUtils.h"
#include "GoTools/geometry/BoundedUtils.h"
#include "GoTools/geometry/HermiteInterpolator.h"
#include "GoTools/geometry/ElementaryCurve.h"
#include "GoTools/geometry/Circle.h"
#include "GoTools/geometry/SplineCurve.h"
#include "GoTools/creators/CoonsPatchGen.h"
#include "sislP.h"
#include "GoTools/geometry/SISLconversion.h"
#include <fstream>

using namespace Go;
using std::vector;
using std::pair;
using std::make_pair;

//==========================================================================
vector<shared_ptr<ftSurface> > 
RegularizeUtils::divideVertex(shared_ptr<ftSurface> face,
			      shared_ptr<Vertex> vx, 
			      vector<shared_ptr<Vertex> >& cand_vx,
			      ftEdge* cand_edge,
			      vector<shared_ptr<Vertex> >& prio_vx,
			      double epsge, double tol2, double angtol,
			      double bend,
			      vector<shared_ptr<Vertex> >& non_corner,
			      const Point& centre, const Point& axis,
			      bool strong)
//==========================================================================
{
  // Perform splitting
  vector<shared_ptr<CurveOnSurface> > trim_segments;
  shared_ptr<BoundedSurface> bd_sf;

  trim_segments = findVertexSplit(face, vx, cand_vx, cand_edge,
				  prio_vx, epsge, tol2, angtol, bend,
				  non_corner, centre, axis, bd_sf,
				  strong);
  // Define faces
  vector<shared_ptr<BoundedSurface> > sub_sfs =
    BoundedUtils::splitWithTrimSegments(bd_sf, trim_segments, epsge);

#ifdef DEBUG_REG
  std::ofstream of("split_surf.g2");
  for (size_t kr=0; kr<sub_sfs.size(); ++kr)
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
vector<shared_ptr<CurveOnSurface> > 
RegularizeUtils::findVertexSplit(shared_ptr<ftSurface> face,
				 shared_ptr<Vertex> vx, 
				 vector<shared_ptr<Vertex> >& cand_vx,
				 ftEdge* cand_edge,
				 vector<shared_ptr<Vertex> >& prio_vx,
				 double epsge, double tol2, double angtol,
				 double bend,
				 vector<shared_ptr<Vertex> >& non_corner,
				 const Point& centre, const Point& axis,
				 shared_ptr<BoundedSurface>& bd_sf,
				 bool strong)
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
  ofvx << "400 1 0 4 0 100 155 255" << std::endl;
  ofvx << prio_vx.size() << std::endl;
  for (size_t kj=0; kj<prio_vx.size(); ++kj)
    ofvx << prio_vx[kj]->getVertexPoint() << std::endl;
}
#endif
  
  shared_ptr<ParamSurface> surf = face->surface();
  RectDomain dom = surf->containingDomain();
  Point vx_point = vx->getVertexPoint();
  Point vx_par = vx->getFacePar(face.get());
  double level_ang = M_PI/3; // M_PI/2.0; // M_PI/4.0; //M_PI/6.0;
  Point centre2 = centre;
  Point axis2;

  // Fetch adjacent vertices
  vector<shared_ptr<Vertex> > next_vxs = vx->getNextVertex(face.get());

  // Get the plane with which to divide the current face to get subdivision
  // information
  Point pnt;
  Point normal;
  if (centre.dimension() > 0)
    {
      pnt = centre;
      normal = (vx_point - centre).cross(axis);
    }
  else if (axis.dimension() > 0)
    {
      pnt = vx_point;
      normal = axis;
    }
   else
    getDivisionPlane(face, vx, epsge, pnt, normal);

  // Fetch a vector in the given vertex pointing into the surface
  Point in_vec = getInVec(vx, face);

  // A patch with 5 sides requires a special treatment to avoid
  // a degenerate solution. Check the number of corners
  // The same applies to a patch with 4 corners where the given vertex
  // is not one of the corners
  vector<shared_ptr<Vertex> > corners = face->getCornerVertices(bend);
  ftEdge* opposite=NULL;  // Pointer to edge with opposite point if set
  double opposite_dist;
  double opposite_par;
  Point opposite_point;
  if (corners.size() == 4 || corners.size() == 5)
    opposite = getOppositeBoundaryPar(face, vx, corners, epsge, 
				      opposite_point, opposite_par, 
				      opposite_dist);
    
  if (opposite)
    {
#ifdef DEBUG_REG
      std::ofstream ofvx2("opposite_vx2.g2");
      ofvx2 << "400 1 0 4 0 100 155 255" << std::endl;
      ofvx2 << "1" << std::endl;
      ofvx2 << opposite_point << std::endl;
#endif
    }

   // Fetch boundary curve information
  vector<shared_ptr<ftEdge> > all_edg = face->getAllEdges();
  vector<shared_ptr<ParamCurve> > all_cvs;
  getSourceCvs(all_edg, all_cvs);

   // Fetch boundary curve information related to selected split vertex
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

  if (centre.dimension() == 0)
    {
      // Check for circular behaviour
      size_t ka;
      for (ka=0; ka<all_cvs.size(); ++ka)
	{
	  if (!all_cvs[ka].get())
	    continue;
	  if (all_cvs[ka]->instanceType() == Class_Circle)
	    {
	      // Check that the circle is not directly connected to
	      // the split vertex
	      size_t kb;
	      for (kb=0; kb<vx_edg.size(); ++kb)
		if (all_edg[ka].get() == vx_edg[kb])
		  break;

	      // @@@ VSK. There is a risk that the curve is the same, but
	      // the edge is different so the test above is probably too simple
	      
	      if (kb == vx_edg.size())
		break;
	    }
	}
      if (ka < all_cvs.size())
	{
	  shared_ptr<Circle> circ = 
	    dynamic_pointer_cast<Circle,ParamCurve>(all_cvs[ka]);
	  centre2 = circ->getCentre();
	  axis2 = circ->getNormal();
	  double axis_ang = axis2.angle(normal);
	  double level_ang = 0.25*M_PI;
	  if (std::min(axis_ang, fabs(M_PI-axis_ang)) < level_ang)
	    {
	      centre2.resize(0);
	      axis2.resize(0);
	    }
	}
    }

  int close_idx;
  double close_dist;
  Point close_par;
  getClosestBoundaryPar(face, vx, vx_cvs, vx_point, epsge, 
			close_idx, close_dist, close_par, 
			strong ? 0 : -1);
  Point close_pt = face->point(close_par[0], close_par[1]);

#ifdef DEBUG_REG
  std::ofstream ofvx3("closest_vx2.g2");
  ofvx3 << "400 1 0 4 0 100 155 255" << std::endl;
  ofvx3 << "1" << std::endl;
  ofvx3 << face->point(close_par[0],close_par[1]) << std::endl;
  if (cand_edge)
    {
      std::ofstream ofvx4("edge_vx2.g2");
      ofvx4 << "400 1 0 4 155 0 100 255" << std::endl;
      ofvx4 << "1" << std::endl;
      ofvx4 << cand_edge->point(0.5*(cand_edge->tMin()+cand_edge->tMax()));
     }
#endif

  // Traverse all candidate vertices and check if any is feasible for 
  // division
  // Compute distance between vertices and found plane
  // Select vertex with minimum distance
  vector<shared_ptr<CurveOnSurface> > trim_segments;
  while (prio_vx.size() > 0) 
    {
      size_t nmb_cand = prio_vx.size();
      double cyl_rad = -1.0;
      int min_idx = selectCandVx(face, vx, in_vec, prio_vx, dom, epsge, 
				 bend, centre2, normal, vx_cvs, close_dist, 
				 close_pt, cyl_rad, strong);
      if (min_idx < 0)
	{
	  if (prio_vx.size() < nmb_cand)
	    continue;  // No vertex is choosen, look for a new
	  else
	    break;     // No legal candidate vertex
	}

#ifdef DEBUG_REG
      if (min_idx >= 0)
	{
	  std::ofstream ofcurr0("curr_prio_vx.g2");
	  ofcurr0 << "400 1 0 4 155 0 100 255" << std::endl;
	  ofcurr0 << 1 << std::endl;
	  ofcurr0 << prio_vx[min_idx]->getVertexPoint() << std::endl;
	}
#endif

      if (min_idx >= 0 && cyl_rad > 0.0)
	{
	  // Perform cylinder intersection
	  trim_segments = BoundedUtils::getCylinderIntersections(surf, centre2, 
								 axis2, cyl_rad,
								 epsge, bd_sf);
      
	  // Remove intersections not connected with the initial point
	  // Remove also segments going through an adjacent vertex
	  Point dummy;
	  checkTrimSeg(trim_segments, next_vxs, vx_point, dummy, tol2 /*epsge*/);
	}
      if (trim_segments.size() == 0 && min_idx >= 0)
	{
	  // Check the feasability of a stright curve in the
	  // parameter domain
	  Point parval2 = prio_vx[min_idx]->getFacePar(face.get());
	  shared_ptr<ParamCurve> pcurve = checkStrightParCv(face, vx, prio_vx[min_idx], 
							    epsge);
	  if (pcurve.get())
	    {
	      trim_segments = BoundedUtils::getTrimCrvsPcrv(surf, pcurve, epsge,
							    bd_sf);
	    }
	  else
	    {
	      // Find division curve between vertices
	      trim_segments = BoundedUtils::getTrimCrvsParam(surf, vx_par,
							     parval2, epsge,
							     bd_sf);
	    }
#ifdef DEBUG_REG
	  std::ofstream ofcv("trim_seg0.g2");
	  for (size_t kj=0; kj<trim_segments.size(); ++kj)
	    {
	      shared_ptr<ParamCurve> cv = trim_segments[kj]->spaceCurve();
	      cv->writeStandardHeader(ofcv);
	      cv->write(ofcv);
	    }
#endif
	  // Check output
	  checkTrimSeg2(trim_segments, vx_par, parval2, epsge);
	      
	  // Remove intersections not connected with the initial point
	  // Remove also segments going through an adjacent vertex
	  if (trim_segments.size() > 0)
	    {
	      Point other_pt = prio_vx[min_idx]->getVertexPoint();
	      checkTrimSeg(trim_segments, next_vxs, vx_point, 
			   other_pt, tol2 /*epsge*/);
	    }
	}
      if (trim_segments.size() == 0)
	{
	  // The choosen vertex did not work. Remove it from the pool and try again
	  prio_vx.erase(prio_vx.begin()+min_idx);
	}
      else 
	break;
    }

  if (trim_segments.size() == 0)
    {
      while (cand_vx.size() > 0) 
	{
	  size_t nmb_cand = cand_vx.size();
	  double cyl_rad = -1.0;
	  int min_idx = selectCandVx(face, vx, in_vec, cand_vx, dom, epsge, 
				     bend, centre2, normal, vx_cvs, close_dist,
				      close_pt, cyl_rad, strong);
	  if (min_idx < 0)
	    {
	      if (cand_vx.size() < nmb_cand)
		continue;  // No vertex is choosen, look for a new
	      else
		break;     // No legal candidate vertex
	    }

#ifdef DEBUG_REG
	  if (min_idx >= 0)
	    {
	      std::ofstream ofcurr("curr_cand_vx.g2");
	      ofcurr << "400 1 0 4 155 0 100 255" << std::endl;
	      ofcurr << 1 << std::endl;
	      ofcurr << cand_vx[min_idx]->getVertexPoint() << std::endl;
	    }
#endif

	  if (min_idx >= 0 && cyl_rad > 0.0)
	    {
	      // Perform cylinder intersection
	      trim_segments = BoundedUtils::getCylinderIntersections(surf, centre2, 
								     axis2, cyl_rad,
								     epsge, bd_sf);
      
	      // Remove intersections not connected with the initial point
	      // Remove also segments going through an adjacent vertex
	      Point dummy;
	      checkTrimSeg(trim_segments, next_vxs, vx_point, dummy, tol2 /*epsge*/);
	    }

	  if (trim_segments.size() == 0 && min_idx >= 0)
	    {
	      // Check the feasability of a stright curve in the
	      // parameter domain
	      Point parval2 = cand_vx[min_idx]->getFacePar(face.get());
	      shared_ptr<ParamCurve> pcurve = checkStrightParCv(face, vx, cand_vx[min_idx], 
								epsge);
	      if (pcurve.get())
		{
		  trim_segments = BoundedUtils::getTrimCrvsPcrv(surf, pcurve, epsge,
								bd_sf);
		}
	      else
		{
		  // Find division curve between vertices
		  trim_segments = BoundedUtils::getTrimCrvsParam(surf, vx_par,
								 parval2, epsge,
								 bd_sf);
		}

	      // Check output
	      checkTrimSeg2(trim_segments, vx_par, parval2, epsge);
	      
	      // Remove intersections not connected with the initial point
	      // Remove also segments going through an adjacent vertex
	      if (trim_segments.size() > 0)
		{
		  Point other_pt = cand_vx[min_idx]->getVertexPoint();
		  checkTrimSeg(trim_segments, next_vxs, vx_point, 
			       other_pt, tol2 /*epsge*/);
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
      // if (trim_segments.size() == 0 && min_frac < fac*max_frac &&
      // 	   edge_par.dimension() == 2)
      //   {
      //     // Let the division curve end at an edge closest point
      //     trim_segments = BoundedUtils::getTrimCrvsParam(surf, vx_par,
      // 						     edge_par, epsge,
      // 						     bd_sf);
      //     // Check output
      //     Point par1, par2;
      //     for (kr=0; kr<trim_segments.size(); ++kr)
      // 	{
      // 	  par1 = 
      // 	    trim_segments[kr]->parameterCurve()->point(trim_segments[kr]->startparam());
      // 	   par2 = 
      // 	    trim_segments[kr]->parameterCurve()->point(trim_segments[kr]->endparam());
      // 	  if (par1.dist(edge_par) < epsge || par2.dist(edge_par) < epsge)
      // 	    break;
      // 	}
      //     if (kr == trim_segments.size() ||
      // 	  (trim_segments.size()>1 && 
      // 	   par1.dist(vx_par)>epsge && par2.dist(vx_par)>epsge))
      // 	trim_segments.clear();
	      
      //   }
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
      if (((ang1 < 0.25*ang2 && close_idx < 0) ||
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

	  // Adjust curves ending very close to a non-significant vertex
	  adjustTrimSeg(trim_segments, NULL, NULL, face, bd_sf, non_corner,
			tol2, epsge);

	  // Remove intersections not connected with the initial point
	  // Remove also segments going through an adjacent vertex
	  Point dummy;
	  checkTrimSeg(trim_segments, next_vxs, vx_point, dummy, tol2 /*epsge*/);

	  // Check configuration to avoid 3-sided surfaces
	  checkTrimConfig(face, trim_segments, vx, corners, tol2 /*epsge*/);
 	}
      if (trim_segments.size() == 0 && 
	  ((ang2 < 0.25*ang1  && close_idx < 0) ||
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

	  // Adjust curves ending very close to a non-significant vertex
	  adjustTrimSeg(trim_segments, NULL, NULL, face, bd_sf, non_corner,
			tol2, epsge);

	  // Remove intersections not connected with the initial point
	  // Remove also segments going through an adjacent vertex
	  Point dummy;
	  checkTrimSeg(trim_segments, next_vxs, vx_point, dummy, tol2 /*epsge*/);

	  // Check configuration to avoid 3-sided surfaces
	  checkTrimConfig(face, trim_segments, vx, corners, tol2 /*epsge*/);
 	}
    }

  if (opposite && trim_segments.size() == 0)
    {
      // Split to opposite edge. First fetch face parameter
      // Make sure to split in the inner of the edge
      // double ta = opposite->tMin();
      // double tb = opposite->tMax();
      // opposite_par = std::max(opposite_par, ta+0.2*(tb-ta));
      // opposite_par = std::min(opposite_par, tb-0.2*(tb-ta));
      Point face_par = opposite->faceParameter(opposite_par);
      opposite_point = opposite->point(opposite_par);
      trim_segments = BoundedUtils::getTrimCrvsParam(surf, vx_par,
						     face_par, epsge,
						     bd_sf);
 
      // Adjust curves ending very close to a non-significant vertex
      adjustTrimSeg(trim_segments, &vx_par, &face_par, face, bd_sf, non_corner,
		    tol2, epsge);

      // Remove intersections not connected with the initial point
      // Remove also segments going through an adjacent vertex
      checkTrimSeg(trim_segments, next_vxs, vx_point, 
		   opposite_point, tol2 /*epsge*/);
    }
      
  if (trim_segments.size() == 0)
    {
      // Connect to closest point
      trim_segments = BoundedUtils::getTrimCrvsParam(surf, vx_par,
						     close_par, epsge,
						     bd_sf);

      // Adjust curves ending very close to a non-significant vertex
      adjustTrimSeg(trim_segments, &vx_par, &close_par, face, bd_sf, non_corner,
		    tol2, epsge);

      // Remove intersections not connected with the initial point
      // Remove also segments going through an adjacent vertex
      checkTrimSeg(trim_segments, next_vxs, vx_point, close_pt, tol2 /*epsge*/);
    }

  if (trim_segments.size() == 0)
    {
      // Find intersections between the face and this plane
      trim_segments = BoundedUtils::getPlaneIntersections(surf, vx_point,
							  normal, epsge,
							  bd_sf);

      // Adjust curves ending very close to a non-significant vertex
      adjustTrimSeg(trim_segments, NULL, NULL, face, bd_sf, non_corner,
		    tol2, epsge);

       // Remove intersections not connected with the initial point
      // Remove also segments going through an adjacent vertex
      Point dummy;
      checkTrimSeg(trim_segments, next_vxs, vx_point, dummy, tol2 /*epsge*/);
    }

  if (trim_segments.size() == 0)
    {
      // No split. 
      return trim_segments;
    }

#ifdef DEBUG_REG
  std::ofstream out_file("split_segments.g2");
  for (size_t kj=0; kj<trim_segments.size(); ++kj)
    {
      shared_ptr<ParamCurve> cv = trim_segments[kj]->spaceCurve();
      cv->writeStandardHeader(out_file);
      cv->write(out_file);
    }
#endif

  return trim_segments;
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
	      shared_ptr<Vertex> v1 = edges[kh]->getVertex(true);
	      shared_ptr<Vertex> v2 = edges[kh]->getVertex(false);

	      if (vx_pt.dist(v1->getVertexPoint()) < epsge ||
		  vx_pt.dist(v2->getVertexPoint()) < epsge)
		continue;  // Non-corner vertex not transferred

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
				       Point& close_par, int loop_idx)
//==========================================================================
{
  close_idx = -1;
  close_dist = 1.0e8;

  // Fetch information about boundary curves
  size_t kr, kh;
  vector<shared_ptr<ftEdge> > all_edg = (loop_idx < 0) ?
    face->getAllEdges() :  face->getAllEdges(loop_idx);
  vector<shared_ptr<ParamCurve> > cvs;
  vector<Point> adj_pnt;  // Vertex points adjacent to the split vertex.
  // Can not serve as vertices to connect to and is thus not feasible as
  // a closest point to compare with
  for (kr=0; kr<all_edg.size(); ++kr)
    {
      shared_ptr<Vertex> next_vx = all_edg[kr]->getOtherVertex(vx.get());
      if (next_vx.get())
	adj_pnt.push_back(next_vx->getVertexPoint());
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

      // Check feasability of closest point
      for (kh=0; kh<adj_pnt.size(); ++kh)
	if (adj_pnt[kh].dist(close_pnt) < epsge)
	  break;
      if (kh < adj_pnt.size())
	continue;

      if (dist < close_dist)
	{
	  close_dist = dist;
	  close_idx = (int)kr;

	  surf->closestPoint(close_pnt, upar, vpar, close_pnt2, dist,
			     epsge);
	  close_par = Point(upar,vpar);
	}
    }

  if (close_idx < 0)
    {
      // No closest point selected. Choose a point
      close_idx = (int)cvs.size()/2;
      close_t = 0.5*(cvs[close_idx]->startparam() + cvs[close_idx]->endparam());
      close_pnt = cvs[close_idx]->point(close_t);
      double upar, vpar, dist;
      Point close_pnt2;
      surf->closestPoint(close_pnt, upar, vpar, close_pnt2, dist,
			 epsge);
      close_par = Point(upar,vpar);
    }
}


//==========================================================================
bool
  RegularizeUtils::checkPath(shared_ptr<Vertex> vx1, shared_ptr<Vertex> vx2,
			     shared_ptr<Vertex> vx, shared_ptr<ftSurface> face,
			     double angtol)
//==========================================================================
{
  // Find the path between vx1 and vx
  ftEdge *edg = vx1->getCommonEdge(vx2.get());
  vector<ftEdge*> path;
  bool found = getPath(edg, vx2, vx, face, path);

  if (!found)
    return true; //false;   // Should be a path, if not it is probably not a good split

  // Fetch corners
  vector<shared_ptr<Vertex> > vx_corners;
  size_t kj;
  vx_corners.push_back(vx);
  for (kj=path.size()-1; kj>0; --kj)
    {
      shared_ptr<Vertex> common_vx = 
	path[kj-1]->getCommonVertex(path[kj]);
      double t1 = path[kj-1]->parAtVertex(common_vx.get());
      double t2 = path[kj]->parAtVertex(common_vx.get());
      Point tan1 = path[kj-1]->tangent(t1);
      Point tan2 = path[kj]->tangent(t2);
      double ang = tan1.angle(tan2);
      if (ang > angtol)
	vx_corners.push_back(common_vx);
    }
  vx_corners.push_back(vx1);

  if (vx_corners.size() > 4)
    return true;   // Does not lead to a regular face, not necessary to 
  // test quality
  
  /*bool OK =*/ checkRegularity(vx_corners, face, false);
  //return OK;
  return true;
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
  int min_nmb_corners = 200;  // A large and unlikely number
  size_t ki, kj;
  for (ki=0; ki<edg1.size(); ++ki)
    {
      shared_ptr<Vertex> other_vx = edg1[ki]->getOtherVertex(vx1.get());
      vector<ftEdge*> path;
      bool found = getPath(edg1[ki], other_vx, vx2, face, path);
      if (found)
	{
	  // Count corners
	  int nmb_corners = 0;
	  for (kj=1; kj<path.size(); ++kj)
	    {
	      shared_ptr<Vertex> common_vx = 
		path[kj-1]->getCommonVertex(path[kj]);
	      double t1 = path[kj-1]->parAtVertex(common_vx.get());
	      double t2 = path[kj]->parAtVertex(common_vx.get());
	      Point tan1 = path[kj-1]->tangent(t1);
	      Point tan2 = path[kj]->tangent(t2);
	      double ang = tan1.angle(tan2);
	      if (ang > angtol)
		nmb_corners++;
	    }
	  if (ki == 0 || nmb_corners < min_nmb_corners)
	    min_nmb_corners = nmb_corners;
	}
    }

  if (min_nmb_corners >= 2)
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
int
  RegularizeUtils::noExtension(shared_ptr<Vertex> vx, ftSurface* face,
			       shared_ptr<Vertex>& vx2, pair<Point, Point>& co_par1, 
			       pair<Point, Point>& co_par2, int& dir1, int& dir2,
			       double& val1, double& val2, double angtol, 
			       bool check_constant_curve)
//==========================================================================
{
  double ptol = 1.0e-8;
  size_t ki, kj, kr;
  int idx = -1;
  shared_ptr<Vertex> last_vx;

 // Fetch faces and remove faces from other bodies
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
  if (vx_faces.size() == 2)
    return 1;  // Not a significant vertex in this body

  if (vx_faces.size() != 3)
    return 0;

  // Get edges
  vector<ftEdge*> edges = vx->uniqueEdges();

  // Look for an edge connecting two T-vertices where 3 edges meet
  for (ki=0; ki<edges.size(); ++ki)
    {
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

      // Dismiss edges going along the initial face
      if (edges[ki]->face() == face ||
	  edges[ki]->twin()->geomEdge()->face() == face)
	continue;

      // Check configuration
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

      if (vx_faces2.size() < 3)
	{
	  // Traverse along the edge until a T-joint is found
#ifdef DEBUG_REG
	  std::ofstream of("traverse_vx.g2");
	  of << "400 1 0 4 255 0 0 255" << std::endl;
	  of << "1" << std::endl;
	  of << vx->getVertexPoint() << std::endl;
	  of << "400 1 0 4 0 255 0 255" << std::endl;
	  of << "1" << std::endl;
	  of << vx2->getVertexPoint() << std::endl;
#endif

	  int status0 = traverseUntilTJoint(vx_faces, vx, vx2, vx_faces2);
	  if (status0 > 1)
	    continue;  // Not a smooth transition
	}

      bool continued_merge = false;
      if (vx_faces2.size() == 4)
	{
	  // Check the continuation of this edge to see if there are
	  // several merge situations in a row
	  continued_merge = mergeSituationContinuation(face, vx, edges[ki], angtol);
	}

      if (vx_faces2.size() != 3 && !continued_merge)
	continue;

      // Remove the face(s) which are adjacent to the initial vertex
      for (kj=0; kj<vx_faces2.size();)
	{
	  for (kr=0; kr<vx_faces.size(); ++kr)
	    if (vx_faces[kr] == vx_faces2[kj])
	      break;
	  if (kr < vx_faces.size())
	    {
	      // vx_faces2.erase(vx_faces2.begin()+kj); Something wrong with this!!
	      std::swap(vx_faces2[kj], vx_faces2[vx_faces2.size()-1]);
	      vx_faces2.pop_back();
	    }
	  else
	    kj++;
	}
      if ( vx_faces2.size() == 0 || vx_faces2.size() > 2)
	continue;  // Not a legal configuration

      // Fetch edges in this vertex
      vector<ftEdge*> edges2 = vx2->uniqueEdges();
      
      // Remove the edges that do not follow the remaining faces only once
      for (kj=0; kj<edges2.size(); )
	{
	  int nmb = 0;
	  for (kr=0; kr<vx_faces2.size(); ++kr)
	    {
	      if (edges2[kj]->face() == vx_faces2[kr] ||
		  (edges2[kj]->twin() && 
		   edges2[kj]->twin()->geomEdge()->face() == vx_faces2[kr]))
		nmb++;
	    }
	  if (nmb == 1)
	    kj++;
	  else
	    edges2.erase(edges2.begin()+kj);
	}
      if (edges2.size() != 2)
	continue;

      // Check angle
      double t1 = edges2[0]->parAtVertex(vx2.get());
      double t2 = edges2[1]->parAtVertex(vx2.get());
      Point tan1 = edges2[0]->tangent(t1);
      Point tan2 = edges2[1]->tangent(t2);
      double ang = tan1.angle(tan2);
      ang = std::min(ang, fabs(M_PI-ang));
      if (ang > angtol)
	continue;  // The faces meet in a corner

      idx = ki;
      last_vx = vx2;
    }
  if (idx < 0)
    return 0;

  vx2 = last_vx;
  ftSurface *f1 = edges[idx]->face()->asFtSurface();
  ftSurface *f2 = edges[idx]->twin()->geomEdge()->face()->asFtSurface();
  Point par1_1 = vx->getFacePar(f1);
  Point par1_2 = vx->getFacePar(f2);
  Point par2_1 = vx2->getFacePar(f1);
  Point par2_2 = vx2->getFacePar(f2);
  co_par1 = make_pair(par1_1, par1_2);
  co_par2 = make_pair(par2_1, par2_2);

  // Update output
  vx2 = edges[idx]->getOtherVertex(vx.get());

  shared_ptr<ParamCurve> cv1 = edges[idx]->geomCurve();
  shared_ptr<ParamCurve> cv2 = edges[idx]->twin()->geomEdge()->geomCurve();
  shared_ptr<CurveOnSurface> sf_cv1 = 
    dynamic_pointer_cast<CurveOnSurface,ParamCurve>(cv1);
  shared_ptr<CurveOnSurface> sf_cv2 = 
    dynamic_pointer_cast<CurveOnSurface,ParamCurve>(cv2);

  // Point pt1 = cv1->point(edges[idx]->tMin());
  // Point pt2 = cv2->point(edges[idx]->twin()->tMin());
  // Point pt3 = cv2->point(edges[idx]->twin()->tMax());
  // if (pt1.dist(pt2) < pt1.dist(pt3))
  //   {
  //     co_par1 = make_pair(sf_cv1->parameterCurve()->point(edges[idx]->tMin()),
  // 			  sf_cv2->parameterCurve()->point(edges[idx]->twin()->tMin()));
  //     co_par2 = make_pair(sf_cv1->parameterCurve()->point(edges[idx]->tMax()),
  // 			  sf_cv2->parameterCurve()->point(edges[idx]->twin()->tMax()));
  //   }
  // else
  //   {
  //     co_par1 = make_pair(sf_cv1->parameterCurve()->point(edges[idx]->tMin()),
  // 			  sf_cv2->parameterCurve()->point(edges[idx]->twin()->tMax()));
  //     co_par2 = make_pair(sf_cv1->parameterCurve()->point(edges[idx]->tMax()),
  // 			  sf_cv2->parameterCurve()->point(edges[idx]->twin()->tMin()));
  //   }

  if (check_constant_curve)
    {

      // Check if the boundary between the
      // faces are constant parameter curves in both parameter directions
      if (!sf_cv1.get() || 
	  !sf_cv1->isConstantCurve(ptol, dir1, val1))
	return 1;
      if (!sf_cv2.get() || 
	  !sf_cv2->isConstantCurve(ptol, dir2, val2))
	return 1;

      dir1--;
      dir2--;
	  
      return 2;   // A possible insignificant splitting curve is found
    }

      // if (false)
      // 	{
      // 	  // Check if the next or previous edge follow the same constant boundary
      // 	  int dir;
      // 	  double val;
      // 	  shared_ptr<ParamCurve> tmp_crv = edges[idx]->next()->geomEdge()->geomCurve();
      // 	  shared_ptr<CurveOnSurface> sf_crv = 
      // 	    dynamic_pointer_cast<CurveOnSurface,ParamCurve>(tmp_crv);
      // 	  if (sf_crv.get())
      // 	    sf_crv->isConstantCurve(ptol, dir, val);
      // 	  else
      // 	    dir = -1;
      // 	  if (dir == dir1)
      // 	    continue;
      // 	  tmp_crv = edges[idx]->prev()->geomEdge()->geomCurve();
      // 	  sf_crv = dynamic_pointer_cast<CurveOnSurface,ParamCurve>(tmp_crv);
      // 	  if (sf_crv.get())
      // 	    sf_crv->isConstantCurve(ptol, dir, val);
      // 	  else
      // 	    dir = -1;
      // 	  if (dir == dir1)
      // 	    continue;
      // 	  tmp_crv = edges[idx]->twin()->next()->geomEdge()->geomCurve();
      // 	  sf_crv = dynamic_pointer_cast<CurveOnSurface,ParamCurve>(tmp_crv);
      // 	  if (sf_crv.get())
      // 	    sf_crv->isConstantCurve(ptol, dir, val);
      // 	  else
      // 	    dir = -1;
      // 	  if (dir == dir2)
      // 	    continue;
      // 	  tmp_crv = edges[idx]->twin()->prev()->geomEdge()->geomCurve();
      // 	  sf_crv = dynamic_pointer_cast<CurveOnSurface,ParamCurve>(tmp_crv);
      // 	  if (sf_crv.get())
      // 	    sf_crv->isConstantCurve(ptol, dir, val);
      // 	  else
      // 	    dir = -1;
      // 	  if (dir == dir2)
      // 	    continue;
      // 	}

  return 1;
}


//==========================================================================
bool
RegularizeUtils::mergeSituationContinuation(ftSurface* init_face, shared_ptr<Vertex> vx,
					    ftEdge* edge, double angtol)
//==========================================================================
{
  ftSurface *face1, *face2;
  
  // Fetch faces adjacent to edge
  face1 = edge->face()->asFtSurface();
  if (!edge->twin())
    return false;
  face2 = edge->twin()->geomEdge()->face()->asFtSurface();

  shared_ptr<Vertex> vx1 = vx;
  ftEdge* edge1 = edge;
  size_t ki;
  while (true)
    {
      // Check continuity
      if (!edge1->hasConnectivityInfo())
	return false;  // No continuity info
      shared_ptr<FaceConnectivity<ftEdgeBase> > info = 
	edge1->getConnectivityInfo();
      int status = info->WorstStatus();
      if (status > 1)
	return false;  // Not G1 or almost G1
      
      // Fetch next vertex
      shared_ptr<Vertex> vx2 = edge1->getOtherVertex(vx1.get());

      // Fetch faces surrounding the vertex
      vector<ftSurface*> vx_faces2 = vx2->faces();
	        
      // Get edges
      vector<ftEdge*> edges = vx2->uniqueEdges();

      if (vx_faces2.size() == 3)
	{
	  // Check continuity at end of edge sequence
	  // Remove the edges that do not follow the recent faces only once
	  for (ki=0; ki<edges.size(); )
	    {
	      int nmb = 0;
	      if (edges[ki]->face() == face1 ||
		  (edges[ki]->twin() && 
		   edges[ki]->twin()->geomEdge()->face() == face1))
		nmb++;
	      if (edges[ki]->face() == face2 ||
		  (edges[ki]->twin() && 
		   edges[ki]->twin()->geomEdge()->face() == face2))
		nmb++;
	      if (nmb == 1)
		ki++;
	      else
		edges.erase(edges.begin()+ki);
	    }
	  if (edges.size() != 2)
	    return false;

	  // Check angle
	  double t1 = edges[0]->parAtVertex(vx2.get());
	  double t2 = edges[1]->parAtVertex(vx2.get());
	  Point tan1 = edges[0]->tangent(t1);
	  Point tan2 = edges[1]->tangent(t2);
	  double ang = tan1.angle(tan2);
	  ang = std::min(ang, fabs(M_PI-ang));
	  if (ang > angtol)
	    return false;  // The faces meet in a corner
	  
	  // Check if the initial face is found
	  int kj;
	  for (kj=0; kj<vx_faces2.size(); ++kj)
	    if (vx_faces2[kj] == init_face)
	      return false;
	  return true;  // The end of the edge sequence is found
	}

      // Remove edges that cannot be a part of a continuation
      if (vx_faces2.size() == 2)
	{
	  for (ki=0; ki<edges.size(); )
	    {
	      if (edges[ki] == edge1 || edges[ki]->twin() == edge1)
		edges.erase(edges.begin()+ki);
	      else
		ki++;
	    }
	}
      else
	{
	  for (ki=0; ki<edges.size(); )
	    {
	      if (edges[ki]->face() == face1 || edges[ki]->face() == face2 ||
		  (!edges[ki]->twin()))
		edges.erase(edges.begin()+ki);
	      else
		{
		  ftFaceBase *tmp_face = edges[ki]->twin()->geomEdge()->face();
		  if (tmp_face == face1 || tmp_face == face2)
		    edges.erase(edges.begin()+ki);
		  else
		    ki++;
		}
	    }
	}
      
      if (edges.size() != 1)
	return false;
      edge1 = edges[0];
      if (edge1 == edge)
	return false;   // A loop is found

      // Fetch faces adjacent to edge
      face1 = edge1->face()->asFtSurface();
      if (!edge1->twin())
	return false;
      face2 = edge1->twin()->geomEdge()->face()->asFtSurface();

      vx1 = vx2;
    }
  return true;  // Should not get here
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

 
//==========================================================================
int
RegularizeUtils::selectCandVx(shared_ptr<ftSurface> face,
			      shared_ptr<Vertex> vx, const Point& in_vec,
			      vector<shared_ptr<Vertex> > cand_vx,
			      RectDomain& dom,
			      double epsge, double angtol, 
			      const Point& centre, const Point& normal,
			      vector<shared_ptr<ParamCurve> >& vx_cvs,
			      double close_dist, const Point& close_pt,
			      double& cyl_rad, bool strong)
//==========================================================================
{
  // Statistics
  double parfrac = getMaxParFrac(face);
  double level_frac = std::max(std::min(0.5, 100.0*parfrac), 0.1);
  Point edge_par;
  //double level_ang = M_PI/3; // M_PI/2.0; // M_PI/4.0; //M_PI/6.0;
  double level_ang = (centre.dimension() >  0) ? M_PI/3 : M_PI/2.0; 
  double fac = 0.5; // 0.2;
  double fac2 = 2.0;
  double fac3 = 0.05;
  double tol = 1.0e-4;
   
  Point vx_point = vx->getVertexPoint();
  Point vx_par = vx->getFacePar(face.get());
  Point close_vec = close_pt - vx_point;
  Point curr_vx_par, curr_vx_par1;
  Point min_deriv;

  int min_idx = -1;
  double min_frac = MAXDOUBLE;
  double max_frac = 0.0;
  double min_ang = 1.0e8;
  double min_close_ang = 1.0e8;
  double min_dist = MAXDOUBLE;
  double curr_rad_dist = MAXDOUBLE;
  double min_close_dist = MAXDOUBLE;
  double d1=-1.0, d2=-1.0;
  if (centre.dimension() > 0)
    d1 = vx_point.dist(centre);

  // Fetch adjacent vertices
  vector<shared_ptr<Vertex> > next_vxs = vx->getNextVertex(face.get());

  // Tangent info
  shared_ptr<ParamSurface> surf = face->surface();
  DirectionCone cone1 = surf->tangentCone(true);
  DirectionCone cone2 = surf->tangentCone(false);

  // Traverse all candidate vertices and check if any is feasible for 
  // division
  // Compute distance between vertices and found plane
  // Select vertex with minimum distance
   for (int ki=0; ki<(int)cand_vx.size(); ++ki)
    {
      size_t kr, kh;
      Point curr_vx_par2 = cand_vx[ki]->getFacePar(face.get());
      vector<Point> der(3);
      surf->point(der, curr_vx_par2[0], curr_vx_par2[1], 1);
      Point cand_vx_pt = cand_vx[ki]->getVertexPoint();
      Point vec = cand_vx_pt - vx_point;
      double dist = vec.length();
      double dist1 = fabs(vx_par[0]-curr_vx_par2[0]);
      double dist2 = fabs(vx_par[1]-curr_vx_par2[1]);
      dist1 /= (dom.umax() - dom.umin());
      dist2 /= (dom.vmax() - dom.vmin());
      double frac = std::min(dist1, dist2);
      double ang = vec.angle(normal);
      ang = fabs(0.5*M_PI - ang);
      double close_dist = close_pt.dist(cand_vx_pt);
      double rad_dist;
      min_close_ang = std::min(min_close_ang, close_vec.angle(vec));
      if (centre.dimension() > 0)
	{
	  d2 = cand_vx[ki]->getVertexPoint().dist(centre);
	  rad_dist = fabs(d1 - d2);
	}
      else 
	rad_dist = MAXDOUBLE;

      // Skip vertices lying in the wrong direction compared to the material
      // of the surface
      if (vec*in_vec < -tol)
	continue;   // The tolerances is arbitrary here, but do not want
      // to rule out orthogonal cases when the face is curved. The
      // test should be made more precise by taking the shape of the
      // face into consideration

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
      if ((!cornerInShortestPath(vx, cand_vx[ki], face, angtol)) && (!strong))
	continue;

      // Avoid vertices which is inline with the neighbouring vertices
      // to the split vertex
      Point vec3 = cand_vx[ki]->getVertexPoint() - vx_point;
     for (kr=0; kr<next_vxs.size(); ++kr)
	{
	  Point vec4 = next_vxs[kr]->getVertexPoint() - vx_point;
	  double vx_ang = vec3.angle(vec4);
	  if (vx_ang < angtol)
	    break;
	}
     if (kr < next_vxs.size() && cyl_rad <= 0.0)
       {
	 // This candidate vertex is not a good coice for splitting
	 continue;
       }
     
     // Check also other angles
     vector<shared_ptr<Vertex> > next_vxs2 = 
       cand_vx[ki]->getNextVertex(face.get());
     for (kr=0; kr<next_vxs2.size(); ++kr)
       {
	 bool OK = checkPath(cand_vx[ki], next_vxs2[kr], vx,
			     face, angtol);
	 if (!OK)
	   break;
       }
     if (kr < next_vxs2.size() && cyl_rad <= 0.0)
       {
	 // This candidate vertex is not a good coice for splitting
	 continue;
       }

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
	      max_frac = std::max(dist1,dist2);
	      min_deriv = (dist1 > dist2) ? der[2] : der[1];
	      min_close_dist = close_dist;
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
	      max_frac = std::max(dist1,dist2);
	      min_deriv = (dist1 > dist2) ? der[2] : der[1];
	      min_close_dist = close_dist;
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
	      max_frac = std::max(dist1,dist2);
	      min_deriv = (dist1 > dist2) ? der[2] : der[1];
	      min_close_dist = close_dist;
	    }
	  else if ((par_ang < par_limit || close_dist < min_close_dist) && 
		   dist < min_dist)
	    {
	      curr_rad_dist = rad_dist;
	      min_ang = ang;
	      min_dist = dist;
	      min_frac = frac;
	      min_idx = ki;
	      curr_vx_par = curr_vx_par2;
	      max_frac = std::max(dist1,dist2);
	      min_deriv = (dist1 > dist2) ? der[2] : der[1];
	      min_close_dist = close_dist;
	    }
	}

      
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
	  curr_vx_par = curr_vx_par1 = curr_vx_par2;
	}
    }
   
   double ang = 0.0;
   double ang2 = 0.0;
   double ang3 = 0.0;
   Point vec;
   if (min_idx >= 0)
     {
       vec = cand_vx[min_idx]->getVertexPoint() - vx_point;
       ang = vec.angle(normal);

       // Compute also the angle in the candidate end point of the split
       Point pnt2, normal2;
       getDivisionPlane(face, cand_vx[min_idx], epsge, pnt2, normal2);
       ang2 = vec.angle(normal2);
       ang3 = normal.angle(normal2);
     }

   if (min_idx >= 0)
     {
       double deriv_ang1 = min_deriv.angle(normal);
       deriv_ang1 = fabs(0.5*M_PI - deriv_ang1);
       double deriv_ang2 = min_deriv.angle(close_vec);
       double close_ang = close_vec.angle(vec);
       // if (min_deriv.angle(vec) > level_frac*level_ang && 
       // 	   close_ang > 2.0*min_close_ang)
       // 	 min_frac = max_frac;   // Not the same parameter direction
       if (curr_rad_dist < epsge)
	 {
	   // A good candidate for a cylinder split
	   cyl_rad = d1;
	 }
       else
	 {
	   double scp = 
	     (vx_point - close_pt)*(vx_point - cand_vx[min_idx]->getVertexPoint());
	   double fac4 = (scp < 0.0) ? 0.5 : 1;
	   if (false
	       /*(fabs(ang - ang2) < epsge || fabs(M_PI-(ang+ang2)) < epsge) &&
		 (ang3 < epsge || fabs(M_PI-ang3) < epsge)*/)
	     {
	       if (!(vx->isCornerInFace(face.get(), angtol) && 
		     cand_vx[min_idx]->isCornerInFace(face.get(), angtol)))
		 level_frac *= 0.1;
	     }
	   if (strong)
	     {
	       // The candidate is selected already
	       ;
	     }
	   else if ((!((fac4*fac*min_dist < close_dist ||
			min_frac < fac3*max_frac) &&
		       (min_frac < level_frac*max_frac || 
			fabs(0.5*M_PI - ang) < level_frac*level_ang))) /*||
		    (deriv_ang1 > level_frac*level_ang && 
		     deriv_ang2 > level_frac*level_ang &&
		     close_ang > 2.0*min_close_ang)*/)
	     {
	       // Not a good candidate. Remove it from the list
	       cand_vx.erase(cand_vx.begin()+min_idx);
	       min_idx = -1;
	     }
	   else
	     {
	       int stop_break = 1;
	       stop_break *= 2;
	     }
	 }
     }

   return min_idx;
}

//==========================================================================
void 
RegularizeUtils::adjustTrimSeg(vector<shared_ptr<CurveOnSurface> >& trim_segments,
			       Point *parval1, Point *parval2,
			       shared_ptr<ftSurface> face,
			       shared_ptr<BoundedSurface>& bd_sf,
			       vector<shared_ptr<Vertex> >& non_corner,
			       double tol, double epsge)
      // Avoid unstability be dividing very close to an insignificant vertex
//==========================================================================
{
  shared_ptr<ParamSurface> surf = face->surface();
   for (size_t kr=0; kr<trim_segments.size(); ++kr)
    {
      Point pos1 = trim_segments[kr]->ParamCurve::point(trim_segments[kr]->startparam());
      Point pos2 = trim_segments[kr]->ParamCurve::point(trim_segments[kr]->endparam());

      double min_dist1=HUGE, min_dist2=HUGE;
      int idx1 = -1, idx2 = -1;
      for (size_t kj=0; kj<non_corner.size(); ++kj)
	{
	  Point pos3 = non_corner[kj]->getVertexPoint();
	  double dist = pos1.dist(pos3);
	  if (dist < min_dist1)
	    {
	      min_dist1 = dist;
	      idx1 = (int)kj;
	    }
	  dist = pos2.dist(pos3);
	  if (dist < min_dist2)
	    {
	      min_dist2 = dist;
	      idx2 = (int)kj;
	    }
	}

      bool replace = false;
      Point par1, par2;
      Point dummy_vec;
      if (idx1 >= 0 && min_dist1 < tol)
	{
	  // Modify trim segment by adjusting the end points and make a
	  // constant parameter curve
	  replace = true;
	  par1 = non_corner[idx1]->getFacePar(face.get());
	  if (parval2 == NULL)
	    {
	      double upar, vpar, dist, edg_par;
	      Point clo_pt;
	      ftEdgeBase *edg = face->closestBoundaryPoint(pos2, dummy_vec, upar, 
							   vpar, clo_pt, 
							   dist, edg_par);
	      par2 = Point(upar, vpar);
	    }
	  else
	    par2 = (*parval2);
	}

      if (idx2 >= 0 && min_dist2 < tol)
	{
	  replace = true;
	  if (parval1 == NULL)
	    {
	      double upar, vpar, dist, edg_par;
	      Point clo_pt;
	      ftEdgeBase *edg = face->closestBoundaryPoint(pos1, dummy_vec, upar, 
							   vpar, clo_pt, 
							   dist, edg_par);
	      par1 = Point(upar, vpar);
	    }
	  else
	    par1 = (*parval1);
	  par2 = non_corner[idx2]->getFacePar(face.get());
	}

      if (replace)
	{
	  // Modified curve
	  vector<shared_ptr<CurveOnSurface> > mod_seg =
	    BoundedUtils::getTrimCrvsParam(surf, par1, par2, epsge, bd_sf);
	  if (mod_seg.size() == 1)
	    {
	      trim_segments[kr] = mod_seg[0];
	      break;
	    }
	}
    }
}


 //==========================================================================
void 
RegularizeUtils::checkTrimSeg(vector<shared_ptr<CurveOnSurface> >& trim_segments,
			      vector<shared_ptr<Vertex> >& next_vxs,
			      const Point& vx_point, const Point& other_pt,
			      double epsge)
      // Remove intersections not connected with the initial point
      // Remove also segments going through an adjacent vertex
//==========================================================================
{
  for (size_t kr=0; kr<trim_segments.size(); )
    {
      Point pos1 = trim_segments[kr]->ParamCurve::point(trim_segments[kr]->startparam());
      Point pos2 = trim_segments[kr]->ParamCurve::point(trim_segments[kr]->endparam());
      if (pos1.dist(vx_point) > epsge && pos2.dist(vx_point) > epsge)
	trim_segments.erase(trim_segments.begin()+kr);
      else if (other_pt.dimension() == vx_point.dimension() &&
	       pos1.dist(other_pt) > epsge && pos2.dist(other_pt) > epsge)
	trim_segments.erase(trim_segments.begin()+kr);
      else
	{
	  double ta = trim_segments[kr]->startparam();
	  double tb = trim_segments[kr]->endparam();
	  size_t k3;
	  for (k3=0; k3<next_vxs.size(); ++k3)
	    {
	      Point next_pt = next_vxs[k3]->getVertexPoint();
	      Point next_close;
	      double next_par, next_dist;
	      trim_segments[kr]->closestPoint(next_pt, ta, tb, next_par, 
					      next_close, next_dist);
	      if (next_dist < epsge && next_par > ta && next_par < tb)
		break;
	    }
	  if (k3 < next_vxs.size())
	    trim_segments.erase(trim_segments.begin()+kr);
	  else
	    kr++;
	}
    }

}

//==========================================================================
void 
RegularizeUtils::checkTrimSeg2(vector<shared_ptr<CurveOnSurface> >& trim_segments,
			       const Point& vx_par1, const Point& vx_par2, 
			       double epsge)
// Remove intersections not connected with the initial points in the 
// parameter domain
//==========================================================================
{
  Point par1, par2;
  size_t kr;
  for (kr=0; kr<trim_segments.size(); ++kr)
    {
      par1 = 
	trim_segments[kr]->parameterCurve()->point(trim_segments[kr]->startparam());
      par2 = 
	trim_segments[kr]->parameterCurve()->point(trim_segments[kr]->endparam());
      if (par1.dist(vx_par2) < epsge || par2.dist(vx_par2) < epsge)
	break;
    }
  if (kr == trim_segments.size() ||
      (trim_segments.size()>1 && 
       par1.dist(vx_par1)>epsge && par2.dist(vx_par1)>epsge))
    trim_segments.clear();
}

//==========================================================================
void 
RegularizeUtils::checkTrimSeg3(vector<shared_ptr<CurveOnSurface> >& trim_segments,
			       const Point& vx_par1, const Point& vx_par2, 
			       double epsge)
// Remove intersections not connected with the initial points in the 
// parameter domain
//==========================================================================
{
  Point par1, par2;
  size_t kr;
  for (kr=0; kr<trim_segments.size();)
    {
      par1 = 
	trim_segments[kr]->parameterCurve()->point(trim_segments[kr]->startparam());
      par2 = 
	trim_segments[kr]->parameterCurve()->point(trim_segments[kr]->endparam());
      if (!((par1.dist(vx_par2) < epsge || par2.dist(vx_par2) < epsge) &&
	    (par1.dist(vx_par1) < epsge || par2.dist(vx_par1) < epsge)))
	trim_segments.erase(trim_segments.begin()+kr);
      else
	++kr;
    }
}

//==========================================================================
void 
RegularizeUtils::checkTrimConfig(shared_ptr<ftSurface> face,
				 vector<shared_ptr<CurveOnSurface> >& trim_segments,
				 shared_ptr<Vertex> vx,
				 vector<shared_ptr<Vertex> >& corners,
				 double epsge)
      // Remove intersections that would lead to 3-sided surface
//==========================================================================
{
  // Fetch the edge corresponding to end points of trim segments and count the
  // number of corners between the given vertex and the end point of the segment
  Point vx_point = vx->getVertexPoint();
  vector<ftEdge*> vx_edges = vx->getFaceEdges(face.get());
  if (vx_edges.size() != 2)
    return;  // An unexpected number of edges meeting in vertex

  Point dummy_vec;
   for (size_t kr=0; kr<trim_segments.size(); )
    {
      Point pos1 = trim_segments[kr]->ParamCurve::point(trim_segments[kr]->startparam());
      Point pos2 = trim_segments[kr]->ParamCurve::point(trim_segments[kr]->endparam());
      
      if (pos1.dist(vx_point) > epsge && pos2.dist(vx_point) > epsge)
	{
	  kr++;  // Not a candidate for this test
	  continue;
	}

      Point pos = (pos1.dist(vx_point) < pos2.dist(vx_point)) ? pos2 : pos1;
      Point close;
      double upar, vpar, par, dist;
      ftEdgeBase* tmp_edge = face->closestBoundaryPoint(pos, dummy_vec, upar, vpar, close,
							dist, par);
      ftEdge* edge2 = tmp_edge->geomEdge();
      if (!edge2)
	{
	  kr++;  // Something mysterious
	  continue;
	}

      // Check if the trim_segment joins in an edge
      double tol = std::max(1.0e-10, 1.0e-6*(pos1.dist(pos2)));
      shared_ptr<Vertex> other_vx;
      shared_ptr<Vertex> vx2_1, vx2_2;
      edge2->getVertices(vx2_1, vx2_2);
      if (pos.dist(vx2_1->getVertexPoint()) < tol)
	other_vx = vx2_1;
      else if (pos.dist(vx2_2->getVertexPoint()) < tol)
	other_vx = vx2_2;

      // Count number of corners in both directions from the given vertex
      ftEdge* edge1 = vx_edges[0];
      bool forward = (vx.get() == edge1->getVertex(true).get());
      int nmbc1 = 0;
      shared_ptr<Vertex> v1 = edge1->getOtherVertex(vx.get());
      while (edge1 != edge2)
	{
	  size_t kj;
	  for (kj=0; kj<corners.size(); ++kj)
	    if (v1.get() == corners[kj].get())
	      {
		nmbc1++;
		break;
	      }

	  edge1 = (forward) ? edge1->next()->geomEdge() : edge1->prev()->geomEdge();
	  v1 = edge1->getOtherVertex(v1.get());
	  if (edge1 == vx_edges[0])
	    break;
	  if (v1.get() == other_vx.get())
	    break;
	}

      edge1 = vx_edges[1];
      forward = (vx.get() == edge1->getVertex(true).get());
      int nmbc2 = 0;
      v1 = edge1->getOtherVertex(vx.get());
      while (edge1 != edge2)
	{
	  size_t kj;
	  for (kj=0; kj<corners.size(); ++kj)
	    if (v1.get() == corners[kj].get())
	      {
		nmbc2++;
		break;
	      }

	  edge1 = (forward) ? edge1->next()->geomEdge() : edge1->prev()->geomEdge();
	  v1 = edge1->getOtherVertex(v1.get());
	  if (edge1 == vx_edges[1])
	    break;
	  if (v1.get() == other_vx.get())
	    break;
	}

       if (nmbc1 < 2 || nmbc2 < 2)
	trim_segments.erase(trim_segments.begin()+kr);
      else
	kr++;
    }
  
}

//==========================================================================
ftEdge* RegularizeUtils::getOppositeBoundaryPar(shared_ptr<ftSurface> face,
					     shared_ptr<Vertex> vx, 
					     vector<shared_ptr<Vertex> >& corners,
					     double epsge, Point& point, 
					     double& par, double& dist)
//==========================================================================
{
  int cx = -1;
  for (cx=0; cx<(int)corners.size(); ++cx)
    if (corners[cx].get() == vx.get())
      break;

  int id1, id2;
  if (cx == (int)corners.size() && corners.size() != 4)
    return NULL;  // vx is not a corner
  else if (cx == (int)corners.size())
    {
      // Fetch adjacent corners
      vector<ftEdge*> edges = vx->getFaceEdges(face.get());
      if (edges.size() != 2)
	return NULL;  // An unexpected number of edges meeting in vertex
      ftEdge* edge1 = edges[0];
      ftEdge* edge2 = edges[1];
      shared_ptr<Vertex> v1 = edge1->getOtherVertex(vx.get());
      shared_ptr<Vertex> v2 = edge2->getOtherVertex(vx.get());

      // Identify corners corresponding to the adjacent edges
      id1 = id2 = -1;
      while (v1.get() != vx.get())
	{
	  cx = -1;
	  for (cx=0; cx<(int)corners.size(); ++cx)
	    if (corners[cx].get() == v1.get())
	      break;
	  if (cx  == (int)corners.size())
	    {
	      if (v1.get() == edge1->getVertex(true).get())
		edge1 = edge1->prev()->geomEdge();
	      else
		edge1 = edge1->next()->geomEdge();
	      v1 = edge1->getOtherVertex(v1.get());
	    }
	  else
	    {
	      id1 = cx;
	      break;
	    }
	}

      while (v2.get() != vx.get())
	{
	  cx = -1;
	  for (cx=0; cx<(int)corners.size(); ++cx)
	    if (corners[cx].get() == v2.get())
	      break;
	  if (cx  == (int)corners.size())
	    {
	      if (v2.get() == edge2->getVertex(true).get())
		edge2 = edge2->prev()->geomEdge();
	      else
		edge2 = edge2->next()->geomEdge();
	      v2 = edge2->getOtherVertex(v2.get());
	    }
	  else
	    {
	      id2 = cx;
	      break;
	    }
	}

      if (id1 == -1 || id2 == -1)
	return NULL;  // No corners found

      // Check consistency
      if (id2 > id1)
	std::swap(id1, id2);
      if (id2 == 0 && id1 == 3)
	std::swap(id1, id2);
      if (!(id1-id2 == 1 || (id1 == 0 && id2 == 3)))
	return NULL;
	  
      // Set indices
      id1 = (int)((id1+1)%corners.size());
      id2--;
     if (id2 < 0)
	id2 = (int)corners.size()+id2;
      }
  else
    {
      // Set indices of corners on each side of the relevant edges. 
      id1 = (int)((cx + 2)%corners.size());
      id2 = (cx - 2);
      if (id2 < 0)
	id2 = (int)corners.size()+id2;
    }

  // Find path between selected corners
  vector<ftEdge*> edges = corners[id1]->getFaceEdges(face.get());
  vector<ftEdge*> path;
  shared_ptr<Vertex> curr_vx = corners[id1];
  size_t ki;
  for (ki=0; ki<edges.size(); ++ki)
    {
      ftEdge *edg = edges[ki];
      shared_ptr<Vertex> other;
      while (true)
	{
	  path.push_back(edg);
	  other = edg->getOtherVertex(curr_vx.get());
	  if (other.get() == vx.get() || other.get() == corners[id2].get() ||
	      other.get() == corners[id1].get())
	    break;
	  vector<ftEdge*> edges2 = other->getFaceEdges(face.get());
	  if (edges2.size() != 2)
	    break;
	  edg = (edges2[0] == edg) ? edges2[1] : edges2[0];
	  curr_vx = other;
	}
      if (other.get() == corners[id2].get())
	break;

      path.clear();
      curr_vx = corners[id1];
    }

  if (path.size() == 0)
    return NULL;  // No edges

  // Perform closest point to the found edges
  Point vx_point = vx->getVertexPoint();
  dist = 1.0e8;  // A huge number
  int idx = -1;
  double min_sc = 1.0e8;
  for (ki=0; ki<path.size(); ++ki)
    {
      double par2, dist2;
      Point point2;
      path[ki]->closestPoint(vx_point, par2, point2, dist2);
      Point tang = path[ki]->tangent(par2);
      tang.normalize();
      Point vec = point2 - vx_point;
      vec.normalize();
      double sc = fabs(tang*vec);
	
      if (sc < min_sc)
	{
	  min_sc = sc;
	  dist = dist2;
	  par = par2; 
	  point = point2;
	  idx = (int)ki;
	}
    }

  return (idx >= 0) ? path[idx] : NULL;
}
//==========================================================================
Point RegularizeUtils::getInVec(shared_ptr<Vertex> vx, 
				shared_ptr<ftSurface> face)
//==========================================================================
{
  vector<ftEdge*> edges = vx->getFaceEdges(face.get());

  // Don't expect more than two edges
  double t1 = edges[0]->parAtVertex(vx.get());
  double t2 = edges[1]->parAtVertex(vx.get());
  Point tan1 = edges[0]->tangent(t1);
  Point tan2 = edges[1]->tangent(t2);
  
  Point par = vx->getFacePar(face.get());
  Point vec = 0.5*(tan1 + tan2);
  Point norm1 = face->normal(par[0], par[1]);
  Point norm2 = norm1.cross(vec);

  return norm2;
 }

//==========================================================================
shared_ptr<ParamCurve> RegularizeUtils::checkStrightParCv(shared_ptr<ftSurface> face,
							  shared_ptr<Vertex> vx1, 
							  shared_ptr<Vertex> vx2,
							  double epsge)
//==========================================================================
{
  vector<ftEdge*> edges1 = vx1->getFaceEdges(face.get());
  vector<ftEdge*> edges2 = vx2->getFaceEdges(face.get());

  // Fetch tangents in the face boundary at the vertices
  // Don't expect more than two edges
  vector<Point> tan(4);
  double t1 = edges1[0]->parAtVertex(vx1.get());
  double t2 = edges1[1]->parAtVertex(vx1.get());
  tan[0] = edges1[0]->tangent(t1);
  tan[1] = edges1[1]->tangent(t2);
  // if (edges1[0]->tMax() - t1 < t1 - edges1[0]->tMin())
  //   tan[0] *= -1;
  // else						
  //   tan[1] *= -1;

  double t3 = edges2[0]->parAtVertex(vx2.get());
  double t4 = edges2[1]->parAtVertex(vx2.get());
  tan[2] = edges2[0]->tangent(t3);
  tan[3] = edges2[1]->tangent(t4);
  // if (edges2[0]->tMax() - t3 < t3 - edges2[0]->tMin())
  //   tan[2] *= -1;
  // else						
  //   tan[3] *= -1;

  // Project into the parameter domain
  Point par1 = vx1->getFacePar(face.get());
  Point par2 = vx2->getFacePar(face.get());

  // Compute partial derivatives in the surface
  shared_ptr<ParamSurface> surf = face->surface();
  vector<Point> sf_der1(3), sf_der2(3);
  surf->point(sf_der1, par1[0], par1[1], 1);
  surf->point(sf_der2, par2[0], par2[1], 1);

  // For each tangent vector describe it as a linear combination of the
  // surface derivatives to find the tangents in the parameter domain
  int ki;
  vector<Point> ptan(4);
  int dim = surf->dimension();
  double coef1, coef2;
  for (ki=0; ki<2; ++ki)
    {
      CoonsPatchGen::blendcoef(&sf_der1[1][0], &sf_der1[2][0], &tan[ki][0], dim, 1, 
			       &coef1, &coef2);
      ptan[ki] = Point(coef1, coef2);
    }
  for (ki=0; ki<2; ++ki)
    {
      CoonsPatchGen::blendcoef(&sf_der2[1][0], &sf_der2[2][0], &tan[2+ki][0], dim, 1, 
			       &coef1, &coef2);
      ptan[2+ki] = Point(coef1, coef2);
    }

  // Vector of stright curve in the parameter domain
  Point vec = par2 - par1;

  // Check if this vector is well within the sector defined by the tangents in the
  // parameter domain
  double ang1 = ptan[0].angle(ptan[1]);
  double ang2 = ptan[0].angle(vec);
  ang2 = std::min(ang2, fabs(M_PI-ang2));
  double ang3 = ptan[1].angle(vec);
  ang3 = std::min(ang3, fabs(M_PI-ang3));

  double ang4 = ptan[2].angle(ptan[3]);
  double ang5 = ptan[2].angle(vec);
  ang5 = std::min(ang5, fabs(M_PI-ang5));
  double ang6 = ptan[3].angle(vec);
  ang6 = std::min(ang6, fabs(M_PI-ang6));
  
  bool make_pcrv = false;
  double fac = 0.9;
  Point d1(0.0, 0.0), d2(0.0, 0.0);
  
  double angtol = 10.0*epsge;
  if (ang1 < angtol && std::max(ang2, ang3) < angtol)
    {
      make_pcrv = true;
      d1[0] = -ptan[0][1];
      d1[1] = ptan[0][0];
      d1.normalize();
      Point tmp = ptan[0];
      tmp.normalize();
     if (tmp*vec < 0)
	tmp *= -1;
      double fac2 = 0.5;
      d1 = (1.0-fac2)*d1 + fac*tmp;
    }
  else if ((ang2 < angtol && vec*ptan[0] < 0.0) || 
	   (ang3 < angtol && vec*ptan[1] > 0.0))
    {
      make_pcrv = true;
      Point tmp = 0.5*(ptan[0]+ptan[1]);
      d1 = Point(-tmp[1], tmp[0]);
      d1.normalize();
      tmp = vec;
      tmp.normalize();
      double fac2 = 0.5;
      d1 = (1.0-fac2)*d1 + fac2*tmp;
    }

  if (ang4 < angtol && std::max(ang5, ang6) < angtol)
    {
      make_pcrv = true;
      d2[0] = -ptan[2][1];
      d2[1] = ptan[2][0];
      d2.normalize();
      Point tmp = ptan[2];
      tmp.normalize();
      if (tmp*vec < 0)
	tmp *= -1;
      double fac2 = 0.5;
      d2 = (1.0-fac2)*d2 + fac2*tmp;
    }
  else if ((ang5 < angtol && vec*ptan[2] > 0.0) || 
	   (ang6 < angtol && vec*ptan[3] < 0.0))
    {
      make_pcrv = true;
      Point tmp = 0.5*(ptan[2]+ptan[3]);
      d2 = Point(-tmp[1], tmp[0]);
      d2.normalize();
      tmp = -vec;
      tmp.normalize();
      double fac2 = 0.5;
      d2 = (1.0-fac2)*d2 + fac*tmp;
    }

  shared_ptr<ParamCurve> pcrv;
  if (make_pcrv && (d1.length() > epsge || d2.length() > epsge))
    {
      // Set length of tangents
      double len_fac = 3.0; //0.1; //5.0; //0.1;
      double len = vec.length();
      if (d1.length() > epsge)
	d1.normalize();
      //d1 *= len_fac*len;
      if (d2.length() > epsge)
	d2.normalize();
      d2 *= -1; //len_fac*len;

      // Prepare for interpolation using sisl
      vector<double> epoint;
      vector<int> ntype;
      epoint.insert(epoint.end(), par1.begin(), par1.end());
      ntype.push_back(1);
      if (d1.length() > epsge)
	{
	  epoint.insert(epoint.end(), d1.begin(), d1.end());
	  ntype.push_back(4);
	}
      epoint.insert(epoint.end(), par2.begin(), par2.end());
      ntype.push_back(1);
      if (d2.length() > epsge)
	{
	  epoint.insert(epoint.end(), d2.begin(), d2.end());
	  ntype.push_back(4);
	}
 

      // Interpolate
      SISLCurve *qc = NULL;
      double *gpar = NULL;
      double endpar;
      int nbpar = 0;
      int status = 0;
      s1356(&epoint[0], (int)ntype.size(), 2, &ntype[0], 0, 0, 1, 3, 0.0,
	    &endpar, &qc, &gpar, &nbpar, &status);

      if (status >= 0)
	pcrv = shared_ptr<ParamCurve>(SISLCurve2Go(qc));

      if (qc) free(qc);
      if (gpar) free(gpar);
    }
  return pcrv;
}

//==========================================================================
shared_ptr<ParamCurve> RegularizeUtils::checkStrightParCv(shared_ptr<ftSurface> face,
							  const Point& pos1, 
							  const Point& pos2,
							  double epsge)
//==========================================================================
{
  // Find end boundary points associated to the input points
  Point dummy_vec;
  double u1, v1, dt1, p1, u2, v2, dt2, p2;
  Point clo1, clo2;
  ftEdgeBase *edge1 = face->closestBoundaryPoint(pos1, dummy_vec, u1, v1, clo1,
						 dt1, p1);
  ftEdgeBase *edge2 = face->closestBoundaryPoint(pos2, dummy_vec, u2, v2, clo2,
						 dt2, p2);

  // Fetch tangents in the face boundary points 
  vector<Point> tan(2);
  tan[0] = edge1->tangent(p1);
  tan[1] = edge2->tangent(p2);

  // Project into the parameter domain
  Point par1(u1, v1);
  Point par2(u2, v2);

  // Compute partial derivatives in the surface
  shared_ptr<ParamSurface> surf = face->surface();
  vector<Point> sf_der1(3), sf_der2(3);
  surf->point(sf_der1, par1[0], par1[1], 1);
  surf->point(sf_der2, par2[0], par2[1], 1);

  // For each tangent vector describe it as a linear combination of the
  // surface derivatives to find the tangents in the parameter domain
  int ki;
  vector<Point> ptan(2);
  int dim = surf->dimension();
  double coef1, coef2;
  CoonsPatchGen::blendcoef(&sf_der1[1][0], &sf_der1[2][0], &tan[0][0], dim, 1, 
			   &coef1, &coef2);
  ptan[0] = Point(coef1, coef2);
  CoonsPatchGen::blendcoef(&sf_der2[1][0], &sf_der2[2][0], &tan[1][0], dim, 1, 
			   &coef1, &coef2);
  ptan[1] = Point(coef1, coef2);

  // Vector of stright curve in the parameter domain
  Point vec = par2 - par1;

  // Check if this vector is well within the sector defined by the tangents in the
  // parameter domain
  double ang1 = ptan[0].angle(vec);
  ang1 = std::min(ang1, fabs(M_PI-ang1));
  vec *= -1;
  double ang2 = ptan[1].angle(vec);
  ang2 = std::min(ang2, fabs(M_PI-ang2));

  
  bool make_pcrv = false;
  double fac = 0.9;
  Point d1(0.0, 0.0), d2(0.0, 0.0);
  double ang_tol = 0.15;
  d1[0] = -ptan[0][1];
  d1[1] = ptan[0][0];
  if (ang1 < ang_tol || d1*(par2-par1) <= 0)
    {
      make_pcrv = true;
      d1.normalize();
      Point tmp = ptan[0];
      tmp.normalize();
      if (tmp*(par2-par1) < 0)
	tmp *= -1;
      double fac2 = 0.5;
      d1 = (1.0-fac2)*d1 + fac2*tmp;
    }
  else
    d1[0] = d1[1] = 0.0;

  d2[0] = -ptan[1][1];
  d2[1] = ptan[1][0];
  if (ang2 < ang_tol || d2*vec <= 0)
    {
      make_pcrv = true;
      d2.normalize();
      Point tmp = ptan[1];
      tmp.normalize();
      if (tmp*vec < 0)
	tmp *= -1;
      double fac2 = 0.5;
      d2 = (1.0-fac2)*d2 + fac2*tmp;
    }
    else
      d2[0] = d2[1] = 0.0;

  shared_ptr<ParamCurve> pcrv;
  if (make_pcrv && (d1.length() > epsge || d2.length() > epsge))
    {
      // Set length of tangents
      double len_fac = 3.0; //10.0; //5.0; //0.1;
      double len = vec.length();
      if (d1.length() > epsge)
	d1.normalize();
      //d1 *= len_fac*len;
      if (d2.length() > epsge)
	d2.normalize();
      d2 *= -1; //len_fac*len;

     // Prepare for interpolation using sisl
      vector<double> epoint;
      vector<int> ntype;
      epoint.insert(epoint.end(), par1.begin(), par1.end());
      ntype.push_back(1);
      if (d1.length() > epsge)
	{
	  epoint.insert(epoint.end(), d1.begin(), d1.end());
	  ntype.push_back(4);
	}
      epoint.insert(epoint.end(), par2.begin(), par2.end());
      ntype.push_back(1);
      if (d2.length() > epsge)
	{
	  epoint.insert(epoint.end(), d2.begin(), d2.end());
	  ntype.push_back(4);
	}
 

      // Interpolate
      SISLCurve *qc = NULL;
      double *gpar = NULL;
      double endpar;
      int nbpar = 0;
      int status = 0;
      s1356(&epoint[0], (int)ntype.size(), 2, &ntype[0], 0, 0, 1, 3, 0.0,
	    &endpar, &qc, &gpar, &nbpar, &status);

      if (status >= 0)
	pcrv = shared_ptr<ParamCurve>(SISLCurve2Go(qc));

      if (qc) free(qc);
      if (gpar) free(gpar);
    }
  return pcrv;
}

//==========================================================================
shared_ptr<ParamCurve> RegularizeUtils::checkStrightParCv(shared_ptr<ftSurface> face,
							  shared_ptr<Vertex> vx1, 
							  const Point& mid,
							  double epsge)
//==========================================================================
{
  vector<ftEdge*> edges1 = vx1->getFaceEdges(face.get());

  // Fetch tangents in the face boundary at the vertices
  // Don't expect more than two edges
  vector<Point> tan(4);
  double t1 = edges1[0]->parAtVertex(vx1.get());
  double t2 = edges1[1]->parAtVertex(vx1.get());
  tan[0] = edges1[0]->tangent(t1);
  tan[1] = edges1[1]->tangent(t2);
  tan[0].normalize();
  tan[1].normalize();
  // if (edges1[0]->tMax() - t1 < t1 - edges1[0]->tMin())
  //   tan[0] *= -1;
  // else						
  //   tan[1] *= -1;

  // Project into the parameter domain
  Point par1 = vx1->getFacePar(face.get());

  // The point will typically lie inside a hole. Thus, we need the underlying
  // surface.
  shared_ptr<ParamSurface> surf = face->surface();
  shared_ptr<BoundedSurface> bd_sf = dynamic_pointer_cast<BoundedSurface,ParamSurface>(surf);
  if (bd_sf.get())
    surf = bd_sf->underlyingSurface();
  Point close;
  double paru, parv, dist;
  surf->closestPoint(mid, paru, parv, close, dist, epsge);
  Point par2(paru, parv);

  // Compute partial derivatives in the surface
  vector<Point> sf_der1(3);
  surf->point(sf_der1, par1[0], par1[1], 1);

  // Describe the tangent vector as a linear combination of the
  // surface derivatives to find the tangents in the parameter domain
  int ki;
  vector<Point> ptan(2);
  int dim = surf->dimension();
  double coef1, coef2;
  for (ki=0; ki<2; ++ki)
    {
      CoonsPatchGen::blendcoef(&sf_der1[1][0], &sf_der1[2][0], &tan[ki][0], dim, 1, 
			       &coef1, &coef2);
      ptan[ki] = Point(coef1, coef2);
    }

  // Vector of stright curve in the parameter domain
  Point vec = par2 - par1;

  // Check if this vector is well within the sector defined by the tangents in the
  // parameter domain
  double ang1 = ptan[0].angle(ptan[1]);
  double ang2 = ptan[0].angle(vec);
  double ang3 = ptan[1].angle(vec);

  
  bool make_pcrv = false;
  double fac = 0.9;
  Point d1(0.0, 0.0);
  double ang_tol = 0.15;
  double angtol = 10.0*epsge;
  if (ang1 < angtol && std::max(ang2, ang3) < angtol)
    {
      make_pcrv = true;
      d1[0] = -ptan[0][1];
      d1[1] = ptan[0][0];
      d1.normalize();
      Point tmp = ptan[0];
      tmp.normalize();
       if (tmp*vec < 0)
	tmp *= -1;
       double fac2 = 0.5;
      d1 = (1.0-fac2)*d1 + fac2*tmp;
   }
  else if ((ang2 < angtol && vec*ptan[0] < 0.0) || 
	   (ang3 < angtol && vec*ptan[1] > 0.0))
    {
      make_pcrv = true;
      Point tmp = 0.5*(ptan[0] + ptan[1]);
      d1[0] = -tmp[1];
      d1[1] = tmp[0];
      d1.normalize();
      tmp = vec;
      tmp.normalize();
      double fac2 = 0.5;
      d1 = (1.0-fac2)*d1 + fac2*tmp;
     }
 
  // Check if the given point lies at a face boundary and must be checked towards
  // the boundary tangent
  Point dummy_vec;
  double u2, v2, dt2, p2;
  Point clo2;
  ftEdgeBase *edg2 = face->closestBoundaryPoint(mid, dummy_vec, u2, v2, clo2,
						dt2, p2);
  Point d2(0.0, 0.0);
  if (dt2 <= epsge)
    {
      Point tan2 = edg2->tangent(p2);
      vector<Point> sf_der2(3);
      surf->point(sf_der2, par2[0], par2[1], 1);

      CoonsPatchGen::blendcoef(&sf_der2[1][0], &sf_der2[2][0], &tan2[0], dim, 1, 
			       &coef1, &coef2);
      Point ptan2 = Point(coef1, coef2);

      double angle = ptan2.angle(vec);

      d2[0] = -ptan2[1];
      d2[1] = ptan2[0];
      if (angle < ang_tol || d2*(par1-par2) <= 0)
	{
	  make_pcrv = true;
	  d2.normalize();
	  Point tmp = ptan2;
	  tmp.normalize();
	  if (tmp*(par1-par2) < 0)
	    tmp *= -1;
	  double fac2 = 0.5;
	  d2 = (1.0-fac2)*d2 + fac2*tmp;
	}
      else
	d2[0] = d2[1] = 0.0;
    }


  shared_ptr<ParamCurve> pcrv;
  if (make_pcrv /*&& d1.length() > epsge*/)
    {
      // Set length of tangent
      double len_fac = 3.0; //10.0; //6.0; //3.0; //0.3;
      double len = vec.length();
      if (d1.length() > epsge)
	d1.normalize();
     //d1 *= len_fac*len;
      if (d2.length() > epsge)
	d2.normalize();
      d2 *= -1; //len_fac*len;
 
      // Prepare for interpolation using sisl
      vector<double> epoint;
      vector<int> ntype;
      epoint.insert(epoint.end(), par1.begin(), par1.end());
      ntype.push_back(1);
      if (d1.length() > epsge)
	{
	  epoint.insert(epoint.end(), d1.begin(), d1.end());
	  ntype.push_back(4);
	}
      epoint.insert(epoint.end(), par2.begin(), par2.end());
      ntype.push_back(1);
      if (d2.length() > epsge)
	{
	  epoint.insert(epoint.end(), d2.begin(), d2.end());
	  ntype.push_back(4);
	}

      // Interpolate
      SISLCurve *qc = NULL;
      double *gpar = NULL;
      double endpar;
      int nbpar = 0;
      int status = 0;
      s1356(&epoint[0], (int)ntype.size(), 2, &ntype[0], 0, 0, 1, 3, 0.0,
	    &endpar, &qc, &gpar, &nbpar, &status);

      if (status >= 0)
	pcrv = shared_ptr<ParamCurve>(SISLCurve2Go(qc));

      if (qc) free(qc);
      if (gpar) free(gpar);
    }
  return pcrv;
}

//==========================================================================
bool RegularizeUtils::checkRegularity(vector<shared_ptr<Vertex> >& cand_vx,
				      shared_ptr<ftSurface> face,
				      bool checkConvex)
//==========================================================================
{
  if (cand_vx.size() != 4)
    return false;

  // For the time being, we apply a simple check that requires the the difference
  // in length between opposite edges to be less than 50% and corner angles to
  // be between pi/4 and 3pi/4.
  vector<Point> pos(4);
  int ki, kj;
  for (ki=0; ki<4; ++ki)
    pos[ki] = cand_vx[ki]->getVertexPoint();

  vector<Point> vec(4);
  for (ki=0; ki<4; ++ki)
    vec[ki] = pos[(ki+1)%4] - pos[ki];

  double frac = 0.1;
  for (ki=0; ki<2; ++ki)
    {
      double l1 = vec[ki].length();
      double l2 = vec[ki+2].length();
      if (std::min(l1,l2) < frac*std::max(l1,l2))
	return false;
      if (vec[ki]*vec[ki+2] >= 0.0)
	return false;
    }

  double level_ang = 0.1*M_PI; //0.25*M_PI;
  for (ki=0; ki<4; ki++)
    {
      kj = (ki+1)%4;
      double ang = vec[ki].angle(vec[kj]);
      if (ang < level_ang || ang > M_PI-level_ang)
	return false;
    }

  if (checkConvex)
    {
      // An extra check for the tangent in the concave vertex
      ftEdge *edge = cand_vx[0]->getCommonEdge(cand_vx[1].get());
      if (!edge)
	return false;

      Point tan = edge->tangent(edge->parAtVertex(cand_vx[0].get()));
      double ang = tan.angle(vec[3]);
      if (ang < level_ang || ang > M_PI-level_ang)
	return false;

      // Make sure that the new edges lies inside the material
      Point invec = getInVec(cand_vx[0], face);
      if (invec*vec[3] > 0.0)
	return false;
    }
  return true;
}

//==========================================================================
vector<shared_ptr<Vertex> > RegularizeUtils::endVxInChain(shared_ptr<ftSurface> face,
							  ftSurface* face1,
							  ftSurface* face2,
							  shared_ptr<Vertex> vx,
							  shared_ptr<Vertex> prev,
							  shared_ptr<Vertex> vx0,
							  vector<shared_ptr<Vertex> >& met_already)
//==========================================================================
{
  vector<shared_ptr<Vertex> > end_vx;

  // Fetch vertex edges not bounding the current face
  vector<ftEdge*> edges = vx->uniqueEdges();
  size_t kr;
  for (kr=0; kr<edges.size(); )
    {
      if (edges[kr]->face() == face1 ||
	  (edges[kr]->twin() && edges[kr]->twin()->geomEdge()->face() == face1))
	edges.erase(edges.begin()+kr);
      else if (face2 && 
	       (edges[kr]->face() == face2 ||
		(edges[kr]->twin() && edges[kr]->twin()->geomEdge()->face() == face2)))
	edges.erase(edges.begin()+kr);
      else 
	kr++;
    }

  // Traverse all candidate edges
  for (kr=0; kr<edges.size(); ++kr)
    {
      shared_ptr<Vertex> vx2 = edges[kr]->getOtherVertex(vx.get());
#ifdef DEBUG_REG
      std::ofstream of("vx2.g2");
      of << "400 1 0 4 0 255 0 255 " << std::endl;
      of << "1" << std::endl;
      of << vx2->getVertexPoint() << std::endl;
#endif

      if (vx2.get() == prev.get())
	continue;  // Do not go back in loop

      if (vx2.get() == vx0.get())
	continue;  // Stop loop 

      if (face1 != face.get() && vx2->hasFace(face.get()))
	{
	  // A candidate endpoint is found
	  end_vx.push_back(vx2);
	  continue;
	}

      // Check if this path is pursued before
      size_t kj;
      for (kj=0; kj<met_already.size(); ++kj)
	if (met_already[kj].get() == vx2.get())
	  break;
      if (kj < met_already.size())
	continue;

      // Fetch faces adjacent to current edge and continue the search
      ftSurface* curr_face1 = edges[kr]->face()->asFtSurface();
      ftSurface* curr_face2 = (edges[kr]->twin()) ? 
	edges[kr]->twin()->face()->asFtSurface() : NULL;

      if (face1 != face.get() && face2 == NULL && curr_face2 == NULL)
	continue;  // Not a good patch

      met_already.push_back(vx2);

      vector<shared_ptr<Vertex> > curr_end_vx = endVxInChain(face, curr_face1, 
							     curr_face2, vx2, 
							     vx, vx0, met_already);
      if (curr_end_vx.size() > 0)
	end_vx.insert(end_vx.end(), curr_end_vx.begin(), curr_end_vx.end());
    }

  return end_vx;
}
      
//==========================================================================
int RegularizeUtils::traverseUntilTJoint(vector<ftSurface*> vx_faces,
					shared_ptr<Vertex> vx,
					shared_ptr<Vertex>& vx2,
					vector<ftSurface*>& vx_faces2)
//==========================================================================
{
  int status = 2;  // Indicate not a smooth transition

  Body *bd0 = vx_faces[0]->getBody();
  shared_ptr<Vertex> vx0 = vx;
  shared_ptr<Vertex> vx1 = vx2;
  size_t ki, kj;
  while (vx_faces2.size() < 3)
    {
      // Get next edge
      vector<ftEdge*> edges2 = vx1->uniqueEdges();
      if (edges2.size() > 2)
	return status;

      for (ki=0; ki<edges2.size(); ++ki)
	{
	  shared_ptr<Vertex> vx3 = edges2[ki]->getOtherVertex(vx1.get());
	  if (vx3.get() == vx0.get())
	    continue;

	  vx2 = vx3;
	  vx_faces2 = vx2->faces();

	  // Remove faces from other bodies
	  for (kj=0; kj<vx_faces2.size();)
	    {
	      if (vx_faces2[kj]->getBody() != bd0)
		vx_faces2.erase(vx_faces2.begin()+kj);
	      else
		kj++;
	    }
	  if (vx_faces2.size() > 3)
	    continue;

	  shared_ptr<FaceConnectivity<ftEdgeBase> > info = 
	    edges2[ki]->getConnectivityInfo();

	  status = info->WorstStatus();
	  if (status > 1)
	    return status;
	}
      vx0 = vx1;
      vx1 = vx2;
    }
  return status;
}

//==========================================================================
void RegularizeUtils::angleInEndpoints(shared_ptr<CurveOnSurface> seg,
				       shared_ptr<Vertex> vx1, 
				       shared_ptr<Vertex> vx2,
				       shared_ptr<ftSurface> face,
				       double& min_ang1, double& min_ang2)
//==========================================================================
{
  vector<ftEdge*> edg1 = vx1->getFaceEdges(face.get());
  vector<ftEdge*> edg2 = vx2->getFaceEdges(face.get());
  vector<Point> der1(2), der2(2);
  seg->point(der1, seg->startparam(), 1);
  seg->point(der2, seg->endparam(), 1);
  Point pos1 = vx1->getVertexPoint();
  Point pos2 = vx2->getVertexPoint();
  if (der1[0].dist(pos1) > der1[0].dist(pos2))
    std::swap(der1, der2);

  // First endpoint
  size_t ki;
  min_ang1 = M_PI;
  for (ki=0; ki<edg1.size(); ++ki)
    {
      double t1 = edg1[ki]->parAtVertex(vx1.get());
      Point tan = edg1[ki]->tangent(t1);
      double ang = der1[1].angle(tan);
      if (fabs(M_PI-ang) < ang)
	ang = fabs(M_PI-ang);
      min_ang1 = std::min(min_ang1, ang);
    }

  // Second endpoint

  min_ang2 = M_PI;
  for (ki=0; ki<edg2.size(); ++ki)
    {
      double t1 = edg2[ki]->parAtVertex(vx2.get());
      Point tan = edg2[ki]->tangent(t1);
      double ang = der2[1].angle(tan);
      if (fabs(M_PI-ang) < ang)
	ang = fabs(M_PI-ang);
      min_ang2 = std::min(min_ang2, ang);
    }
}

//==========================================================================
void RegularizeUtils::getSourceCvs(vector<shared_ptr<ftEdge> >& all_edg,
				   vector<shared_ptr<ParamCurve> >& all_cvs)
//==========================================================================
{
  all_cvs.resize(all_edg.size());

  for (size_t ki=0; ki<all_edg.size(); ++ki)
    {
      shared_ptr<ParamCurve> curr = all_edg[ki]->geomCurve();

      // Make sure to have access to the geometry space curve
      shared_ptr<CurveOnSurface> sf_cv = 
	dynamic_pointer_cast<CurveOnSurface,ParamCurve>(curr);
      if (sf_cv.get())
	curr = sf_cv->spaceCurve();

      // Check if the curve is a spline representation of an elementary curve
      shared_ptr<SplineCurve> spline_cv =
	dynamic_pointer_cast<SplineCurve,ParamCurve>(curr);
      if (spline_cv.get() && spline_cv->isElementaryCurve())
	all_cvs[ki] = spline_cv->getElementaryCurve();
      else
	all_cvs[ki] = curr;
    }
}

