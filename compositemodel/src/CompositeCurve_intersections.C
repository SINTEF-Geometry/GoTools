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

#include "GoTools/utils/Point.h"
#include "GoTools/compositemodel/ftLine.h"
#include "GoTools/geometry/PointOnCurve.h"
#include "GoTools/compositemodel/IntResultsCompCv.h"
#include "GoTools/compositemodel/CompositeCurve.h"
#include "sislP.h"
#include "GoTools/geometry/SISLconversion.h"


using std::vector;

namespace Go
{

//===========================================================================
  /// Intersection with a plane. 
  shared_ptr<IntResultsModel> CompositeCurve::intersect_plane(const ftPlane& plane)
//===========================================================================
  {
    shared_ptr<IntResultsCompCv> intersections = 
      shared_ptr<IntResultsCompCv>(new IntResultsCompCv(this,
						      plane)); // Empty storage f

  // Check if any intersection is possible
  BoundingBox bigbox = boundingBox();
  if (!plane.intersectsBox(bigbox))
    return intersections;

  // For each curve in the composite curve, check if an intersection is possible  and
  // compute eventual intersections
  for (size_t ki=0; ki<curves_.size(); ++ki)
    {
      BoundingBox box = curves_[ki]->boundingBox();
      if (!plane.intersectsBox(box))
	continue;
      
      vector<double> pt_par;  // Parameter values of intersection points
      vector<double> line_par; // Start and end parameter values of intersection line segments
      localIntersect(plane, curves_[ki].get(), pt_par, line_par);

      size_t kj;
      for (kj=0; kj<pt_par.size(); kj++)
	intersections->addIntPt(curves_[ki], &pt_par[kj]);

      for (kj=0; kj<line_par.size(); kj+=2)
	intersections->addIntCv(curves_[ki], &line_par[kj],&line_par[kj+1]);
    }

     return intersections;
  }


//===========================================================================
  void 
  CompositeCurve::localIntersect(const ftPlane& plane, ParamCurve* cv, 
				 vector<double>& pt_par,
				 vector<double>& line_par)
//===========================================================================
  {
    // Get spline curve
    SplineCurve* spl = cv->geometryCurve();
    if (!spl)
      return;   // Cannot intersect

    // Convert to SISL format
    SISLCurve *sislcv = Curve2SISL(*spl, false);
    int dim = spl->dimension();
    double epsco = 1e-15; // Not used
    double epsge = 1e-6;
    int numintpt;  // number of single intersection points
    double* pointpar = 0; // array containing the parameter values of single intersect. pt.
    int numintcr; // number of intersection curves
    SISLIntcurve** intcurves = 0;
    int stat = 0;
    
    Point pnt = plane.point();
    Point normal = plane.normal();

    // Compute plane-curve intersections
    s1850(sislcv, pnt.begin(), normal.begin(), dim, epsco, epsge, &numintpt, &pointpar, 
	  &numintcr, &intcurves, &stat);
    MESSAGE_IF(stat!=0, "s1850 returned code: " << stat);
	
    int ki;
    for (ki=0; ki<numintpt; ++ki)
      pt_par.push_back(pointpar[ki]);

    for (ki=0; ki<numintcr; ++ki)
      {
	line_par.push_back(intcurves[ki]->epar1[0]);
	line_par.push_back(intcurves[ki]->epar1[intcurves[ki]->ipoint-1]);
      }
      
    free(pointpar);
    freeIntcrvlist(intcurves, numintcr);
    freeCurve(sislcv);
  }

//===========================================================================
// Intersection with a line. Expected output is points, probably one point. Curves 
// can occur in special configurations
  shared_ptr<IntResultsModel> 
  CompositeCurve::intersect(const ftLine& line)  // Line class, consist of one point and one direction
			// represented by Go::Point. Just storage, not much content (yet)
//===========================================================================
{
  shared_ptr<IntResultsCompCv> intersections = 
    shared_ptr<IntResultsCompCv>(new IntResultsCompCv(this,
						      line)); // Empty storage for output

  // Check if any intersection is possible
  BoundingBox bigbox = boundingBox();
  if (!line.intersectsBox(bigbox))
    return intersections;

  // In case of a 3D intersection, represent the line as a linear curve
  if (curves_.size() == 0)
    return intersections;
  int dim = curves_[0]->dimension();
  shared_ptr<SplineCurve> line_cv;
  if (dim == 3)
    {
      Point line_pt = line.point();
      Point dir = line.direction();
      double len = (bigbox.high() - bigbox.low()).length();
      Point pt1 = line_pt - 1.5*len*dir;
      Point pt2 = line_pt + 1.5*len*dir;

      line_cv = shared_ptr<SplineCurve>(new SplineCurve(pt1, pt2));
    }
  
  // For each curve in the composite curve, check if an intersection is possible  and
  // compute eventual intersections
  for (size_t ki=0; ki<curves_.size(); ++ki)
    {
      BoundingBox box = curves_[ki]->boundingBox();
      if (!line.intersectsBox(box))
	continue;

      vector<double> pt_par;  // Parameter values of intersection points
      vector<double> line_par; // Start and end parameter values of intersection line segments
      int idx_del;
      if (dim == 2)
	{
	  localIntersect(line, curves_[ki].get(), pt_par, line_par);
	  idx_del = 1;
	}
      else if (dim == 3)
	{
	  localIntersect(curves_[ki].get(), line_cv.get(), pt_par, line_par);
	  idx_del = 2;
	}

      size_t kj;
      for (kj=0; kj<pt_par.size(); kj+=idx_del)
	intersections->addIntPt(curves_[ki], &pt_par[kj]);

      for (kj=0; kj<line_par.size(); kj+=2*idx_del)
	intersections->addIntCv(curves_[ki], &line_par[kj],&line_par[kj+idx_del]);
    }

  return intersections;
}
  
//===========================================================================
bool CompositeCurve::hit(const Point& point, const Point& dir, 
			 PointOnCurve& result) 
//===========================================================================
{
  bool hit = false;

  // Check if any intersection is possible
  BoundingBox bigbox = boundingBox();
  ftLine line(dir, point);  // Represent beam as line
  if (!line.intersectsBox(bigbox))
    return hit;

  Point box_mid = 0.5*(bigbox.low() + bigbox.high());  // Midpoint in the box
  double rad = box_mid.dist(bigbox.low());          // Radius in surronding sphere
  double min_dist = point.dist(box_mid) + rad;   // A long distance

  // In case of a 3D intersection, represent the line as a linear curve
  if (curves_.size() == 0)
    return hit;
  int dim = curves_[0]->dimension();
  shared_ptr<SplineCurve> line_cv;
  if (dim == 3)
    {
      Point line_pt = line.point();
      Point dir = line.direction();
      double len = (bigbox.high() - bigbox.low()).length();
      Point pt1 = line_pt - 1.5*len*dir;
      Point pt2 = line_pt + 1.5*len*dir;

      line_cv = shared_ptr<SplineCurve>(new SplineCurve(pt1, pt2));
    }
  
  // For each curve in the composite curve, check if an intersection is possible  and
  // compute eventual intersections
  for (size_t ki=0; ki<curves_.size(); ++ki)
    {
      BoundingBox box = curves_[ki]->boundingBox();
      if (!line.intersectsBox(box))
	continue;

      Point mid = 0.5*(box.low()+box.high());
      double dist = point.dist(mid) - mid.dist(box.low());
      if (dist > min_dist)
	continue;  // No minimum distance can be found

      vector<double> pt_par;  // Parameter values of intersection points
      vector<double> line_par; // Start and end parameter values of intersection line segments
      if (dim == 2)
	{
	  localIntersect(line, curves_[ki].get(), pt_par, line_par);
	}
      else if (dim == 3)
	{
	  localIntersect(curves_[ki].get(), line_cv.get(), pt_par, line_par);
	}

      // Find closest intersction and update smallest distance
      size_t kd;
      for (kd=0; kd<pt_par.size(); ++kd)
	{
	  hit = true;
	  Point pos = curves_[ki]->point(pt_par[kd]);
	  double dist = point.dist(pos);
	  if (dist < min_dist)
	    {
	      result = PointOnCurve(curves_[ki],pt_par[kd]);
	      min_dist = dist;
	    }
	}

      for (kd=0; kd<line_par.size(); ++kd)
	{
	  hit = true;
	  for (int kr=0; kr<2; ++kr)
	    {
	      Point pos = curves_[ki]->point(line_par[2*kd+kr]);
	      double dist = point.dist(pos);
	      if (dist < min_dist)
		{
		  result = PointOnCurve(curves_[ki],line_par[2*kd+kr]);
		  min_dist = dist;
		}
	    }
	}
    }

  return hit;
}

//===========================================================================
  void 
  CompositeCurve::localIntersect(const ftLine& line, ParamCurve* cv, 
				 vector<double>& pt_par,
				 vector<double>& line_par)
//===========================================================================
  {
    // Get spline curve
    SplineCurve* spl = cv->geometryCurve();
    if (!spl)
      return;   // Cannot intersect

    // Convert to SISL format
    SISLCurve *sislcv = Curve2SISL(*spl, false);
    int dim = spl->dimension();
    double epsco = 1e-15; // Not used
    double epsge = 1e-6;
    int numintpt;  // number of single intersection points
    double* pointpar = 0; // array containing the parameter values of single intersect. pt.
    int numintcr; // number of intersection curves
    SISLIntcurve** intcurves = 0;
    int stat;
    
    // Check dimension
    ASSERT(dim == 2);
    Point pnt = line.point();
    Point dir = line.direction();
    Point normal(dir[1], -dir[0]);

    // Compute line-curve intersections
    s1850(sislcv, pnt.begin(), normal.begin(), dim, epsco, epsge, &numintpt, &pointpar, 
	  &numintcr, &intcurves, &stat);
    MESSAGE_IF(stat!=0, "s1850 returned code: " << stat);
	
    int ki;
    for (ki=0; ki<numintpt; ++ki)
      pt_par.push_back(pointpar[ki]);

    for (ki=0; ki<numintcr; ++ki)
      {
	line_par.push_back(intcurves[ki]->epar1[0]);
	line_par.push_back(intcurves[ki]->epar1[intcurves[ki]->ipoint-1]);
      }
      
    free(pointpar);
    freeIntcrvlist(intcurves, numintcr);
    freeCurve(sislcv);
  }

  void 
  CompositeCurve::localIntersect(ParamCurve *cv1, ParamCurve* cv2, 
				 vector<double>& pt_par,
				 vector<double>& line_par)
  {
    // Get spline curve2
    SplineCurve* spl1 = cv1->geometryCurve();
    SplineCurve* spl2 = cv2->geometryCurve();
    if (!spl1 || !spl2)
      return;   // Cannot intersect

    // Convert to SISL format
    SISLCurve *sislcv1 = Curve2SISL(*spl1, false);
    SISLCurve *sislcv2 = Curve2SISL(*spl2, false);
    double epsco = 1e-15; // Not used
    double epsge = 1e-6;
    int numintpt;  // number of single intersection points
    double *pointpar1 = 0, *pointpar2 = 0; // arrays containing the parameter values of 
                                           // single intersect. pts.
    int numintcr; // number of intersection curves
    SISLIntcurve** intcurves = 0;
    int stat;
    
    // Compute curve-curve intersections
    s1857(sislcv1,sislcv2, epsco, epsge, &numintpt, &pointpar1, &pointpar2, 
	  &numintcr, &intcurves, &stat);
    MESSAGE_IF(stat!=0, "s1857 returned code: " << stat);
	
    int ki;
    for (ki=0; ki<numintpt; ++ki)
      {
	pt_par.push_back(pointpar1[ki]);
	pt_par.push_back(pointpar2[ki]);
      }

    for (ki=0; ki<numintcr; ++ki)
      {
	line_par.push_back(intcurves[ki]->epar1[0]);
	line_par.push_back(intcurves[ki]->epar2[0]);
	line_par.push_back(intcurves[ki]->epar1[intcurves[ki]->ipoint-1]);
	line_par.push_back(intcurves[ki]->epar2[intcurves[ki]->ipoint-1]);
      }
      
    free(pointpar1);
    free(pointpar2);
    freeIntcrvlist(intcurves, numintcr);
    freeCurve(sislcv1);
    freeCurve(sislcv2);
  }

//===========================================================================
void 
CompositeCurve::extremalPoint(Point& dir, Point& ext_pnt, int& idx, 
			      double ext_par[])
//===========================================================================
{
  idx = -1;  // No candidate found yet

  // Check all curves
  for (size_t ki=0; ki<curves_.size(); ++ki)
    {
      if (idx < 0 || boxExtreme(curves_[ki]->boundingBox(), dir, ext_pnt))
	{
	  shared_ptr<SplineCurve> spl(curves_[ki]->geometryCurve());
	  if (!spl.get())
	    continue;  // Cannot compute extremal points

	  SISLCurve *sislcv = Curve2SISL(*spl, false);
	  int dim = spl->dimension();
	  double epsco = 1e-15; // Not used
	  double epsge = 1e-6;
	  int numintpt;  // number of single intersection points
	  double* pointpar = 0; // array containing the parameter values of single intersect. pt.
	  int numintcr; // number of intersection curves
	  SISLIntcurve** intcurves = 0;
	  int stat = 0;

	  s1920(sislcv, dir.begin(), dim, epsco, epsge, &numintpt, 
		&pointpar, &numintcr, &intcurves, &stat);
	  
	  MESSAGE_IF(stat!=0, "s1920 returned code: " << stat); 

	  // Check current extremal points
	  int kj;
	  for (kj=0; kj<numintpt; ++kj)
	    {
	      // Evaluate curve
	      Point pos = spl->ParamCurve::point(pointpar[kj]);
	      if (idx < 0 || pos*dir > ext_pnt*dir)
		{
		  idx = (int)ki;
		  ext_pnt = pos;
		  ext_par[0] = pointpar[kj];
		}
	    }

	  for (kj=0; kj<numintcr; ++kj)
	    {
	      Point pos = spl->ParamCurve::point(intcurves[kj]->epar1[0]);
	      if (idx < 0 || pos*dir > ext_pnt*dir)
		{
		  idx = (int)ki;
		  ext_pnt = pos;
		  ext_par[0] = intcurves[kj]->epar1[0];
		}

	      double *pp = intcurves[kj]->epar1 + intcurves[kj]->ipoint-1;
	      pos = spl->ParamCurve::point(pp[0]);
	      if (idx < 0 || pos*dir > ext_pnt*dir)
		{
		  idx = (int)ki;
		  ext_pnt = pos;
		  ext_par[0] = pp[0];
		}
	    }
	  if (sislcv)
	    freeCurve(sislcv);
	  if (pointpar)
	    free(pointpar);
	  if (intcurves)
	    freeIntcrvlist(intcurves, numintcr);
	}
    }
}




} // namespace Go
