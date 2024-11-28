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

#ifndef _REVENGUTILS_H
#define _REVENGUTILS_H

#include "GoTools/utils/Point.h"
#include "GoTools/compositemodel/RevEngPoint.h"
#include "GoTools/geometry/ParamSurface.h"
#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/geometry/SplineCurve.h"
#include "GoTools/utils/BoundingBox.h"
#include <vector>

namespace Go {

  class RevEngPoint;
  class Line;
  class Plane;
  class Cylinder;
  class Torus;
  class Sphere;
  class Cone;
  
  /** RevEngUtils -  Utility functionality for everse engineering.
   * 
   */
  
  namespace RevEngUtils
  {
    /// Principal analysis based on a group of points
    void principalAnalysis(Point& curr, std::vector<Point>& points, 
			   double lambda[3], double eigenvec[3][3]);

    void principalAnalysis(std::vector<RevEngPoint*>& points, 
			   double lambda[3], double eigenvec[3][3]);

    /// Rotate point group according to specified coordinate system (x=vec1, y=vec2),
    /// and approximate points with a bicubic function. Compute surface normal and
    /// principal curvatures
    void computeLocFunc(Point& curr, std::vector<Point>& points,
			Point& vec1, Point& vec2, Point& normal, Point& mincvec,
			double& minc, Point& maxcvec, double& maxc,
			double& currdist, double& avdist);

    /// Interface to the class SmoothSurf
    void smoothSurf(shared_ptr<SplineSurface>& surf, int fixed);

    /// Approximate parameterized data points given information of polynomial degrees
    /// and initial number of coefficients. Iterate until the tolerance tol or the
    /// maximum number of iterations (max_iter) is reached. Parameter iteration is performed.
    /// Output of accuracy (maximum distance, average distance and number of points outside
    /// the tolerance) and final parameter pairs associated to the points
    shared_ptr<SplineSurface> surfApprox(std::vector<double>& data, int dim,
					 std::vector<double>& param, int order1,
					 int order2, int nmb_coef1, int nmb_coef2,
					 bool close1, bool close2,
					 int max_iter, double tol, double& maxd, 
					 double& avd, int& num_out,
					 std::vector<double>& parvals,
					 double del=0.0);
    
    /// Surface approximation given parameterized data points and information to create
    /// spline space. Limits of parameter intervals is given by the parameterized points
    shared_ptr<SplineSurface> surfApprox(std::vector<double>& data, int dim,
					 std::vector<double>& param, int order1,
					 int order2, int nmb_coef1, int nmb_coef2,
					 double del=0.0);

    /// Surface approximation given parameterized data points and information to create
    /// spline space. Limits of parameter intervals is given as input
    shared_ptr<SplineSurface> surfApprox(std::vector<double>& data, int dim,
					 std::vector<double>& param, int order1,
					 int order2, int nmb_coef1, int nmb_coef2,
					 double umin, double umax, double vmin,
					 double vmax);

    /// Paramerize RevEngPoints by projecting them onto a given plane (spanned by vec1 and
    /// vec2). Extracted xyz values in data and parameters (uv) in param
    void parameterizeWithPlane(std::vector<RevEngPoint*>& pnts, const BoundingBox& bbox,
			       const Point& vec1, const Point& vec2,
			       std::vector<double>& data, std::vector<double>& param);
    
    /// Paramerize points by projecting them onto a given plane (spanned by vec1 and
    /// vec2). xyz values in data and parameters (uv) in param
    void parameterizeWithPlane(std::vector<Point>& pnts, const BoundingBox& bbox,
			       const Point& vec1, const Point& vec2,
			       std::vector<double>& data, std::vector<double>& param);

    /// Parameterize RevEngPoints on the surface surf. Returns information to specify
    /// a spline space in which the points can be approximated (inner1, inner2, close1, close2)
    bool parameterizeOnPrimary(std::vector<RevEngPoint*>& points,
			       shared_ptr<ParamSurface> surf,
			       std::vector<double>& data, 
			       std::vector<double>& param,
			       int& inner1, int& inner2, bool& close1, bool& close2);

    /// Cylinder axis and x- and y-vectors from a point group
    void computeAxis(std::vector<Point>& points,
		     Point& axis, Point& Cx, Point& Cy);

    /// Cylinder axis and x- and y-vectors from a number of RevEngPoint groups
    void computeAxis(std::vector<std::pair<std::vector<RevEngPoint*>::iterator,
		     std::vector<RevEngPoint*>::iterator> >& points,
		     Point& axis, Point& Cx, Point& Cy);

    /// Sphere centre and radius from a number of RevEngPoint groups
    void
    computeSphereProp(std::vector<std::pair<std::vector<RevEngPoint*>::iterator,
		      std::vector<RevEngPoint*>::iterator> >& points,
		      Point& centre, double& radius);

    /// Cone axis and x- and y-vectors from a number of RevEngPoint groups. If the
    /// point group corresponds to a small cone sector, the function for computing
    /// a cylinder axis may be more appropriate
    void coneAxis(std::vector<std::pair<std::vector<RevEngPoint*>::iterator,
		     std::vector<RevEngPoint*>::iterator> >& points,
		     Point& axis, Point& Cx, Point& Cy);

    /// Compute cone apex and opening angle
    void coneApex(std::vector<std::pair<std::vector<RevEngPoint*>::iterator,
		  std::vector<RevEngPoint*>::iterator> >& points,
		  Point axis, Point& apex, double& phi);

    /// Compute cylinder position and radius given point and axis
    void computeCylPosRadius(std::vector<std::pair<std::vector<RevEngPoint*>::iterator,
			     std::vector<RevEngPoint*>::iterator> >& points,
			     Point& low, Point& high, Point& axis, Point& Cx, 
			     Point& Cy, Point& pos, double& radius);

    /// Given 3D points in a plane given with normal and x- and y-axis, compute
    /// circle center and radius
    void computeCircPosRadius(std::vector<Point>& points,
			      const Point& axis, const Point& Cx, const Point& Cy,
			      Point& pos, double& radius);
    
    /// As above, given 3D points in a plane given with normal and x- and y-axis, compute
    /// radius
    void computeRadius(std::vector<Point>& points,
		       Point& axis, Point& Cx, Point& Cy, double& radius);

    /// Compute line approximating an unordered group of points
    void computeLine(vector<Point>& points, Point& pos, Point& dir);
    
    /// Approximate a group of points with a plane  using implicit approximation.
    /// An initial surface normal guess must be provided. Returns point in plane and
    /// surface normal
    void computePlane(std::vector<Point>& points, Point normal, Point mainaxis[3],
		      Point& pos, Point& norm, Point& Cx, Point& Cy);

    /// Rotate a number of RevEngPoint groups around the axis axis with point mid and 
    /// return the resulting points in the plane spanned by xvec and axis
    void rotateToPlane(std::vector<std::pair<std::vector<RevEngPoint*>::iterator,
		       std::vector<RevEngPoint*>::iterator> >& points,
		       Point& xvec, Point& axis, Point& mid, std::vector<Point>& rotated);

    /// Project a number of RevEngPoint groups into the plane represented by the point mid and the
    /// normal axis and return the resulting points along with information about the distance
    /// between the initial points and the plane
    void projectToPlane(std::vector<std::pair<std::vector<RevEngPoint*>::iterator,
			std::vector<RevEngPoint*>::iterator> >& points,
			Point& axis, Point& mid, std::vector<Point>& projected,
			double& maxdist, double& avdist, double dlen=-1.0);

    /// Project a number of RevEngPoints into the plane represented by the point mid and the
    /// normal axis and return the resulting points along with information about the distance
    /// between the initial points and the plane
    void projectToPlane(std::vector<RevEngPoint*>& points,
			Point& axis, Point& mid, std::vector<Point>& projected,
			double& maxdist, double& avdist, double dlen=-1.0);

    /// Rotate points around the axis axis width point mid and return the resulting
    /// points in the plane spanned by xvec and axis
     void rotateToPlane(std::vector<Point>& points,
			Point& xvec, Point& axis, Point& mid, std::vector<Point>& rotated);

    /// Compute distance between RevEngPoints and a surface, and summarize results
    /// \param start Start iterator to vector of points
    /// \param end Iterator to end of vector
    /// \param surf Surface with whom the points will be compared
    /// \param tol Approximation tolerance
    /// \param maxdist Maximum distance between the point cloud and the surface
    /// \param avdist Average absolute distance between the point cloud and the surface
    /// \param inside Number of points where the distance is less than tol and the angle
    /// between the surface normal and the normal in the point (least of triangle normal and
    /// local function normal) is less than angtol if given
    /// \param inside2 Number of points where the distance is less than tol
    /// \param in Points being inside
    /// \param out Points not being inside
    /// \param parvals Parameter values of points projected onto the surface (u1, v1, u2, v2, ...)
    /// \param distang Distance and angle difference between normal vectors for each point
   void distToSurf(std::vector<RevEngPoint*>::iterator start,
		    std::vector<RevEngPoint*>::iterator end,
		    shared_ptr<ParamSurface> surf, double tol,
		   double& maxdist, double& avdist,
		   int& inside, int& inside2,
		   std::vector<RevEngPoint*>& in,
		   std::vector<RevEngPoint*>& out,
		   std::vector<double>& parvals,
		   std::vector<std::pair<double,double> >& distang,
		   double angtol=-1.0);

    /// As above. Slightly simplified interface
   void distToSurf(std::vector<RevEngPoint*>::iterator start,
		   std::vector<RevEngPoint*>::iterator end,
		   shared_ptr<ParamSurface> surf, double tol,
		   double& maxdist, double& avdist,
		   int& inside, int& inside2,
		   std::vector<double>& parvals,
		   std::vector<std::pair<double,double> >& distang,
		   double angtol=-1.0);

    /// As above. Simplified interface
    void distToSurf(std::vector<RevEngPoint*>& points,
		    shared_ptr<ParamSurface> surf, double tol,
		    double angtol, double& maxdist, double& avdist, 
		    int& inside, int& inside2,
		    std::vector<std::pair<double,double> >& dist_ang);
    
    /// As above. Taking points as input. Simplified interface. Only distance
    /// between points and surface is considered (normals not available)
    void distToSurf(std::vector<Point>& points,
		    shared_ptr<ParamSurface> surf, double tol,
		    double& maxdist, double& avdist, int& inside,
		    std::vector<double>& distance);

    /// Distance between points and curve.
    /// \param points Input points
    /// \param curve Check distance to the curve
    /// \param tol Approximation tolerance
    /// \param maxdist Maximum distance between the point cloud and the curve
    /// \param avdist Average absolute distance between the point cloud and the curve
    /// \param inside Number of points where the distance is less than tol 
    /// \param parvals Parameter values of points projected onto the curve (t1, t2, ...)
    /// \param dist Distance for each point
    void distToCurve(std::vector<Point>& points,
		     shared_ptr<ParamCurve> curve, double tol,
		     double& maxdist, double& avdist, int& inside,
		     std::vector<double>& parvals, std::vector<double>& dist);
    
    /// Distance between points and curve. Slightly simplified intervace
     void distToCurve(std::vector<Point>& points,
		     shared_ptr<ParamCurve> curve, double tol,
		     double& maxdist, double& avdist, int& inside,
		     std::vector<double>& dist);
    

    /// Distance between points and curve. Simplified intervace
    void distToCurve(std::vector<Point>& points,
		     shared_ptr<ParamCurve> curve, double tol,
		     double& maxdist, double& avdist, int& inside);

    /// Modify given primary surface (sf_in) to fit in the coordinate system given
    /// by mainaxis
    /// \param points Point group associated with the surface
    /// \diag Used to bound the surface
    shared_ptr<ElementarySurface> elemsurfWithAxis(shared_ptr<ElementarySurface> sf_in,
						 std::vector<RevEngPoint*>& points,
						   Point mainaxis[3], double diag);
    
    /// Compute plane with surface normal axis to best fit the point group points
    /// \param init_loc Initial point in plane
    /// \param mainaxis Used to define the plane parameterization
   shared_ptr<Plane> planeWithAxis(std::vector<RevEngPoint*>& points,
				    Point axis, Point init_loc,
				    Point mainaxis[3]);
    
    /// Compute cylinder approximating the RevEngPoints points. The cylinder axis
    /// axis is input. mainaxis is used to compute x- and y-vectors. The points low and high
    /// are used to assist locating the cylinder point
    shared_ptr<Cylinder> cylinderWithAxis(std::vector<RevEngPoint*>& points,
					  Point axis, Point low, 
					  Point high, Point mainaxis[3]);
    
    /// Compute cylinder approximating the RevEngPoints points. The cylinder axis,
    /// x-vector and point is given as input. The cylinder radius is computed
    shared_ptr<Cylinder> cylinderWithAxis(std::vector<RevEngPoint*>& points,
					  Point axis, Point Cx, Point pos);
    
    /// Compute torus approximating the RevEngPoints points. The torus axis,
    /// and mid point is given as input. Minor and major radius is computed. mainaxix
    /// is used to define x- and y-vector
    shared_ptr<Torus> torusWithAxis(std::vector<RevEngPoint*>& points,
				    Point axis, Point loc, 
				    Point mainaxis[3]);
    
    /// Compute torus approximating the RevEngPoints points. The torus axis,
    /// x-vector and mid point is given as input. Minor and major radius is computed. 
    shared_ptr<Torus> torusWithAxis(std::vector<RevEngPoint*>& points,
				    Point axis, Point Cx, Point pos);
    
    /// Approximate sphere from points. The sphere axis is given as input and the remaining
    /// coordinate axes are computed from mainaxis
    shared_ptr<Sphere> sphereWithAxis(std::vector<RevEngPoint*>& points,
				      Point axis, 
				      Point mainaxis[3]);

    /// Compute sphere radius from points. Center and coordinate axes are given
    shared_ptr<Sphere> sphereWithAxis(std::vector<RevEngPoint*>& points,
				      Point axis, Point Cx, Point pos);
    
    /// Compute cone approximating the RevEngPoints points. The cone axis
    /// axis is input. mainaxis is used to compute x- and y-vectors. The points low and high
    /// are used to assist locating a point on the axis
     shared_ptr<Cone> coneWithAxis(vector<RevEngPoint*>& points,
				  Point axis, Point low, 
				  Point high, Point mainaxis[3]);
    
    /// Compute cone approximating the RevEngPoints points. The  axis,
    /// x-vector and point on axis is given as input. The cone radius at the point and
    /// the cone opening angle is computed
    shared_ptr<Cone> coneWithAxis(std::vector<RevEngPoint*>& points,
				  Point axis, Point Cx, Point pos,
				  double len);

    /// Compute plane approximating a number of RevEngPoint groups (points)
    /// \param bbox BoundingBox containing all points
    /// \param nmbpts Number of points for each point group
    /// \param set_bound Whether the plane should be bounded (from bbox information)
    shared_ptr<ParamSurface> doMergePlanes(std::vector<std::pair<std::vector<RevEngPoint*>::iterator,
					   std::vector<RevEngPoint*>::iterator> > points,
					   const BoundingBox& bbox,
					   std::vector<int>& nmbpts,
					   bool set_bound = true);
    /// Compute cylinder approximating a number of RevEngPoint groups (points)
    /// Parameters as for doMergePlanes
    shared_ptr<ParamSurface> doMergeCylinders(std::vector<std::pair<std::vector<RevEngPoint*>::iterator,
					      std::vector<RevEngPoint*>::iterator> > points,
					      const BoundingBox& bbox,
					      std::vector<int>& nmbpts,
					      bool set_bound = true);
    /// Compute sphere approximating a number of RevEngPoint groups (points)
    /// Parameters as for doMergePlanes
    /// \param normal Axis through sphere poles. The sphere  is initially bounded
     shared_ptr<ParamSurface> doMergeSpheres(std::vector<std::pair<std::vector<RevEngPoint*>::iterator,
					    std::vector<RevEngPoint*>::iterator> > points,
					    const BoundingBox& bbox,
					    std::vector<int>& nmbpts, Point& normal);
    
    /// Compute torus approximating a number of RevEngPoint groups (points)
    /// Parameters as for doMergePlanes. The torus is initially bounded
     shared_ptr<ParamSurface> doMergeTorus(std::vector<std::pair<std::vector<RevEngPoint*>::iterator,
					  std::vector<RevEngPoint*>::iterator> > points,
					  const BoundingBox& bbox,
					  std::vector<int>& nmbpts);

    /// Middle curve between two spline curves
    shared_ptr<SplineCurve> midCurve(shared_ptr<SplineCurve>& cv1,
				     shared_ptr<SplineCurve>& cv2);

    /// Approximate points with spline curve. The points are parameterized on the curve
    /// cvin. ik=degree+1, in=number of coefficients
    void  curveApprox(std::vector<Point>& points,
		      shared_ptr<ParamCurve> cvin,
		      int ik, int in, 
		      shared_ptr<SplineCurve>& curve);
    
    /// Approximate parameterized points with spline curve. ik=degree+1, in=number of coefficients
    void  curveApprox(std::vector<Point>& points,
		      std::vector<double>& param,
		      int ik, int in, 
		      shared_ptr<SplineCurve>& curve);

    /// Approximate an ordered sequence of RevEngPoints by a spline curve with degree degree
    /// Refinement of the intial spline space and reparameterization is performed at most maxiter
    /// times
    /// \param tol Approximation tolerance
    shared_ptr<SplineCurve> createCurve(std::vector<RevEngPoint*>& points, int degree,
					double tol, int maxiter);

    /// Extract points at the start or end of a point sequence where corresponding
    /// points rotated into a given plane (input) can be approximated by a line.
    /// Special function accessed from RevEngRegion
    void extractLinearPoints(std::vector<RevEngPoint*>& points, 
			     std::vector<Point>& rotated, double len,
			     Point& pos, Point& axis, double rad,
			     Point& axis2, bool plane,
			     double tol, double angtol,
			     std::vector<std::pair<double,double> >& dist_ang,
			     std::vector<RevEngPoint*>& linear, bool start,
			     std::vector<RevEngPoint*>& remaining);

    bool extractLinearPoints(std::vector<Point>& points,
			     shared_ptr<Line>& line, Point norm, double tol, 
			     bool start, double& splitpar,
			     std::vector<Point>& linear, 
			     std::vector<Point>& remaining);
    
    void identifyEndPoints(std::vector<RevEngPoint*> edge_pts, shared_ptr<CurveOnSurface>& sfcv,
			   RevEngPoint*& first_pt, double& t1, RevEngPoint*& last_pt, double& t2);

    /// Reorder curves to form a sequence. The parameter direction of curves may be changed
    void setLoopSeq(std::vector<shared_ptr<CurveOnSurface> >& cvs);

    /// Identifyed groups of connected RevEngPoints (by triangle edges) in a given point group
    void identifyConGroups(std::vector<RevEngPoint*>& init,
			   std::vector<std::vector<RevEngPoint*> >& groups);
  }
  
} // namespace Go


#endif // _REVENGUTILS_H

