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

#ifndef LRAPPROXAPP_H
#define LRAPPROXAPP_H

#include "GoTools/lrsplines2D/LRSplineSurface.h"

namespace Go
{
  /// Functionality related to approximation of a point cloud by an LR B-spline
  /// surface
  
  namespace LRApproxApp
  {
    /// Approximate point cloud with LR B-spline surface.
    /// Interface to LRSurfApprox with some default choices. The approach
    /// is iterative approximation within a presecribed number of
    /// iterations (max_iter). Refinement is performed in areas where the
    /// accuracy threshold (eps) is not met. Information about the maximum
    /// and average distance in each point and the number of points
    /// being more distant than this threshold is reported back. The function
    /// uses two approximation methods, starting with least squares 
    /// approximation and turning then turning to multilevel B-spline
    /// approximation adapted for LR B-splines depending on the parameter
    /// initmba and tomba. At most 5 iterations will be performed with
    /// least squares approximation
    void pointCloud2Spline(std::vector<double>& points, int dim,
			   double domain[], double reduced_domain[],
			   double eps, int max_iter,
			   shared_ptr<LRSplineSurface>& surf,
			   double& maxdist, double& avdist, 
			   double& avdist_out, int& nmb_out,
			   int mba=0, int initmba=1, int tomba=5);

    /// Approximate point cloud with LR B-spline surface by updating an
    /// initial LR B-spline surface. Interface to LRSurfApprox with some
    /// default choices. The approach
    /// is iterative approximation within a presecribed number of
    /// iterations (max_iter). Refinement is performed in areas where the
    /// accuracy threshold (eps) is not met. Information about the maximum
    /// and average distance in each point and the number of points
    /// being more distant than this threshold is reported back. The function
    /// uses two approximation methods, starting with least squares 
    /// approximation and turning then turning to multilevel B-spline
    /// approximation adapted for LR B-splines depending on the parameter
    /// initmba and tomba. 
    void pointCloud2Spline(std::vector<double>& points, 
			   shared_ptr<LRSplineSurface>& init_surf,
			   std::vector<double>& extent,	    
			   double eps, int max_iter,
			   shared_ptr<LRSplineSurface>& surf,
			   double& maxdist, double& avdist, 
			   double& avdist_out, int& nmb_out,
			   int mba=1, int tomba=0);

    /// Compute point cloud distance with respect to an LR B-spline surface
    void computeDistPointSpline(std::vector<double>& points,
				shared_ptr<LRSplineSurface>& surf,
				double& max_above, double& max_below, 
				double& avdist, int& nmb_points,
				std::vector<double>& pointsdist,
				int use_proj = 0);

    /// Compute point cloud distance with respect to an LR B-spline surface
    /// Multi-threaded version
    void computeDistPointSpline_omp(std::vector<double>& points,
				     shared_ptr<LRSplineSurface>& surf,
				     double& max_above, double& max_below, 
				     double& avdist, int& nmb_points,
				     std::vector<double>& pointsdist,
				    int use_proj = 0);

    /// Compute point cloud distance with respect to an LR B-spline surface
    /// and group points according to this distances
    void classifyCloudFromDist(std::vector<double>& points,
			       shared_ptr<LRSplineSurface>& surf,
			       std::vector<double>& limits,
			       double& max_above, double& max_below, 
			       double& avdist, int& nmb_points,
			       std::vector<std::vector<double> >& level_points,
			       std::vector<int>& nmb_group,
			       int use_proj = 0);

    /// Compute point cloud distance with respect to an LR B-spline surface
    /// and group points according to this distances
    /// Multi-threaded version
    void classifyCloudFromDist_omp(std::vector<double>& points,
				   shared_ptr<LRSplineSurface>& surf,
				   std::vector<double>& limits,
				   double& max_above, double& max_below, 
				   double& avdist, int& nmb_points,
				   std::vector<std::vector<double> >& level_points,
				   std::vector<int>& nmb_group,
				   int use_proj = 0);

    /// Compute point cloud distance with respect to an LR B-spline surface
    /// and classify each point according to this distances
    // classification.size() == points.size()/3
    void categorizeCloudFromDist(std::vector<double>& points,
				 shared_ptr<LRSplineSurface>& surf,
				 std::vector<double>& limits,
				 double& max_above, double& max_below, 
				 double& avdist, int& nmb_points,
				 std::vector<int>& classification,
				 std::vector<int>& nmb_group,
				 int use_proj = 0);

    /// Compute point cloud distance with respect to an LR B-spline surface
    /// and classify each point according to this distances
    /// Multi-threaded version
    void categorizeCloudFromDist_omp(std::vector<double>& points,
				     shared_ptr<LRSplineSurface>& surf,
				     std::vector<double>& limits,
				     double& max_above, double& max_below, 
				     double& avdist, int& nmb_points,
				     std::vector<int>& classification,
				     std::vector<int>& nmb_group,
				     int use_proj = 0);

    /// Compute surface bounding a point cloud approximated by an LR B-spline
    /// surface. The bounding surfaces has the same LR mesh as the initial
    /// surface. Applied to LR B-spline surfaces and bounded surfaces having
    /// an LR B-spline surface as its underlying surface
    void limitingSurfs(std::vector<double>& points,  // The points are modified!!!
		       shared_ptr<ParamSurface>& surf,
		       int nmb_iter,
		       shared_ptr<ParamSurface>& limsf1,
		       shared_ptr<ParamSurface>& limsf2);

   /// Compute surface bounding a point cloud approximated by an LR B-spline
    /// surface. The bounding surfaces has the same LR mesh as the initial
    /// surface. Applied to LR B-spline surfaces
    void limitingSurfs(std::vector<double>& points,  // The points are modified!!!
		       shared_ptr<LRSplineSurface>& surf,
		       int nmb_iter,
		       shared_ptr<LRSplineSurface>& limsf1,
		       shared_ptr<LRSplineSurface>& limsf2);

    /// Update surface with information already maintained in the surface.
    /// Aimed at improving the approximation accuracy with respect to
    /// points classified as significant. A post process to pointCloud2Spline
    /// or similar functionality in LRSurfApprox
    void updateSurfWithSignificantPts(shared_ptr<LRSplineSurface>& surf,
				      double tol, double tol_sign,
				      double fac1, double fac2,
				      double& maxdist, double& avdist,
				      double& avdist_out, int& nmb_out,
				      double& maxsign, double& avsign,
				      int& nmb_out_sign);
  };
};

#endif
