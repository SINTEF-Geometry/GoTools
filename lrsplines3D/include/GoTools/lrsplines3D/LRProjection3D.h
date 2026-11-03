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

#ifndef LR_PROJECTION3D_H
#define LR_PROJECTION3D_H

#include "GoTools/lrsplines3D/LRSplineVolume.h"
#include "GoTools/lrsplines3D/LRBSpline3D.h"
#include "GoTools/utils/Point.h"

namespace Go
{
  /// Update coefficients of LR B-spline surface using quasi interpolation
  namespace LRProjection3D
  {
    double computeCoef(LRSplineVolume *vol, 
		       LRBSpline3D *bspl,
		       int proj_type, bool apply_smooth, double dlim,
		       Point& coef, int num_points = -1);

    void RDataSet(LRSplineVolume *vol, LRBSpline3D *bspl, int nmb_pts,
		  double& rad, std::vector<double>& data,
		  int del, int& nmb, int& nmb_out,
		  double& max_dist, double& av_dist);
    
    bool TPproject(LRBSpline3D *bspl,
		   std::vector<double>& data,
		   int del,
		   double smoothwgt,
		   Point& coef);
    
    void IDWproject(LRBSpline3D *bspl, double rad,
		    std::vector<double>& data,
		    int del,
		    Point& coef);

    bool TriLinProject(LRBSpline3D *bspl, std::vector<double>& data, int del,
		      double rad, Point& coef, bool apply_smooth = false);
    
    bool LinProject(LRBSpline3D *bspl, std::vector<double>& data, int del,
		      double rad, Point& coef, bool apply_smooth = true);
    
   bool TriQuadProject(LRBSpline3D *bspl, std::vector<double>& data, int del,
		       double rad, Point& coef, bool apply_smooth = false);
    
    bool QuadProject(LRBSpline3D *bspl, std::vector<double>& data, int del,
		       double rad, Point& coef, bool apply_smooth = false);
    
    bool CubicProject(LRBSpline3D *bspl, std::vector<double>& data, int del,
		      double rad, Point& coef, bool apply_smooth = false);
    
    bool TriCubicProject(LRBSpline3D *bspl, std::vector<double>& data, int del,
			double rad, Point& coef, bool apply_smooth = false);
    
    
    bool PolynomialProject(int degree1, int degree2, int degree3, 
			   int tot_degree, LRBSpline3D *bspl,
			   std::vector<double>& data, int del, double rad, 
			   Point& coef, bool apply_smooth);

    Point Polynomial2Coef(std::vector<double>& pol, int degree1, int degree2,
			  int degree3, int tot_degree, double domain[],
			  LRBSpline3D *bspl);
    
  }; // end namespace LRProjection3D

}; // end namespace Go

#endif
