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

#ifndef _ADAPTSURFACE_H
#define _ADAPTSURFACE_H

#include "GoTools/geometry/BsplineBasis.h"
#include "GoTools/compositemodel/ftPointSet.h"
#include <vector>

namespace Go
{
  class SplineSurface;
  class ParamSurface;
  class SplineCurve;
  class ParamCurve;

/// This namespace contains functions for generating SplineSurfaces by
/// approximation, etc.
  namespace AdaptSurface
  {
    /// Approximate surf1 within the tolerance tol in the spline space of 
    /// surf2 or in a refinement of this spline space.
    shared_ptr<SplineSurface>
      approxInSplineSpace(shared_ptr<ParamSurface> surf,
			  shared_ptr<SplineSurface> surf2,
			  double tol);

    /// Approximate the surfaces surf1 and surf2 in the same spline space
    /// with approximation errors less than tol
    std::vector<shared_ptr<SplineSurface> >
      expressInSameSplineSpace(shared_ptr<ParamSurface> surf1,
			       shared_ptr<ParamSurface> surf2,
			       double tol);

    /// Compute the correspondance between corners in the surfaces
    /// surf1 and surf2. This function is used from approxInSplineSpace
    /// and expressInSameSplineSpace
    bool
      getCornerCorrespondance(shared_ptr<ParamSurface> surf1,
			      shared_ptr<ParamSurface> surf2,
			      int& idx, bool& turned);

    
    /// Approximate a number of curves in the same spline space given
    /// an initial spline space that may be refined. The approximation
    /// tolerance is tol.
    /// This function is an interface to the function
    /// CurveCreators::curveApprox in gotools-core creators
    std::vector<shared_ptr<SplineCurve> > 
      curveApprox(shared_ptr<ParamCurve> cvs[], int nmb_cvs,
		  const BsplineBasis& init_basis, double tol);

    /// Approximate a number of curves in the same spline space. The 
    /// approximation tolerance is tol
    /// This function is an interface to the function
    /// CurveCreators::curveApprox in gotools-core creators
    std::vector<shared_ptr<SplineCurve> > 
      curveApprox(shared_ptr<ParamCurve> cvs[], int nmb_cvs,
		  double tol, double degree=3);

    /// Adapt the initial surface, init_surf, to the surface, surf, within the
    /// given tolerance, tol.
    /// The spline space of the initial surface may be refined
    shared_ptr<SplineSurface> 
      adaptSurface(shared_ptr<ParamSurface> surf, 
		   shared_ptr<SplineSurface> init_surf, double tol);

    /// Parameterize the point set points with respect to the surface init_surf.
    /// To be used in surface approximation.
    double 
      parameterizePoints(shared_ptr<SplineSurface> init_surf,
			 shared_ptr<ftPointSet> points,
			 std::vector<int> corner);

    /// Parameterize by projection
    double projectPoints(shared_ptr<SplineSurface> surf,
		       shared_ptr<ftPointSet> points);

    /// Given an initial spline space and a point set, perform maximum
    /// max_iter iterations to approximate the points within the tolerance
    /// tol. The maximum and average error in the given points are returned.
    shared_ptr<SplineSurface>
      doApprox(shared_ptr<SplineSurface> init_surf, int max_iter,
	       shared_ptr<ftPointSet> points, double tol,
	       double& max_error, double& mean_error);

    /// Fetch data and create triangulation
    void createTriangulation(shared_ptr<ParamSurface> surf, 
			     const RectDomain& dom,
			     shared_ptr<ftPointSet>& points, 
			     vector<int>& corner, 
			     bool consider_joint = true, int nmb=-1);
    
    /// Fetch data points at the boundaries a the surface surf. To be used
    /// in surface approximation
    void
      getBoundaryData(shared_ptr<ParamSurface> surf, 
		      const RectDomain& dom,
		      int nmb_sample, 
		      shared_ptr<ftPointSet> points, 
		      std::vector<int>& corner);

    /// Fetch data points in the inner of the surface surf. To be used
    /// in surface approximation
    void
      getInnerData(shared_ptr<ParamSurface> surf, 
		   const RectDomain& dom,
		   int nmb_sample,
		   shared_ptr<ftPointSet> points,
		   bool consider_joint = true);

    /// Compute point set topology
    void updatePointTopology(shared_ptr<ParamSurface> surf, 
			     ftPointSet& points);
  }
}

#endif    // #ifndef _ADAPTSURFACE_H
