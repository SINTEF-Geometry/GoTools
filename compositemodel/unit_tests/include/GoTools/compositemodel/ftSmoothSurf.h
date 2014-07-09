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

#ifndef _FTSMOOTHSURF_H
#define _FTSMOOTHSURF_H

//===========================================================================
//===========================================================================

#include <vector>             // Standard library STL vector
#include <string>             // Standard library string
#include "GoTools/compositemodel/ftPointSet.h"       // Points to approximate
#include "GoTools/geometry/SplineSurface.h"

namespace Go
{

/** ftSmoothSurf - Surface smoothing and approximation. Interface to SmoothSurf in
    gotools-core/creators. An iterative process with parameter iteration and 
    refinement of the initial surface with respect to approximation errors at each
    step.
 * 
 */

class ftSmoothSurf
{
public:

  /// Constructor
  /// \param surf initial surface
  /// \param approxtol tolerance in point approximation
  /// \param approx_orig_tol allowed deviance from the original surface
  /// \param ccw_edge_derivs number of derivatives to keep fixed along the
  /// boundaries, sequence: umin, vmax, umax, vmin
  /// \param maxiter maximum allowed number of iterations in approximation
  /// \param lock_corner_points indicates if the surface corners are fixed
    ftSmoothSurf(shared_ptr<SplineSurface> surf, double approxtol,
		 double approx_orig_tol,
		 std::vector<int> ccw_edge_derivs, int maxiter,
		 bool lock_corner_points = false);

    /// Destructor
    ~ftSmoothSurf();

    /// Set expected continuity of seem in 1. parameter direction
    /// 0 <= k <= 1, C0 and C1 continuity possible
    void setSmoothU(int k);
    /// Set expected continuity of seem in 2. parameter direction
    /// 0 <= k <= 1, C0 and C1 continuity possible
    void setSmoothV(int k);


    /// Express the surface to be modified on a refined knot vector
    void refineSurf(ftPointSet& points, bool reparam = true);

    /// Modify the surface. User may turn off reparametrization of points.
    /// If smoothing was a success, true is returned.
    bool update(ftPointSet& points, double gapeps, bool reparam = true);

    /// Weight on point approximation
    void setApproxWeight(double weight)
      { init_approx_weight_ = weight; }

    /// Fetch weight on point approximation
    double getApproxWeight()
      { return init_approx_weight_; }

    /// Fetch maximum and average error in approximation
    void getError(double& max_error, double& mean_error)
      {
	max_error = max_error_;
	mean_error = mean_error_;
      }

protected:
    shared_ptr<SplineSurface> surf_;
    shared_ptr<SplineSurface> orig_surf_; // We save the original surface.
    double approxtol_;
    double approx_orig_tol_;
    double init_approx_weight_;
    std::vector<int> ccw_edge_derivs_; // size = # edges of surf = 4,
                              // orientation: ccw. # rows to fix.
    double max_error_;
    double mean_error_;
    int maxiter_;
    bool lock_corner_points_;
    int seem_[2];
private:

    // Based on ccw_edge_derivs_ & lock_corner_points_, mark coefs not to be altered.
    std::vector<int> getCoefKnown();

};


} // namespace Go

#endif // _FTSMOOTHSURF_H
 
