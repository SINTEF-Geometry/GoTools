//===========================================================================
//                                                                           
// File: ftSmoothSurf.h                                                   
//                                                                           
// Created: Mon Mar 19 2001                                         
//                                                                           
// Author: Vibeke Skytt 
//                                                                           
// Revision: 
//                                                                           
// Description: 
//
// Implementation files: ftSmoothSurf.C
//                                                                           
//===========================================================================

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
 
