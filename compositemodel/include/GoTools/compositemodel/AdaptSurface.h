//===========================================================================
//
// File : AdaptSurface.h
//
// Created: September 2010
//
// Author: Vibeke Skytt
//
// Revision: 
//
// Description:
//
//===========================================================================


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
    /// tolerance is tol
    std::vector<shared_ptr<SplineCurve> > 
      curveApprox(shared_ptr<ParamCurve> cvs[], int nmb_cvs,
		  const BsplineBasis& init_basis, double tol);

    /// Approximate a number of curves in the same spline space. The 
    /// approximation tolerance is tol
    std::vector<shared_ptr<SplineCurve> > 
      curveApprox(shared_ptr<ParamCurve> cvs[], int nmb_cvs,
		  double tol);

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

    /// Given an initial spline space and a point set, perform maximum
    /// max_iter iterations to approximate the points within the tolerance
    /// tol. The maximum and average error in the given points are returned.
    shared_ptr<SplineSurface>
      doApprox(shared_ptr<SplineSurface> init_surf, int max_iter,
	       shared_ptr<ftPointSet> points, double tol,
	       double& max_error, double& mean_error);

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
		   shared_ptr<ftPointSet> points);

    /// Compute point set topology
    void updatePointTopology(shared_ptr<ParamSurface> surf, 
			     ftPointSet& points);
  }
}

#endif    // #ifndef _ADAPTSURFACE_H
