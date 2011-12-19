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
    shared_ptr<SplineSurface>
      approxInSplineSpace(shared_ptr<ParamSurface> surf,
			  shared_ptr<SplineSurface> surf2,
			  double tol);

    std::vector<shared_ptr<SplineSurface> >
      expressInSameSplineSpace(shared_ptr<ParamSurface> surf1,
			       shared_ptr<ParamSurface> surf2,
			       double tol);

    bool
      getCornerCorrespondance(shared_ptr<ParamSurface> surf1,
			      shared_ptr<ParamSurface> surf2,
			      int& idx, bool& turned);

    std::vector<shared_ptr<SplineCurve> > 
      curveApprox(shared_ptr<ParamCurve> cvs[], int nmb_cvs,
		  const BsplineBasis& init_basis, double tol);

    std::vector<shared_ptr<SplineCurve> > 
      curveApprox(shared_ptr<ParamCurve> cvs[], int nmb_cvs,
		  double tol);

    shared_ptr<SplineSurface> 
      adaptSurface(shared_ptr<ParamSurface> surf, 
		   shared_ptr<SplineSurface> init_surf, double tol);

    double 
      parameterizePoints(shared_ptr<SplineSurface> init_surf,
			 shared_ptr<ftPointSet> points,
			 std::vector<int> corner);

    shared_ptr<SplineSurface>
      doApprox(shared_ptr<SplineSurface> init_surf, int max_iter,
	       shared_ptr<ftPointSet> points, double tol,
	       double& max_error, double& mean_error);

    void
      getBoundaryData(shared_ptr<ParamSurface> surf, 
		      const RectDomain& dom,
		      int nmb_sample, 
		      shared_ptr<ftPointSet> points, 
		      std::vector<int>& corner);

    void
      getInnerData(shared_ptr<ParamSurface> surf, 
		   const RectDomain& dom,
		   int nmb_sample,
		   shared_ptr<ftPointSet> points);

    void updatePointTopology(shared_ptr<ParamSurface> surf, 
			     ftPointSet& points);
  }
}

#endif    // #ifndef _ADAPTSURFACE_H
