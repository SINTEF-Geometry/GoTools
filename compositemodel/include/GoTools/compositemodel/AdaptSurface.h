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
    std::shared_ptr<SplineSurface>
      approxInSplineSpace(std::shared_ptr<ParamSurface> surf,
			  std::shared_ptr<SplineSurface> surf2,
			  double tol);

    std::vector<std::shared_ptr<SplineSurface> >
      expressInSameSplineSpace(std::shared_ptr<ParamSurface> surf1,
			       std::shared_ptr<ParamSurface> surf2,
			       double tol);

    bool
      getCornerCorrespondance(std::shared_ptr<ParamSurface> surf1,
			      std::shared_ptr<ParamSurface> surf2,
			      int& idx, bool& turned);

    std::vector<std::shared_ptr<SplineCurve> > 
      curveApprox(std::shared_ptr<ParamCurve> cvs[], int nmb_cvs,
		  const BsplineBasis& init_basis, double tol);

    std::vector<std::shared_ptr<SplineCurve> > 
      curveApprox(std::shared_ptr<ParamCurve> cvs[], int nmb_cvs,
		  double tol);

    std::shared_ptr<SplineSurface> 
      adaptSurface(std::shared_ptr<ParamSurface> surf, 
		   std::shared_ptr<SplineSurface> init_surf, double tol);

    double 
      parameterizePoints(std::shared_ptr<SplineSurface> init_surf,
			 std::shared_ptr<ftPointSet> points,
			 std::vector<int> corner);

    std::shared_ptr<SplineSurface>
      doApprox(std::shared_ptr<SplineSurface> init_surf, int max_iter,
	       std::shared_ptr<ftPointSet> points, double tol,
	       double& max_error, double& mean_error);

    void
      getBoundaryData(std::shared_ptr<ParamSurface> surf, 
		      const RectDomain& dom,
		      int nmb_sample, 
		      std::shared_ptr<ftPointSet> points, 
		      std::vector<int>& corner);

    void
      getInnerData(std::shared_ptr<ParamSurface> surf, 
		   const RectDomain& dom,
		   int nmb_sample,
		   std::shared_ptr<ftPointSet> points);

    void updatePointTopology(std::shared_ptr<ParamSurface> surf, 
			     ftPointSet& points);
  }
}

#endif    // #ifndef _ADAPTSURFACE_H
