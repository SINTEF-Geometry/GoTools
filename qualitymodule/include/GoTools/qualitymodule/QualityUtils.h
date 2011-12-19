//===========================================================================
//
// File : QualityUtils.h
//
// Created: Mon Jul 28 12:43:28 2008
//
// Author: Kjell Fredrik Pettersen
//
// Revision: $Id: QualityUtils.h,v 1.4 2009-03-20 07:51:53 vsk Exp $
//
// Description:
//
//===========================================================================


#ifndef __CMUTILS_H
#define __CMUTILS_H


#include "GoTools/geometry/ParamSurface.h"
#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/geometry/BoundedSurface.h"
#include "GoTools/geometry/CurveLoop.h"


namespace Go
{

namespace qualityUtils
{

  // Determine whether a surface is a sliver face or not
  // Uses different method for SplineSurface and BoundedSurface
  // A surface is a sliver suface if its maximum length m1 in one
  // parameter direction is less than a thickness given as a parameter,
  // while its minimum length in the other parameter direction is
  // more than m1 * a factor, typically 2.0
  bool isSliverFace(shared_ptr<ParamSurface>,
		    double thickness,
		    double factor = 2.0);

  bool isSliverFace(const SplineSurface& sf,
		    double thickness,
		    double factor = 2.0);

  bool isSliverFace(const BoundedSurface& sf,
		    double thickness,
		    double factor = 2.0);

  bool isSliverFace2(const BoundedSurface& sf,
		     double thickness,
		     double factor = 2.0);

  bool hasIndistinctKnots(shared_ptr<ParamSurface> surf, double tol,
			  std::vector<shared_ptr<ParamCurve> >& trim_cv_knots);

  double estimateArea(shared_ptr<ParamSurface> surf);

  double estimateLoopArea(shared_ptr<CurveLoop> loop);

}   // namespace qualityUtils

}   // namespace Go


#endif    // #ifndef __CMUTILS_H
