//===========================================================================
//
// File: LinDepUtils.h
//
// Created: Fri Oct 26 12:00:00 2012
//
// Author: Peter Nortoft
//
// Revision: $Id: $
//
// Description: Contains methods for test of linear dependence among
//   LR B-splines.
//
// Revised by: Vibeke Skytt, April 2013. Adapted to GoTools
//
//===========================================================================

#ifndef LINDEP_UTILS_H
#define LINDEP_UTILS_H

#include "GoTools/lrsplines2D/LRSplineSurface.h"

//==============================================================================
namespace Go
//==============================================================================
{
  namespace LinDepUtils
  {
  //============================================================================
  // Tests for potential linear dependence among the LR B-splines of a
  // given LR spline. To be more precise: Tests whether a given LR
  // spline is peelable. We say an LR spline is peelable if the
  // incidence matrix contains NO non-zero entries after peeling it,
  // i.e., if ALL the overloaded LR B-splines are peelable,
  // cf. [Dokken, Lyche & Pettersen, 2012]. The LR spline PEELABILITY
  // is a SUFFICIENT condition for linear INDEPENDENCE of the LR
  // B-splines, or equivalently, the LR spline UNPEELABILITY is a
  // NECESSARY condition for linear DEPENDENCE of the LR B-splines. If
  // an LR spline IS peelable, the LR B-splines are linearly
  // INdependent. If the LR spline is NOT peelable, the LR B-splines
  // MAY be linearly dependent, and further investigations are
  // required to determine whether some of the LR B-splines are
  // ACTUALLY part of a linear dependence relation.
  //============================================================================
  bool isPeelable( const LRSplineSurface& );

  //============================================================================
  // Finds unpeelable overloaded LR B-splines of a given LR spline. An
  // LR B-spline is unpeelable, if it cannot be peeled as described in
  // [Dokken, Lyche & Pettersen, 2012]. These unpeeplable LR B-splines
  // corresponds to rows in the incidence matrix with non-zero
  // entries.
  //============================================================================
  std::vector<LRBSpline2D*> unpeelableBasisFunctions ( const LRSplineSurface& );

  } // end namespace LinDepUtils

} // end namespace Go

#endif
