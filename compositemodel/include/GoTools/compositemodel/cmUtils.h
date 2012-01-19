//===========================================================================
//
// File : cmUtils.h
//
// Created: Thu Feb 21 11:32:41 2008
//
// Author: Kjell Fredrik Pettersen
//
// Revision: $Id: cmUtils.h,v 1.2 2008-04-02 07:18:27 vsk Exp $
//
// Description:
//
//===========================================================================


#ifndef __CMUTILS_H
#define __CMUTILS_H


#include "GoTools/compositemodel/ftFaceBase.h"


namespace Go
{

  class ParamSurface;

  /// Various utility functions for the compositemodel module
namespace cmUtils
{

  // Return 2-dimensional point representing parameter values of input point.
  // Point faceParameter(ftEdgeBase* edge, double t);

  /// Estimate the curve length corresponding to an edge
  double estimatedCurveLength(ftEdgeBase* edge, int nmb_samples = 4);

  /// Domain of surface is rescaled according to surface lengths in geometry space.
  RectDomain geometricParamDomain(ParamSurface* sf);


}    // Namespace cmUtils

} // namespace Go

#endif    // #ifndef __CMUTILS_H
