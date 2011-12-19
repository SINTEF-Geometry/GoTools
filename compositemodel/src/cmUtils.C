//===========================================================================
//
// File : cmUtils.C
//
// Created: Thu Feb 21 12:16:45 2008
//
// Author: Kjell Fredrik Pettersen
//
// Revision: $Id: cmUtils.C,v 1.2 2008-04-02 07:18:29 vsk Exp $
//
// Description:
//
//===========================================================================


#include "GoTools/compositemodel/cmUtils.h"
#include "GoTools/compositemodel/ftEdgeBase.h"
#include "GoTools/geometry/BoundedSurface.h"
#include "GoTools/geometry/GeometryTools.h"


using std::vector;


namespace Go
{

//===========================================================================
double cmUtils::estimatedCurveLength(ftEdgeBase* edge, int nmb_samples)
//===========================================================================
{
    ALWAYS_ERROR_IF(nmb_samples < 2,
		"At least 2 points needed to estimate curve length!");

    double total_length = 0.0;
    double tmin = edge->tMin();
    double tmax = edge->tMax();
    double step = (tmax - tmin) / (nmb_samples - 1);
    Point curr_pt = edge->point(tmin);
    Point next_pt = edge->point(tmin + step);
    total_length += (next_pt - curr_pt).length();
    for (int i = 1; i < nmb_samples - 1; ++i) {
	curr_pt = next_pt;
	next_pt = edge->point(tmin + i*step);
	total_length += (next_pt - curr_pt).length();
    }

    return total_length;
}

//===========================================================================
RectDomain cmUtils::geometricParamDomain(ParamSurface* sf)
//===========================================================================
{
  double len_u, len_v;
  estimateSurfaceSize(*sf, len_u, len_v);

  RectDomain domain = sf->containingDomain();

  double u_from = domain.umin();
  double u_to = u_from + len_u;
  double v_from = domain.vmin();
  double v_to = v_from + len_v;
  
  Vector2D ll(u_from, v_from);
  Vector2D ur(u_to, v_to);
  
  return RectDomain(ll, ur);
}

} // namespace Go
