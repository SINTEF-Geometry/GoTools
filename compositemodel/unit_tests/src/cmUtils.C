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
  GeometryTools::estimateSurfaceSize(*sf, len_u, len_v);

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
