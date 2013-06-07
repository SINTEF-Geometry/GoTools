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

#include "GoTools/tesselator/TesselatorUtils.h"
#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/geometry/BoundedSurface.h"
#include "GoTools/geometry/RectDomain.h"
#include "GoTools/geometry/GeometryTools.h"

using namespace Go;
using std::vector;

//===========================================================================
void TesselatorUtils::getResolution(const ParamSurface *surf, 
				    int& u_nmb, int& v_nmb, 
				    int uv_nmb)
//===========================================================================
{
  int min_nmb = 8;

  // Check if the surface is trimmed and trimmed along iso curves. In that
  // case the size of the surface is computed from the smallest possible
  // underlying surface

  double tol2d = 1.0e-4;
  const ParamSurface *sf = surf;
  vector<shared_ptr<ParamSurface> > sfs;
  if (surf->instanceType() == Class_BoundedSurface)
    {
      const BoundedSurface *bd_surf
	  = dynamic_cast<const BoundedSurface*>(surf);
      if (bd_surf->isIsoTrimmed(tol2d))
	{
	  RectDomain domain = bd_surf->containingDomain();
	  RectDomain dom2 = sf->containingDomain();
	  double umin = std::max(domain.umin(), dom2.umin());
	  double umax = std::min(domain.umax(), dom2.umax());
	  double vmin = std::max(domain.vmin(), dom2.vmin());
	  double vmax = std::min(domain.vmax(), dom2.vmax());
    
	  sfs = sf->subSurfaces(umin, vmin, umax, vmax);
	  if (sfs.size() == 1)
	    sf = sfs[0].get();
	}
    }

  // Estimate the size of the surface in the two parameter directions
  double len_u, len_v;
  GeometryTools::estimateSurfaceSize(*sf, len_u, len_v);
  double fac = len_u/len_v;
  double len = sqrt((double)uv_nmb/fac);
  u_nmb = std::max(min_nmb, (int)(fac*len));
  v_nmb = std::max(min_nmb, (int)len);

}

//===========================================================================
shared_ptr<LineCloud> TesselatorUtils::getCtrPol(GeomObject* obj)
//===========================================================================
{
    int dim = obj->dimension();
    ASSERT(dim == 3);
    int ki;
    vector<double> lines;
    if (obj->instanceType() == Class_SplineSurface ||
	obj->instanceType() == Class_BoundedSurface) {
	SplineSurface* spline_sf =
	    (obj->instanceType() == Class_SplineSurface) ?
	    dynamic_cast<SplineSurface*>(obj) :
	    dynamic_cast<SplineSurface*>
	    ((dynamic_cast<BoundedSurface*>(obj))->
	     underlyingSurface().get());
	if (spline_sf != 0) {
	    int kj;
	    int in1 = spline_sf->numCoefs_u();
	    int in2 = spline_sf->numCoefs_v();
	    vector<double>::iterator iter = spline_sf->coefs_begin();
	    for (ki = 0; ki < in1; ++ki) {
		for (kj = 0; kj < in2; ++kj) {
		    if (ki < in1 - 1) {
			lines.insert(lines.end(), iter + (kj*in1 + ki)*dim,
				     iter + (kj*in1 + ki + 1)*dim);
			lines.insert(lines.end(), iter + (kj*in1 + ki + 1)*dim,
				     iter + (kj*in1 + ki + 2)*dim);
		    }
		    if (kj < in2 - 1) {
			lines.insert(lines.end(), iter + (kj*in1 + ki)*dim,
				     iter + (kj*in1 + ki + 1)*dim);
			lines.insert(lines.end(),
				     iter + ((kj + 1)*in1 + ki)*dim,
				     iter + ((kj + 1)*in1 + ki + 1)*dim);
		    }
		}
	    }
	}
    } else if (obj->instanceType() == Class_SplineCurve ||
	       obj->instanceType() == Class_CurveOnSurface) {
      SplineCurve* spline_cv = dynamic_cast<ParamCurve*>(obj)->geometryCurve();

      int in = spline_cv->numCoefs();
      vector<double>::iterator iter = spline_cv->coefs_begin();
      lines.insert(lines.end(), iter, iter + dim);
      for (ki = 0; ki < in - 1; ++ki) {
	lines.insert(lines.end(), iter + ki*dim, iter + (ki + 1)*dim);
	lines.insert(lines.end(), iter + ki*dim, iter + (ki + 1)*dim);
      }
      lines.insert(lines.end(), iter + ki*dim, iter + (ki + 1)*dim);
    }


    int nmb_lines = (int)(lines.size())/(2*dim);
    shared_ptr<LineCloud> line_cloud(new LineCloud(lines.begin(), nmb_lines));

    return line_cloud;
}
