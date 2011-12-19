//===========================================================================
//                                                                           
// File: TesselatorUtils.C                                                     
//                                                                           
// Created: Mars 2009 
//                                                                           
// Author: Vibeke Skytt
//                                                                           
// Revision: $Id: TesselatorUtils.C,v 1.4 2009-06-12 08:56:22 vsk Exp $
//                                                                           
// Description:
//                                                                           
//===========================================================================


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
  estimateSurfaceSize(*sf, len_u, len_v);
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
