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
#include "GoTools/compositemodel/SurfaceModelUtils.h"
#include "GoTools/intersections/Identity.h"
#include "GoTools/geometry/BoundedSurface.h"
#include "GoTools/geometry/ParamSurface.h"
#include "GoTools/geometry/ParamCurve.h"
#include "GoTools/geometry/RectDomain.h"

#include <fstream>

//#define DEBUG

using std::vector;

namespace Go
{
//===========================================================================
vector<shared_ptr<ParamSurface> > 
SurfaceModelUtils::checkClosedFaces(shared_ptr<ParamSurface> surface, double tol)
//===========================================================================
  {
    vector<shared_ptr<ParamSurface> > sfs;

#ifdef DEBUG
    std::ofstream of("close_sf.g2");
    surface->writeStandardHeader(of);
    surface->write(of);
#endif

    // Fetch non-trimmed surface to test
    shared_ptr<ParamSurface> sf;
    RectDomain dom = surface->containingDomain();
    shared_ptr<BoundedSurface> bd_sf = 
      dynamic_pointer_cast<BoundedSurface, ParamSurface>(surface);
    if (bd_sf.get())
      {
	shared_ptr<ParamSurface> tmp = bd_sf->underlyingSurface();
	vector<shared_ptr<ParamSurface> > sub_sfs;
	try {
	  sub_sfs = tmp->subSurfaces(dom.umin(), dom.vmin(), 
				     dom.umax(), dom.vmax());
	}
	catch (...)
	  {
	    sf = tmp;
	  }

	if (sub_sfs.size() > 1)
	  sf = tmp;
	else if (sub_sfs.size() == 1)
	  sf = sub_sfs[0];
      }
    else
      sf = surface;
  
    // Fetch opposite boundary curves and check coincidence
    vector<shared_ptr<ParamCurve> > cvs_u1 = sf->constParamCurves(dom.umin(), false);
    vector<shared_ptr<ParamCurve> > cvs_u2 = sf->constParamCurves(dom.umax(), false);
    vector<shared_ptr<ParamCurve> > cvs_v1 = sf->constParamCurves(dom.vmin(), true);
    vector<shared_ptr<ParamCurve> > cvs_v2 = sf->constParamCurves(dom.vmax(), true);

    Identity ident;
    int coinc1 = ident.identicalCvs(cvs_u1[0], cvs_u1[0]->startparam(), cvs_u1[0]->endparam(),
					cvs_u2[0], cvs_u2[0]->startparam(), cvs_u2[0]->endparam(),
					tol);
    int coinc2 = ident.identicalCvs(cvs_v1[0], cvs_v1[0]->startparam(), cvs_v1[0]->endparam(),
					cvs_v2[0], cvs_v2[0]->startparam(), cvs_v2[0]->endparam(),
					tol);
    vector<shared_ptr<ParamSurface> > sub_sfs1;
    if (coinc1)
      {
	// Split in the first parameter direction
	vector<shared_ptr<ParamSurface> > sub_sfs2;
	double mid = 0.5*(dom.umin()+dom.umax());
	try {
	  sub_sfs1 = surface->subSurfaces(dom.umin(), dom.vmin(), 
					  mid, dom.vmax());
	}
	catch (...)
	  {
	    sfs.push_back(surface);
	    return sfs;
	  }
	
	try {
	  sub_sfs2 = surface->subSurfaces(mid, dom.vmin(), 
					  dom.umax(), dom.vmax());
	}
	catch(...)
	  {
	    sfs.push_back(surface);
	    return sfs;
	  }
	sub_sfs1.insert(sub_sfs1.end(), sub_sfs2.begin(), sub_sfs2.end());
      }
    else
      sub_sfs1.push_back(surface);
	
    if (coinc2)
      {
	for (size_t ki=0; ki<sub_sfs1.size(); ++ki)
	  {
	    // Split in the second parameter direction
	    RectDomain dom2 = sub_sfs1[ki]->containingDomain();
	    vector<shared_ptr<ParamSurface> > sub_sfs2;
	    double mid = 0.5*(dom2.vmin()+dom2.vmax());
	    try {
	      sub_sfs2 = sub_sfs1[ki]->subSurfaces(dom2.umin(), dom2.vmin(), 
						   dom2.umax(), mid);
	    }
	    catch (...)
	      {
		return sub_sfs1;
	      }
	    
	    sfs.insert(sfs.end(), sub_sfs2.begin(), sub_sfs2.end());
	    sub_sfs2.clear();
	    try {
	      sub_sfs2 = sub_sfs1[ki]->subSurfaces(dom2.umin(), mid,
						   dom2.umax(), dom2.vmax());
	    }
	    catch (...)
	      {
		return sub_sfs1;
	      }
	    sfs.insert(sfs.end(), sub_sfs2.begin(), sub_sfs2.end());
	  }
#ifdef DEBUG
	for (size_t kj=0; kj<sfs.size(); ++kj)
	  {
	    sfs[kj]->writeStandardHeader(of);
	    sfs[kj]->write(of);
	  }
#endif
	return sfs;
      }
    else
      {
#ifdef DEBUG
	for (size_t kj=0; kj<sub_sfs1.size(); ++kj)
	  {
	    sub_sfs1[kj]->writeStandardHeader(of);
	    sub_sfs1[kj]->write(of);
	  }
#endif
      return sub_sfs1;
      }
  }

}
