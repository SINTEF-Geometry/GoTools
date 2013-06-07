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

#include "GoTools/geometry/GeometryTools.h"
#include <functional>


using std::vector;
using std::max;
using std::min;


namespace Go {


//==========================================================================
void GeometryTools::splitCurveIntoSegments(const SplineCurve& cv,
			    vector<SplineCurve>& seg)
//==========================================================================
{
    SplineCurve orig = cv;
    orig.makeBernsteinKnots();

    int n = orig.numCoefs();
    int order = orig.order();
    int numseg = n / order;

    seg.resize(numseg);
    vector<double>::const_iterator it = orig.basis().begin();
    for (int i = 0; i < numseg; ++i) {
	shared_ptr<ParamCurve> new_cv(orig.subCurve(*it, *(it+order)));
	seg[i] = dynamic_cast<SplineCurve&>(*new_cv);
	it += order;
    }

    return;
}


//==========================================================================
void GeometryTools::splitSurfaceIntoPatches(const SplineSurface& sf,
			     vector<SplineSurface>& pat)
//==========================================================================
{
    SplineSurface orig = sf;
    orig.makeBernsteinKnotsU();
    orig.makeBernsteinKnotsV();

    int num_u = orig.numCoefs_u();
    int num_v = orig.numCoefs_v();
    int order_u = orig.order_u();
    int order_v = orig.order_v();
    int numpat_u = num_u / order_u;
    int numpat_v = num_v / order_v;

    pat.resize(numpat_u * numpat_v);
    typedef vector<double>::const_iterator const_iter;
    const_iter itu = orig.basis_u().begin();
    const_iter itv;
    for (int i = 0; i < numpat_u; ++i) {
	itv = orig.basis_v().begin();
	for (int j = 0; j < numpat_v; ++j) {
	    shared_ptr<SplineSurface>
		new_sf(orig.subSurface(*itu, *itv,
				       *(itu+order_u), *(itv+order_v)));
	    pat[numpat_u*j + i] = *new_sf;
	    itv += order_v;
	}
	itu += order_u;
    }

    return;
}


//==========================================================================


} // namespace Go
