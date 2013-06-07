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

#include "GoTools/implicitization/BernsteinUtils.h"
#include "GoTools/geometry/SplineCurve.h"
#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/utils/errormacros.h"


using namespace std;


namespace Go {


//==========================================================================
void spline_to_bernstein(const SplineCurve& seg, int dd,
			 BernsteinPoly& bp)
//==========================================================================
{
    int order = seg.order();
    ALWAYS_ERROR_IF(order != seg.numCoefs(),
		    "The curve is not Bezier");

    int dim = seg.dimension();
    bool rational = seg.rational();
    int effdim = rational ? dim+1 : dim;
    vector<double>::const_iterator it;
    if (rational)
	it = seg.rcoefs_begin();
    else
	it = seg.coefs_begin();
    vector<double> coefs(order);
    for (int i = 0; i < order; ++i)
	coefs[i] = it[effdim * i + dd];

    bp = BernsteinPoly(coefs);
    return;
}


//==========================================================================
void spline_to_bernstein(const SplineSurface& pat, int dd,
			 BernsteinMulti& bm)
//==========================================================================
{
    int order_u = pat.order_u();
    int order_v = pat.order_v();
    ALWAYS_ERROR_IF((order_u != pat.numCoefs_u())
		|| (order_v != pat.numCoefs_v()),
		"The surface is not Bezier");

    int dim = pat.dimension();
    bool rational = pat.rational();
    int effdim = rational ? dim+1 : dim;
    vector<double>::const_iterator it;
    if (rational)
	it = pat.rcoefs_begin();
    else
	it = pat.coefs_begin();
    vector<double> coefs(order_u * order_v);
    for (int iv = 0; iv < order_v; ++iv) {
	for (int iu = 0; iu < order_u; ++iu) {
	    int offset = order_u*iv + iu;
	    coefs[offset] = it[effdim * offset + dd];
	}
    }

    bm = BernsteinMulti(order_u - 1, order_v - 1, coefs);
    return;
}


//==========================================================================
void spline_to_bernstein(const SplineCurve& seg,
			 vector<BernsteinPoly>& seg_bp)
//==========================================================================
{
    int dim = seg.dimension();
    int effdim = seg.rational() ? dim+1 : dim;
    seg_bp.resize(effdim);

    for (int dd = 0; dd < effdim; ++dd) {
	spline_to_bernstein(seg, dd, seg_bp[dd]);
    }

    return;
}


//==========================================================================
void spline_to_bernstein(const SplineSurface& pat,
			 vector<BernsteinMulti>& pat_bm)
//==========================================================================
{
    int dim = pat.dimension();
    int effdim = pat.rational() ? dim+1 : dim;
    pat_bm.resize(effdim);

    for (int dd = 0; dd < effdim; ++dd) {
	spline_to_bernstein(pat, dd, pat_bm[dd]);
    }

    return;
}


} // namespace Go
