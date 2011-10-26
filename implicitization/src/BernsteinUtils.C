//==========================================================================
//                                                                          
// File: BernsteinUtils.C                                                    
//                                                                          
// Created: Thu Jan 15 16:43:33 2004                                         
//                                                                          
// Author: Jan B. Thomassen <jbt@sintef.no>
//                                                                          
// Revision: $Id: BernsteinUtils.C,v 1.6 2006-01-27 12:53:23 jbt Exp $
//                                                                          
// Description:
//                                                                          
//==========================================================================


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
