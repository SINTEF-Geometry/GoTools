//==========================================================================
//                                                                          
// File: GGUsplit.C                                                          
//                                                                          
// Created: Wed Jan 22 17:30:30 2003                                         
//                                                                          
// Author: Jan B. Thomassen <jbt@sintef.no>
//                                                                          
// Revision: $Id: GGUsplit.C,v 1.10 2006-05-03 10:50:20 jbt Exp $
//                                                                          
// Description:
//                                                                          
//==========================================================================


#include "GoTools/geometry/GeometryTools.h"
#include <functional>


using std::vector;
using std::max;
using std::min;


namespace Go {


//==========================================================================
void splitCurveIntoSegments(const SplineCurve& cv,
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
void splitSurfaceIntoPatches(const SplineSurface& sf,
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
