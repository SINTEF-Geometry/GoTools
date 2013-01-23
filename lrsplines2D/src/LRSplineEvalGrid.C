//===========================================================================
//                                                                           
// File: LRSplineEvalGrid.C                                                  
//                                                                           
// Created: Tue Jan 22 10:13:31 2013                                         
//                                                                           
// Author: Sverre Briseid <Sverre.Briseid@sintef.no>
//                                                                           
// Revision: $Id$
//                                                                           
// Description:
//                                                                           
//===========================================================================


#include "GoTools/lrsplines2D/LRSplineEvalGrid.h"



//==============================================================================
namespace Go
//==============================================================================
{


LRSplineEvalGrid::LRSplineEvalGrid()
  : dim_(0)
{
}


LRSplineEvalGrid::LRSplineEvalGrid(LRSplineSurface& lr_spline)
    : dim_(lr_spline.dimension())
{
    assert(!lr_spline.rational());

	orig_dom_ = lr_spline.parameterDomain();

    order_u_ = 1 + lr_spline.degree(XFIXED);
    order_v_ = 1 + lr_spline.degree(YFIXED);

    // We run through the sf and extract the elements.
    auto iter = lr_spline.elementsBegin();
    while (iter != lr_spline.elementsEnd())
    {
	elements_.push_back(*iter->second);
	++iter;
    }

    mesh_ = lr_spline.mesh(); 

}

} // end namespace Go
