//==========================================================================
//                                                                          
// File: ElementaryCurve.C                                                   
//                                                                          
// Created: Mon Sep  8 13:30:47 2008                                         
//                                                                          
// Author: Jan B. Thomassen <jbt@sintef.no>
//                                                                          
// Revision: $Id: ElementaryCurve.C,v 1.1 2008-09-09 11:24:23 jbt Exp $
//                                                                          
// Description:
//                                                                          
//==========================================================================


#include "GoTools/geometry/ElementaryCurve.h"
#include <limits>


using std::numeric_limits;


namespace Go {


// Constructor
//===========================================================================
ElementaryCurve::ElementaryCurve()
    : isReversed_(false)
//===========================================================================
{
}


// Destructor
//===========================================================================
ElementaryCurve::~ElementaryCurve()
//===========================================================================
{
}


//===========================================================================
void ElementaryCurve::reverseParameterDirection(bool switchparam)
//===========================================================================
{
    if (switchparam) {
        swapParameters2D();
        return;
    }
    
    isReversed_ = !isReversed_;
}


//===========================================================================
bool ElementaryCurve::isReversed() const
//===========================================================================
{
    return isReversed_;
}


//===========================================================================
void ElementaryCurve::getReversedParameter(double& t) const
//===========================================================================
{
    // startparam() and endparam() might be reversed, but their sum is
    // invariant.
    if (isReversed()) {
        double inf = numeric_limits<double>::infinity();
        bool isBounded = (startparam() > -inf && endparam() < inf);
        if (isBounded) {
            t = startparam() + endparam() - t;
        }
        else {
            t = -t;
        }
    }
}


//===========================================================================


} // namespace Go
