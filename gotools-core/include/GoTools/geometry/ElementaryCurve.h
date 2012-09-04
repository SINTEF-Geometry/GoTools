//==========================================================================
//                                                                          
// File: ElementaryCurve.h                                                   
//                                                                          
// Created: Mon Sep  8 11:16:43 2008                                         
//                                                                          
// Author: Jan B. Thomassen <jbt@sintef.no>
//                                                                          
// Revision: $Id: ElementaryCurve.h,v 1.1 2008-09-09 11:24:23 jbt Exp $
//                                                                          
// Description:
//                                                                          
//==========================================================================

#ifndef _ELEMENTARYCURVE_H
#define _ELEMENTARYCURVE_H


#include "GoTools/geometry/ParamCurve.h"
#include "GoTools/geometry/SplineCurve.h"


namespace Go {


/// \brief ElementaryCurve is a base class for elementary curves like
/// lines and circles. Such curves have natural parametrizations and
/// ElementaryCurve is therefore a subclass of ParamCurve.
class ElementaryCurve : public ParamCurve
{
public:
    ElementaryCurve();

    /// Virtual destructor, enables safe inheritance.
    virtual ~ElementaryCurve();

    // --- Functions inherited from GeomObject ---
    virtual ElementaryCurve* clone() const = 0;

    // --- Functions inherited from ParamCurve ---
    virtual void reverseParameterDirection(bool switchparam = false);

    // --- Functions specific to ElementaryCurve ---
    virtual SplineCurve* createSplineCurve() const = 0;

    /// Set bounds for the parametrization of the curve
    /// \param startpar start parameter
    /// \param endpar end parameter
    virtual void setParamBounds(double startpar, double endpar) = 0;

    // Translate the curve along a given vector
    virtual void translateCurve(const Point& dir) = 0;

    virtual bool isReversed() const;

    // Helper function for reverseParameterDirection() when switchparam
    // is true. Used when curve is a parameter curve and x and y coordinates
    // should be swapped.
    virtual void swapParameters2D() = 0;

protected:
    // Returns reversed parameter in [tmin, tmax] if isReversed_ is true
    void getReversedParameter(double& t) const;

    bool isReversed_;
};


} // namespace Go


#endif // _ELEMENTARYCURVE_H

