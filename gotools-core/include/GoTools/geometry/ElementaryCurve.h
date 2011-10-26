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


namespace Go {


/// \brief ElementaryCurve is a base class for elementary curves like
/// lines and circles. Such curves have natural parametrizations and
/// ElementaryCurve is therefore a subclass of ParamCurve.
class ElementaryCurve : public ParamCurve
{
public:
    /// Virtual destructor, enables safe inheritance.
    virtual ~ElementaryCurve();

    // --- Functions inherited from GeomObject ---
    virtual ElementaryCurve* clone() const = 0;

    // --- Functions inherited from ParamCurve ---

};


} // namespace Go


#endif // _ELEMENTARYCURVE_H

