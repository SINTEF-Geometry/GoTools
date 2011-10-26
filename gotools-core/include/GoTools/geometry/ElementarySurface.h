//==========================================================================
//                                                                          
// File: ElementarySurface.h                                                 
//                                                                          
// Created: Wed Sep  3 17:40:46 2008                                         
//                                                                          
// Author: Jan B. Thomassen <jbt@sintef.no>
//                                                                          
// Revision: $Id: ElementarySurface.h,v 1.3 2008-12-10 14:52:05 vsk Exp $
//                                                                          
// Description:
//                                                                          
//==========================================================================

#ifndef _ELEMENTARYSURFACE_H
#define _ELEMENTARYSURFACE_H



#include "GoTools/geometry/ParamSurface.h"


namespace Go
{


class SplineSurface;


/// \brief ElementarySurface is a base class for elementary surfaces
/// like planes and cylinders. Such surfaces have natural
/// parametrizations and ElementartSurface is therefore a subclass of
/// ParamSurface. These surfaces are non-self-intersecting.
class ElementarySurface : public ParamSurface
{
public:
    /// Virtual destructor, enables safe inheritance.
    virtual ~ElementarySurface();

    // --- Functions inherited from GeomObject ---
    virtual ElementarySurface* clone() const = 0;

    // --- Functions inherited from ParamSurface ---

    RectDomain containingDomain() const;

    virtual Point closestInDomain(double u, double v) const;

    /// Check if a parameter pair lies inside the domain of this surface
    virtual bool inDomain(double u, double v) const;

    virtual SplineSurface* geometrySurface() const = 0;

    virtual double area(double tol) const;

    /// Return surface corners, geometric and parametric points
    /// in that sequence
    virtual void 
      getCornerPoints(std::vector<std::pair<Point,Point> >& corners) const;

    // --- Functions specific to ElemtarySurface ---
    virtual bool isBounded() const;
};


} // namespace Go



#endif // _ELEMENTARYSURFACE_H

