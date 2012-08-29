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
class ElementaryCurve;


/// \brief ElementarySurface is a base class for elementary surfaces
/// like planes and cylinders. Such surfaces have natural
/// parametrizations and ElementartSurface is therefore a subclass of
/// ParamSurface. These surfaces are non-self-intersecting.
class ElementarySurface : public ParamSurface
{
public:
    ElementarySurface();

    /// Virtual destructor, enables safe inheritance.
    virtual ~ElementarySurface();

    // --- Functions inherited from GeomObject ---
    virtual ElementarySurface* clone() const = 0;

    // --- Functions inherited from ParamSurface ---

    virtual CurveLoop outerBoundaryLoop(double degenerate_epsilon
					  = DEFAULT_SPACE_EPSILON) const;

    RectDomain containingDomain() const;

    virtual Point closestInDomain(double u, double v) const;

    /// Check if a parameter pair lies inside the domain of this surface
    virtual bool inDomain(double u, double v) const;

    virtual double area(double tol) const;

    /// Return surface corners, geometric and parametric points
    /// in that sequence
    virtual void 
      getCornerPoints(std::vector<std::pair<Point,Point> >& corners) const;

    virtual SplineSurface* asSplineSurface();

    // --- Functions specific to ElementarySurface ---
    virtual bool isBounded() const;

    virtual bool isClosed(bool& closed_dir_u, bool& closed_dir_v) const;

    virtual SplineSurface* geometrySurface() const = 0;
    /// Create a SplineSurface representation of the elementary surface
    virtual SplineSurface* createSplineSurface() const = 0;

    /// Limit the surface by limiting the parameter domain
    virtual void setParameterBounds(double from_upar, double from_vpar,
				    double to_upar, double to_vpar) = 0;

    /// Fetch the parameter curve in the domain of the elementary surface
    /// corresponding to a given elementary curve in geometry space
    /// if this curve has a simpler elementary representation.
    /// Otherwise, nothing is returned
    virtual shared_ptr<ElementaryCurve> 
      getElementaryParamCurve(ElementaryCurve* space_crv, double tol) const;

    virtual void turnOrientation();
    virtual void reverseParameterDirection(bool direction_is_u);
    virtual void swapParameterDirection();

    //virtual bool isReversedU() const;
    //virtual bool isReversedV() const;
    virtual bool isSwapped() const;


    
protected:
    //bool isReversedU_;
    //bool isReversedV_;
    bool isSwapped_;

    // Helper function to be used in functions like Cylinder::point(), where
    // we need to take the isSwapped_ flag into account.
    void getOrientedParameters(double& u, double&v) const;

};


} // namespace Go



#endif // _ELEMENTARYSURFACE_H

