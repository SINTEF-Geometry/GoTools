//===========================================================================
//                                                                           
// File: Spline2FunctionInt.h                                                
//                                                                           
// Created: Mon Sep 27 16:12:40 2004                                         
//                                                                           
// Author: Sverre Briseid <Sverre.Briseid@sintef.no>
//                                                                           
// Revision: $Id: Spline2FunctionInt.h,v 1.17 2006-12-19 13:49:02 jbt Exp $
//                                                                           
// Description:
//                                                                           
//===========================================================================

#ifndef _SPLINE2FUNCTIONINT_H
#define _SPLINE2FUNCTIONINT_H


#include "GoTools/intersections/Param2FunctionInt.h"


namespace Go {


 /// This class represents the "intersection object" of a spline
 /// surface of dimension 1.

class Spline2FunctionInt : public Param2FunctionInt {
public:
    /// Constructor.
    /// \param surf the parametric 1-dimensional surface defining the
    /// object.
    explicit Spline2FunctionInt(std::shared_ptr<SplineSurface> surf);

    /// Constructor.
    /// \param surf the parametric 1-dimensional surface defining the
    /// object.
    /// \param parent the parent object to this object.
    explicit Spline2FunctionInt(std::shared_ptr<SplineSurface> surf, 
				Param2FunctionInt *parent);

    /// Destructor.
    virtual ~Spline2FunctionInt(){};

    /// Return an intersection object for the input surface.  This
    /// object is used as the parent for the intersection object.
    /// \param surf the parametric surface defining the intersection
    /// object.
    virtual std::shared_ptr<Param2FunctionInt>
    makeIntFunction(std::shared_ptr<ParamSurface> surf);
    
    /// Return true if the object has any inner knots in the specified
    /// parameter direction.
    /// \param pardir the parameter direction in question. Indexing
    /// starts at 0.
    virtual bool hasInnerKnots(int pardir) const;

    /// Return true if the object has any critical parameter values or
    /// inner knots in the specified parameter direction.
    /// \param pardir the parameter direction in question. Indexing
    /// starts at 0.
    virtual bool hasCriticalValsOrKnots(int pardir) const;

    /// Return the inner knot values in the specified direction.
    /// \param pardir the parameter direction in question. Indexing
    /// starts at 0.
    /// \param sort the returned values may be sorted by the function.
    virtual std::vector<double> getInnerKnotVals(int pardir,
						 bool sort=false) const;

    /// Return the critical parameter values and inner knots for
    /// object.
    /// \param pardir the parameter direction in question. Indexing
    /// starts at 0.
    virtual std::vector<double> getCriticalValsAndKnots(int pardir) const;

//     void getLengthAndWiggle(double *length, double *wiggle);

    // pardir = 0 || 1
    // @@sbr A trimmed sf need not have a startparam/endparam ...
    /// Return the start parameter value in the specified direction.
    /// \param pardir the parameter direction in question. Indexing
    /// starts at 0.
    /// \return The start parameter in the specified direction.
    virtual double startParam(int pardir) const;

    /// Return the end parameter in the specified direction.
    /// \param pardir the parameter direction in question. Indexing
    /// starts at 0.
    /// \return The end parameter in the specified direction.
    virtual double endParam(int pardir) const;

    // We want to know if there exists par direction in which the
    // function is monotone (guarantees at most one solution for a
    // continuous function).  We also would like to know in what
    // direction the function is monotone (will then move in
    // orthogonal direction).
    /// Return true if the surface is monotone in any direction.
    /// \param dir the direction in which the surface is
    /// monotone. Pertains only if surface is monotone.
    /// \return Whether or not the surface is monotone.
    virtual bool monotone(Point& dir, double tol=1.0e-15) const;

    // Functions used for coincidence testing

    // Find the interval in which t lies, if t is within tol of an
    // existing knot, t will be changed to that knot.
    /// Return the knot interval for which t lies inside, moving the
    /// value t if it lies close to a knot.
    /// \param pardir the parameter direction in question. Indexing
    /// starts at 0.
    /// \param t the parameter value
    /// \param tol the parametric tolerance deciding if the input
    /// parameter t should be moved.
    virtual int knotIntervalFuzzy(int pardir, double& t, double tol) const;

    /// Return the value of the knot next to the input parameter par.
    /// \param pardir the parameter direction in question. Indexing
    /// starts at 0.
    /// \param par the parameter value
    /// \param forward if true we return the closest knot to the
    /// right, otherwise the closest knot to the left.
    /// \return The knot closest to the input parameter.
    virtual double nextSegmentVal(int pardir, double par, bool forward) const;

    /// Return the boundary objects of this object.
    /// \param bd_objs the boundary objects of this object.
    virtual void 
    getBoundaryObjects(std::vector<std::
		       shared_ptr<BoundaryFunctionInt> >& bd_objs);

    /// Return a 3-dimensional visualization spline surface.
    /// \return The pointer to the 3-dimensional spline surface (u, v,
    /// surf_(u, v)).
    std::shared_ptr<SplineSurface> surface3D();

    /// Return the gradient of the function S as a two-dimensional
    /// spline surface: (dS/du, dS/dv).
    /// \return Pointer to the gradient spline surface
    std::shared_ptr<SplineSurface> createGradSurface() const;

protected:
    // Data members
    std::shared_ptr<SplineSurface> spsf_; // shared_ptr to this
					    // curve

};


} // namespace Go


#endif // _SPLINE2FUNCTIONINT_H

