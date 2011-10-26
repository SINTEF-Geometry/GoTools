//===========================================================================
//                                                                           
// File: ParamObjectInt.h
//                                                                           
// Created: 
//                                                                           
// Author: vsk
//                                                                           
// Revision: $Id: ParamObjectInt.h,v 1.30 2007-06-22 16:33:21 vsk Exp $
//                                                                           
// Description:
//                                                                           
//===========================================================================

#ifndef _PARAMOBJECTINT_H
#define _PARAMOBJECTINT_H


#include "GoTools/intersections/GeomObjectInt.h"
#include "GoTools/utils/Point.h"
#include <vector>


namespace Go {


class CompositeBox;
class ParamPointInt;
class ParamCurveInt;
class ParamSurfaceInt;


/// This class is a base class providing an interface to the
/// "intersection objects". The object is parametric.

class ParamObjectInt : public GeomObjectInt {
public:
    /// Constructor.
    /// \param parent is a parametric object for which this object
    /// constitues only a subpart (0 if there is no parent).
    ParamObjectInt(ParamObjectInt* parent = 0) : parent_(parent) {}
    
    /// Destructor.
    virtual ~ParamObjectInt() {}
    
    /// Return pointer to a parametric intersection point.
    /// \return Pointer to a parametric intersection point. If that is
    /// not the class type for this object a NULL pointer is returned.
    virtual ParamPointInt* getParamPointInt()
    { return 0; }

     /// Return pointer to a parametric intersection curve.
    /// \return Pointer to a parametric intersection curve. If that is
    /// not the class type for this object a NULL pointer is returned.
    virtual ParamCurveInt* getParamCurveInt()
    { return 0; }

    /// Return pointer to a parametric intersection surface.
    /// \return Pointer to a parametric intersection surface. If that
    /// is not the class type for this object a NULL pointer is
    /// returned.
    virtual ParamSurfaceInt* getParamSurfaceInt()
    { return 0; }

   /// Evaluate the object in the input parameter.
    /// \param pt the Point to be returned.
    /// \param tpar the parameter in which to evaluate. The size of
    /// the array should be equal to numParams().
    virtual void point(Point& pt, const double* tpar) const = 0;

    /// Evaluate the object in the input parameter, with the specified
    /// number of derivatives.
    /// \param pt the Point to be returned.
    /// \param tpar the parameter in which to evaluate. The size of
    /// the array should be equal to numParams().
    /// \param derivs the number of derivatives to calculate.
    /// \param from_right if true the evaluation is to be performed
    /// from the right of the parameter value.
    /// \param resolution tolerance used when determining whether
    /// parameters are located at special values of the parameter
    /// domain (in particualar: knot values in case of spline
    /// objects).
    virtual void point(std::vector<Point>& pt, 
		       const double* tpar, 
		       int derivs,  
		       const bool* from_right = 0,
		       double resolution = 1.0e-12) const = 0;

    /// The number of parameters in the object.
    virtual int numParams() const = 0;

    /// Return an estimate on the size and wiggle of the object.
    /// \param length the approximative length of the object in the
    /// corresponding parameter directions.  The size of the array
    /// should be equal to numParams().
    /// \param wiggle a scalar representing the wiggle of the object
    /// in the corresponding parameter directions. The size of the
    /// array should be equal to numParams().
    virtual void getLengthAndWiggle(double *length, double *wiggle) = 0;

    /// Return true if the object has any inner knots in the specified
    /// parameter direction.
    /// \param pardir the parameter direction in question. Indexing
    /// starts at 0.
    virtual bool hasInnerKnots(int pardir) const = 0;

    /// Return true if the object has any critical parameter values in
    /// the specified parameter direction.
    /// \param pardir the parameter direction in question. Indexing
    /// starts at 0.
    virtual bool hasCriticalVals(int pardir) const = 0;

    /// Set info about subdivision value
    virtual void setCriticalVal(int pardir, double par) 
	{;}  // Currently no action

    /// Return true if the object has any critical parameter values or
    /// inner knots in the specified parameter direction.
    /// \param pardir the parameter direction in question. Indexing
    /// starts at 0.
    virtual bool hasCriticalValsOrKnots(int pardir) const = 0;

    /// Return true if we are allowed to divide in the specified
    /// parameter direction.
    /// \param pardir the parameter direction in question.
    virtual bool canDivide(int pardir) = 0;

    virtual bool canDivideTinyTriang(int pardir)
	{ return true;  // Mostly not an issule 
	}

    /// Return the critical parameter values in the specified
    /// direction.
    /// \param pardir the parameter direction in question. Indexing
    /// starts at 0.
    virtual std::vector<double> getCriticalVals(int pardir) const = 0;

    /// Return the inner knot values in the specified direction.
    /// \param pardir the parameter direction in question. Indexing
    /// starts at 0.
    /// \param sort the returned values may be sorted by the function.
    virtual std::vector<double> getInnerKnotVals(int pardir,
						 bool sort=false) const = 0;

    /// Return the critical parameter values and inner knots for
    /// object.
    /// \param pardir the parameter direction in question. Indexing
    /// starts at 0.
    virtual std::vector<double> getCriticalValsAndKnots(int pardir) const = 0;

    /// Return the start parameter value in the specified direction.
    /// \param pardir the parameter direction in question. Indexing
    /// starts at 0.
    virtual double startParam(int pardir) const = 0;

    /// Return the end parameter in the specified direction.
    /// \param pardir the parameter direction in question. Indexing
    /// starts at 0.
    virtual double endParam(int pardir) const = 0;

    /// Return true if the object is degenerate in the specified
    /// direction and parameter.
    /// \param epsge the geometric tolerance defining degeneracy.
    /// \param dir the parameter direction in question.
    /// \param par the parameter in which to evaluate. Size of array
    /// should be equal to numParams().
    virtual bool isDegenerate(double epsge, int dir, double *par)
    { return false; }

    /// Return true if the specified point lies within eps from the
    /// boundary.
    /// \param par the parameter in which to evaluate. Size of array
    /// should be equal to numParams().
    /// \param eps the tolerance defining boundary neighbourhood.
    virtual bool boundaryPoint(const double* par, double eps) const = 0;

    /// Set the parent of this object.
    /// \param parent the parent of this object.
    void setParent(ParamObjectInt* parent)
    { parent_ = parent; }

    /// Return the parent of this object.
    ParamObjectInt* getParent() const 
    { return parent_; }

    /// Return the CompositeBox for the parametric object.
    virtual CompositeBox compositeBox() const = 0;

    /// Return an ancestor of the same type.  If the object has no
    /// parent, or if its parent has a higher number of parameters,
    /// return a pointer to itself.  Else, call this function
    /// recursively on the parent.
    const ParamObjectInt* getSameTypeAncestor() const
    {
	if (parent_ && (parent_->numParams() == numParams())) {
	    return parent_->getSameTypeAncestor();
	} else {
	    return this;
	}
    }

    /// Check if parameter point lies close to a corner of the
    /// parameter domain.
    /// \param par the input parameter point.
    /// \param epspar the parametric tolerance defining the
    /// neighbourhood.
    virtual bool inCorner(const double *par, double epspar) const
	{
	    return false;  // Relevant for surfaces and overriden for those
	}

protected:

    ParamObjectInt* parent_;

};


} // namespace Go


#endif // _PARAMOBJECTINT_H
