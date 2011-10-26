//===========================================================================
//                                                                           
// File: Spline1FunctionInt.h                                                
//                                                                           
// Created: Fri Sep 24 15:51:57 2004                                         
//                                                                           
// Author: Sverre Briseid <Sverre.Briseid@sintef.no>
//                                                                           
// Revision: $Id: Spline1FunctionInt.h,v 1.15 2006-09-01 12:25:26 vsk Exp $
//                                                                           
// Description:
//                                                                           
//===========================================================================

#ifndef _SPLINE1FUNCTIONINT_H
#define _SPLINE1FUNCTIONINT_H


#include "GoTools/intersections/Param1FunctionInt.h"


namespace Go {


class SplineCurve;


/// Class that represents the "intersection object" of a spline curve
/// of dimension 1.

class Spline1FunctionInt : public Param1FunctionInt {
public:
    /// Constructor.
    /// Input curve must be of type SplineCurve. This is not checked
    /// run-time, so we rely on the user to obey this rule.
    /// \param curve the parametric 1-dimensional curve defining the
    /// object.
    explicit Spline1FunctionInt(std::shared_ptr<ParamCurve> curve);

    /// Constructor.
    /// Input curve must be of type SplineCurve. This is not checked
    /// run-time, so we rely on the user to obey this rule.
    /// \param curve the parametric 1-dimensional curve defining the
    /// object. Can be either a curve or a surface.
    /// \param parent the parent object to this object.
    explicit Spline1FunctionInt(std::shared_ptr<ParamCurve> curve, 
				ParamFunctionInt *parent);

    /// Destructor.
    virtual ~Spline1FunctionInt() {};

    /// Return an intersection object for the input curve.
    /// Input curve must be of type SplineCurve. This is not checked
    /// run-time, so we rely on the user to obey this rule.  This
    /// object is used as the parent for the intersection object.
    /// \param curve the parametric curve defining the intersection
    /// object.
    virtual std::shared_ptr<Param1FunctionInt> 
    makeIntFunction(std::shared_ptr<ParamCurve> curve);
    
    /// Return true if the object has any inner knots in the specified
    /// parameter direction.
    /// \param pardir the parameter direction in question. Indexing
    /// starts at 0.
    virtual bool hasInnerKnots(int pardir) const;

    /// Return true if the object has any critical parameter values or
    /// inner knots.
    /// \param pardir the parameter direction in question. Indexing
    /// starts at 0.
    virtual bool hasCriticalValsOrKnots(int pardir) const;

    /// Return the inner knot values in the specified direction.
    /// \param pardir the parameter direction in question. Indexing
    /// starts at 0.
    /// \param sort the returned values may be sorted by the function.
    virtual std::vector<double> getInnerKnotVals(int pardir,
						 bool sort = false) const;

    /// Return the critical parameter values and inner knots for
    /// object.
    /// \param pardir the parameter direction in question. Indexing
    /// starts at 0.
    virtual std::vector<double> getCriticalValsAndKnots(int pardir) const; 

    /// Return the size of the geometric sample mesh in the specified
    /// direction.
    /// \param dir the parameter direction in question. Indexing
    /// starts at 0.
    virtual int getMeshSize(int dir);

    /// Return the corresponding mesh parameter.
    /// \param dir the parameter direction in question. Indexing
    /// starts at 0.
    /// \param idx the mesh idx in the specified direction. Indexing
    /// starts at 0.
    virtual double paramFromMesh(int dir, int idx);

    /// Return the geometric sample mesh for the spline function.
    /// \return The geometric sample mesh for the spline function.
    virtual std::vector<double>::iterator getMesh();

    // Accepts non-strict monotonicity, as long as sf is not totally flat.
    /// Return true if the curve is monotone.
    /// \param dir the direction in which the object is monotone. Is
    /// not of interest here as the curve has only 1 parameter
    /// direction.
    /// \return Whether or not the curve is monotone.
    virtual bool monotone(Point& dir, double tol=1.0e-15) const;

    // Functions used for coincidence testint

    // Find the interval in which t lies, if t is within tol of an
    // existing knot, t will be changed to that knot.
    /// Return the knot interval for which t lies inside, moving the
    /// value t if it lies close to a knot.
    /// \param t the parameter value
    /// \param tol the parametric tolerance deciding if the input
    /// parameter t should be moved.
    virtual int knotIntervalFuzzy(double& t, double tol) const;

    /// Return the value of the knot next to the input parameter par.
    /// \param par the parameter value
    /// \param forward if true we return the closest knot to the
    /// right, otherwise the closest knot to the left.
    /// \return The knot closest to the input parameter.
    virtual double nextSegmentVal(double par, bool forward) const;

protected:

    // Data members
    std::shared_ptr<SplineCurve> spcv_; // Same cv as in
					  // Param1FunctionInt, casted
					  // to spline.

};


} // namespace Go


#endif // _SPLINE1FUNCTIONINT_H

    
