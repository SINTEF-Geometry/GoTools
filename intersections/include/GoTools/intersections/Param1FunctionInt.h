//===========================================================================
//                                                                           
// File: Par1FuncInt.h                                                       
//                                                                           
// Created: Tue Sep 21 11:45:24 2004                                         
//                                                                           
// Author: Sverre Briseid <Sverre.Briseid@sintef.no>
//                                                                           
// Revision: $Id: Param1FunctionInt.h,v 1.24 2006-11-03 14:15:05 jbt Exp $
//                                                                           
// Description: Used to compute zeros of a one-parameter function
//                                                                           
//===========================================================================

#ifndef _PAR1FUNCINT_H
#define _PAR1FUNCINT_H


#include "GoTools/intersections/ParamFunctionInt.h"
#include "GoTools/geometry/ParamCurve.h"


namespace Go {


/// Class that represents an "intersection object" of a parametric
/// curve of dimension 1.

class Param1FunctionInt : public ParamFunctionInt {
public:
    /// Constructor.
    /// \param curve the parametric 1-dimensional curve defining the
    /// intersection object.
    explicit Param1FunctionInt(shared_ptr<ParamCurve> curve);

    /// Constructor.
    /// \param curve the parametric 1-dimensional curve defining the
    /// intersection object.
    /// \param parent the parent object to this object. Can be either
    /// a curve or a surface.
    explicit Param1FunctionInt(shared_ptr<ParamCurve> curve, 
			       ParamFunctionInt *parent);

    /// Destructor.
    virtual ~Param1FunctionInt();

    /// Evaluate the object in the input parameter.
    /// \param res the Point to be returned.
    /// \param par the parameter in which to evaluate. The size of the
    /// array should be equal to numParams().
    virtual void point(Point& res, const double *par) const
    { curve_->point(res, par[0]); }

    /// Evaluate the object in the input parameter, with the specified
    /// number of derivatives.
    /// \param res the Point to be returned.
    /// \param par the parameter in which to evaluate. The size of the
    /// array should be equal to numParams().
    /// \param der the number of derivatives to calculate.
    virtual void point(std::vector<Point>& res, double par, int der)
    { curve_->point(res, par, der); }

    /// Evaluate the object in the input parameter, with the specified
    /// number of derivatives.
    /// \param pt the Point to be returned.
    /// \param tpar the parameter in which to evaluate. The size of
    /// the array should be equal to numParams().
    /// \param derivs the number of derivatives to calculate.
    /// \param from_right if true the evaluation is to be performed
    /// from the right side of the parameter value.
    /// \param resolution tolerance used when determining whether
    /// parameters are located at special values of the parameter
    /// domain (in particualar: knot values in case of spline objects.
    virtual void point(std::vector<Point>& pt, 
		       const double* tpar, 
		       int derivs,  
		       const bool* from_right = 0,
		       double resolution = 1.0e-12) const
    {
	bool fr = (from_right != 0) ? from_right[0] : true;
	pt = curve_->point(*tpar, derivs, fr);
    }

    /// Return pointer to this object.
    virtual Param1FunctionInt* getParam1FunctionInt();

    /// Return pointer to the parametric curve.
    shared_ptr<ParamCurve> getParamCurve();
    /// Return pointer to the parametric curve.
    shared_ptr<const ParamCurve> getParamCurve() const;

    /// Return pointer to a subpart of the parent curve of this
    /// object.  If a parent curve does not exist, return pointer to
    /// curve in this object.  To reduce numerical noise we go
    /// straight to the source (undivided) curve.
    shared_ptr<ParamCurve> getParentParamCurve(double& start,
						      double& end);
    /// Return pointer to the parent curve of this object.  If a
    /// parent curve does not exist, return pointer to curve in this
    /// object.  To reduce numerical noise we go straight to the
    /// source (undivided) curve.
    shared_ptr<ParamCurve> getParentParamCurve();

    /// Return an intersection object for the input curve.  This
    /// object is used as the parent for the intersection object.
    /// \param curve the parametric curve defining the intersection
    /// object.
    virtual shared_ptr<Param1FunctionInt> 
    makeIntFunction(shared_ptr<ParamCurve> curve);
    
    /// The number of parameters in the object.
    virtual int numParams() const;

    /// Return an estimate on the size and wiggle of the object.
    /// \param length the approximative length of the object in the
    /// corresponding parameter directions.  The size of the array
    /// should be equal to numParams().
    /// \param wiggle a scalar representing the wiggle of the object
    /// in the corresponding parameter directions. The size of the
    /// array should be equal to numParams().
    virtual void getLengthAndWiggle(double *length, double *wiggle);

    /// Return true if the object has any inner knots in the specified
    /// parameter direction.
    /// \param pardir the parameter direction in question. Indexing
    /// starts at 0.
    virtual bool hasInnerKnots(int pardir) const;

    /// Return true if the object has any critical parameter values in
    /// the specified parameter direction.
    /// \param pardir the parameter direction in question. Indexing
    /// starts at 0.
    virtual bool hasCriticalVals(int pardir) const;

    /// Return true if the object has any critical parameter values or
    /// inner knots in the specified parameter direction.
    /// \param pardir the parameter direction in question. Indexing
    /// starts at 0.
    virtual bool hasCriticalValsOrKnots(int pardir) const;

    /// Return true if we are allowed to divide in the specified
    /// parameter direction.
    /// \param pardir the parameter direction in question.
    virtual bool canDivide(int pardir);

    /// Return the critical parameter values in the specified
    /// direction.
    /// \param pardir the parameter direction in question. Indexing
    /// starts at 0.
    virtual std::vector<double> getCriticalVals(int pardir) const;

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

    /// Return the geometric sample mesh for the parametric function.
    virtual std::vector<double>::iterator getMesh();

    /// Return the start parameter value in the specified direction.
    /// \param pardir the parameter direction in question. Indexing
    /// starts at 0.
    virtual double startParam(int pardir) const
    { return curve_->startparam(); }

    /// Return the end parameter in the specified direction.
    /// \param pardir the parameter direction in question. Indexing
    /// starts at 0.
    virtual double endParam(int pardir) const
    { return curve_->endparam(); }

    /// Return true if the specified point lies within eps from the
    /// boundary.
    /// \param par the parameter in which to evaluate. Size of array
    /// should be equal to numParams().
    /// \param eps the tolerance defining boundary neighbourhood.
    virtual bool boundaryPoint(const double* par, double eps) const;

    /// Subdivide the object in the specified parameter direction and
    /// parameter value.
    /// \param pardir direction in which to subdive. Indexing starts
    /// at 0.
    /// \param par parameter in which to subdivide.
    /// \param subdiv_objs The subparts of this object. Of the same
    /// geometric dimension as this object.
    /// \param bd_objs the boundaries between the returned \a
    /// subdiv_objs. Of geometric dimension 1 less than this object.
    virtual void subdivide(int pardir, double par, 
			   std::vector<shared_ptr<ParamFunctionInt> >& subdiv_objs,
			   std::vector<shared_ptr<ParamFunctionInt> >& bd_objs);

    // Accepts non-strict monotonicity, as long as sf is not totally
    // flat.
    /// Return true if the curve is monotone.
    /// \param dir the direction in which the object is monotone. Is
    /// not of interest here as the curve has only 1 parameter
    /// direction.
    /// \return Whether or not the curve is monotone.
    virtual bool monotone(Point& dir, double tol=1.0e-15) const; // = 0;

    /// Return the CompositeBox for the parametric object.
    virtual CompositeBox compositeBox() const;

    /// Return the boundary objects of this object.
    /// \param bd_objs the boundary objects of this object.
    virtual void 
    getBoundaryObjects(std::vector<shared_ptr<BoundaryFunctionInt> >& bd_objs);

    /// Return the dimension of the geometric space.
    int dimension()
    { return dim_; }

    /// Return the start parameter of the curve.
    double startparam() const
    { return curve_->startparam(); }
    /// Return the end parameter of the curve.
    double endparam() const
    { return curve_->endparam(); }

    // Functions used from coincidence checking

    /// Make sure that the input parameter lies inside the range of
    /// the parametric curve.  Set t equal to tmin if it lies below
    /// tmin, or tmax if it lies above tmax.
    /// \param t the input parameter
    void assureInRange(double& t);
 
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

    /// Return true if the object is degenerate in the specified
    /// direction and parameter.
    /// \param epsge the geometric tolerance defining degeneracy.
    /// \param dir the parameter direction in question.
    /// \param par the parameter in which to evaluate. Size of array
    /// should be equal to numParams().
    virtual bool isDegenerate(double epsge, int dir, double *par);

protected:

    // Data members
    // The curve defining the parameter domain.
    // However it need not be 1d, could be part of a composite
    // function.
    // Other parts of the curve (as well as evaluators) will then lie
    // in a subclass.
    shared_ptr<ParamCurve> curve_; // Our param curve (R->R).

    int dim_; // Space dimension. 1.

    Param1FunctionInt *parentcurve_; // @@sbr Or surface ... Currently cv.

    std::vector<std::pair<double, int> > segment_;

    // Approximating polygon
    mutable std::vector<double> mesh_;

    double deg_tol_; // Object should not be subdivided if it is
		     // already degenerate.

private:
    // Make approximating polygon
    void makeMesh(int size);

};


} // namespace Go


#endif // _PAR1FUNCINT_H

