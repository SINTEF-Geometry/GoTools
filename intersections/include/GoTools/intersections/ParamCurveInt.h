//===========================================================================
//                                                                           
// File: ParamCurveInt.h 
//                                                                           
// Created: 
//                                                                           
// Author: vsk
//                                                                           
// Revision: $Id: ParamCurveInt.h,v 1.41 2007-01-15 10:12:38 vsk Exp $
//                                                                           
// Description:
//                                                                           
//===========================================================================

#ifndef _PARAMCURVEINT_H
#define _PARAMCURVEINT_H

#include "GoTools/intersections/ParamGeomInt.h"
#include "GoTools/geometry/ParamCurve.h"


namespace Go {


class RotatedBox;


/// This class represents the "intersection object" of a parametric
/// curve.

class ParamCurveInt : public ParamGeomInt {
public:
    /// Constructor.
    /// \param curve the parametric curve defining the intersection
    /// object.
    /// \param parent the parent object to this object. Can be either
    /// a curve or a surface.
    explicit ParamCurveInt(shared_ptr<ParamCurve> curve,
			   ParamGeomInt* parent = 0); 

    /// Destructor.
    virtual ~ParamCurveInt() {};

    /// Evaluate the object in the input parameter.
    /// \param res the Point to be returned.
    /// \param par the parameter in which to evaluate. The size of the
    /// array should be equal to numParams().
    virtual void point(Point& res, const double *par) const
    { curve_->point(res, par[0]); }

    /// Evaluate the object in the input parameter, with the specified
    /// number of derivatives.
    /// \param res the Point to be returned.
    /// \param par the parameter in which to evaluate. The size of
    /// the array should be equal to numParams().
    /// \param der the number of derivatives to calculate.
    /// \param from_right if true the evaluation is to be performed
    /// from the right side of the parameter value.
    /// \param resolution tolerance used when determining whether
    /// parameters are located at special values of the parameter
    /// domain (in particualar: knot values in case of spline objects.
    virtual void point(std::vector<Point>& res, 
		       const double* par, 
		       int der,
		       const bool* from_right = 0,
		       double resolution = 1.0e-12) const 
    {
	(from_right == 0)
	    ? curve_->point(res, *par, der, true)
	    : curve_->point(res, *par, der, from_right[0]);
    }

    /// Return a pointer to this object.
    /// \return Pointer to this object.
    virtual ParamCurveInt* getParamCurveInt();

    /// Return a NULL pointer (as this object is not of the correct
    /// type).
    /// \return NULL pointer (as this object is not of the correct
    /// type).
    virtual ParamSurfaceInt* getParamSurfaceInt()
    { return NULL; }

    /// Return pointer to the parametric curve defining this object.
    /// \return Pointer to the parametric curve defining this object.
    shared_ptr<ParamCurve> getParamCurve();

    /// Return pointer to the parametric curve defining this object.
    /// \return Pointer to the parametric curve defining this object.
    shared_ptr<const ParamCurve> getParamCurve() const;

    /// Return pointer to a subcurve of the parent curve for this
    /// object.  If no such curve exist, we use the curve of this
    /// object instead.
    /// \param start the start parameter of the subcurve.
    /// \param end the end parameter of the subcurve.
    /// \return Pointer to a subcurve of the parent curve for this
    /// object.
    shared_ptr<ParamCurve> 
    getParentParamCurve(double& start, double& end);

    /// Return pointer to the parent curve for this object.  If no
    /// such curve exist, we use the curve of this object instead.
    /// \return Pointer to a subcurve of the parent curve for this
    /// object.
    shared_ptr<ParamCurve> getParentParamCurve();

    /// Return an intersection object for the input curve.  This
    /// object is used as the parent for the intersection object.
    /// \param curve the parametric curve defining the intersection
    /// object.
    virtual shared_ptr<ParamCurveInt> 
    makeIntObject(shared_ptr<ParamCurve> curve);
    
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

    /// Set info about subdivision value
    virtual void setCriticalVal(int pardir, double par);

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
    /// \param eps the tolerance defining the boundary neighbourhood.
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
    virtual void
    subdivide(int pardir, double par, 
	      std::vector<shared_ptr<ParamGeomInt> >& subdiv_objs,
	      std::vector<shared_ptr<ParamGeomInt> >& bd_objs);

    /// Create the CompositeBox for the parametric object.
    /// \return The CompositeBox for the parametric object.
    virtual CompositeBox compositeBox() const;

    /// A cone which contains all normals of the object.
    /// \return A cone which contains all normals of the object.
    virtual DirectionCone directionCone() const;

    /// Return the boundary objects of this object.
    /// \param bd_objs the boundary objects of this object.
    virtual void 
    getBoundaryObjects(std::vector<shared_ptr<BoundaryGeomInt> >& bd_objs);

    /// Check if the object is periodic.  Analyze periodicity of curve
    /// based on number of repeating knots and control points. The
    /// return value is -1 if the curve ends are disjoint, otherwise k
    /// if cv is C^k continuous. These are sufficient but not
    /// necessary conditions for periodicity, so it is possible that a
    /// higher degree of periodicity exists.  Should not be called on
    /// this layer, should be overruled by inherited class.
    /// \param pardir the parameter direction in question. Does not
    /// pertain to for a curve.
    /// \return -1 if the curve ends are disjoint, or k if the curve
    /// is proven to be C^k continuous.
    virtual int checkPeriodicity(int pardir = 0) const;

    /// The dimension of the geometric space.
    virtual int dimension() const
    { return dim_; }

    /// The start parameter of the curve.
    /// \return The start parameter of the curve.
    double startparam() const 
    { return curve_->startparam(); }

    /// The end parameter of the curve.
    /// \return The end parameter of the curve.
    double endparam() const
    { return curve_->endparam(); }

    // Functions used from coincidence checking

    /// Make sure that the input parameter lies inside the range of
    /// the parametric curve.  Set t equal to tmin if it lies below
    /// tmin, or tmax if it lies above tmax.
    /// \param t the input parameter
    void assureInRange(double& t);
 
    /// Find the knot interval for which t lies inside, moving the
    /// value t if it lies close to a knot.
    /// \param t the parameter value
    /// \param tol the parametric tolerance deciding if the input
    /// parameter t should be moved.
    virtual int knotIntervalFuzzy(double& t, double tol) const; 
    
    /// Return the value of the knot next to the input parameter par.
    /// \param par the parameter value
    /// \param forward if true we return the closest knot to the
    /// right, otherwise the closest knot to the left.
    /// \param tol the tolerance to determine if \a par is already
    /// located on the start of the next segment.
    /// \return The knot closest to the input parameter.
    virtual double nextSegmentVal(double par,
				  bool forward, double tol) const;

    /// Verfify whether the object is degenerate in the specified
    /// direction and parameter.
    /// \param epsge the geometric tolerance defining degeneracy.
    /// \param dir the parameter direction in question.
    /// \param par the parameter in which to evaluate. Size of array
    /// should be equal to numParams().
    /// \return True if the object is degenerate.
    virtual bool isDegenerate(double epsge, int dir, double *par);

    /// Verify whether the object is a spline.
    /// \return True if the object is a spline.
    virtual bool isSpline();

    virtual const SplineCurve* getSpline()
	{ return 0; }

    /// We try to treat problems which will never result in a simple
    /// case by shrinking the domain slightly, resulting in smaller
    /// cones.  This is useful for scenarios where the normals are
    /// parallell in a boundary point.
    /// \param axis1 first vector defining a projection plane
    /// \param axis2 second vector defining a projection plane
    /// \return The optimized cone angle
    virtual double getOptimizedConeAngle(Point& axis1, Point& axis2);

    /// Create a box containing the geometric sample mesh in the input
    /// coordinate system.
    /// \param axis the axis defining the coordinate system. The size
    /// of vector axis is 2, the last point is given as the cross
    /// product axis[0]%axis[1].
    /// \return The rotated box
    virtual RotatedBox getRotatedBox(std::vector<Point>& axis) const;

    /// Create an axis by interpolating the end points.
    /// \return The created axis.
    void axisFromEndpts(Point& axis) const;

protected:
    // Data members
    shared_ptr<ParamCurve> curve_;

    int dim_; // Space dimension.

    std::vector<std::pair<double, int> > segment_;

    // Approximating polygon
    mutable std::vector<double> mesh_;

    mutable bool lw_set_;
    mutable double length_, wiggle_;

private:
    // Make approximating polygon
    void makeMesh(int size) const;

};


} // namespace Go


#endif // _PARAMCURVEINT_H
