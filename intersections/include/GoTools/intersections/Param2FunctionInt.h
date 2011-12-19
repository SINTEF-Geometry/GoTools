//===========================================================================
//                                                                           
// File: Param2FunctionInt.h                                                 
//                                                                           
// Created: Mon Sep 27 14:26:41 2004                                         
//                                                                           
// Author: Sverre Briseid <Sverre.Briseid@sintef.no>
//                                                                           
// Revision: $Id: Param2FunctionInt.h,v 1.25 2006-11-03 14:15:05 jbt Exp $
//                                                                           
// Description:
//                                                                           
//===========================================================================

#ifndef _PARAM2FUNCTIONINT_H
#define _PARAM2FUNCTIONINT_H


#include "GoTools/intersections/ParamFunctionInt.h"
#include "GoTools/geometry/ParamSurface.h"


namespace Go {


class SplineSurface;


/// Class that represents the "intersection object" of a parametric
/// surface of dimension 1.

class Param2FunctionInt : public ParamFunctionInt {
public:
    /// Constructor.
    /// \param surf the parametric 1-dimensional surface defining the
    /// object.
    explicit Param2FunctionInt(shared_ptr<ParamSurface> surf);

    /// Constructor.
    /// \param surf the parametric 1-dimensional surface defining the
    /// object.
    /// \param parent the parent object to this object.
    explicit Param2FunctionInt(shared_ptr<ParamSurface> surf,
			       Param2FunctionInt *parent);

    /// Destructor.
    virtual ~Param2FunctionInt() {};

    /// Evaluate the object in the input parameter.
    /// \param res the Point to be returned.
    /// \param par the parameter in which to evaluate. The size of the
    /// array should be equal to numParams().
    virtual void point(Point& res, const double *par) const
    { surf_->point(res, par[0], par[1]); }

    /// Evaluate the object in the input parameter, with the specified
    /// number of derivatives.
    /// \param pt the vector of points to be returned.
    /// \param tpar the parameter in which to evaluate. The size of
    /// the array should be equal to numParams().
    /// \param der the number of derivatives to calculate.
    /// \param from_right if true the evaluation is to be performed
    /// from the right side of the parameter value.
    /// \param resolution tolerance used when determining whether
    /// parameters are located at special values of the parameter
    /// domain (in particualar: knot values in case of spline objects.
    virtual void point(std::vector<Point>& pt, 
		       const double* tpar, 
		       int der,
		       const bool* from_right = 0,
		       double resolution = 1.0e-12) const
    {
	if (from_right == 0) { // User did not specify direction of
			       // differentiation
	    surf_->point(pt, tpar[0], tpar[1], der, true, true);
	} else {
	    surf_->point(pt, tpar[0], tpar[1], der,
			 from_right[0], from_right[1], resolution);
	} 
    }

    /// Return pointer to this object.
    virtual Param2FunctionInt* getParam2FunctionInt();

    /// Return pointer to the parametric surface.
    shared_ptr<ParamSurface> getParamSurface();
    /// Return pointer to the parametric surface.
    shared_ptr<const ParamSurface> getParamSurface() const;

    /// Return pointer to a subpart of the parent surface of this
    /// object.  If a parent surface does not exist, return pointer to
    /// surface in this object.  To reduce numerical noise we go
    /// straight to the source (undivided) surface.
    shared_ptr<ParamSurface>
    getParentParamSurface(RectDomain& domain);
    /// Return pointer to the parent-curve of this object.  If a
    /// parent surface does not exist, return pointer to surface in
    /// this object.  To reduce numerical noise we go straight to the
    /// source (undivided) surface.
    shared_ptr<ParamSurface> getParentParamSurface();

    /// Return an intersection object for the input surface.
    /// This object is used as the parent for the intersection object.
    /// \param surf the parametric surface defining the intersection
    /// object.
    virtual shared_ptr<Param2FunctionInt> 
    makeIntFunction(shared_ptr<ParamSurface> surf);
    
    /// The number of parameters in the object.
    virtual int numParams() const;

    /// Returns the specified isocurve.
    /// \param param_start start parameter for the isocurve.
    /// \param param_end end parameter for the isocurve.
    /// \param isoval the value for the isoparameter.
    /// \param pardir_is_u if 'pardir_is_u' is 'true', then the first
    /// parameter is the running direction and the second parameter is
    /// the isoparameter; vice versa for 'pardir_is_u' equal to
    /// 'false'.
    shared_ptr<ParamCurve> 
    getIsoCurve(double param_start, double param_end,
		double isoval, bool pardir_is_u) const;

    /// Return a curveOnSurface along the current surface in the given
    /// direction and parameter value
    shared_ptr<ParamCurve>
	getConstantParameterCurve(int dir, double par);

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

    /// Return the critical parameter values and inner knots for the
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

    // @@sbr A trimmed sf need not have a startparam/endparam ...
    /// Return the start parameter value in the specified direction.
    /// \param pardir the parameter direction in question. Indexing
    /// starts at 0.
    virtual double startParam(int pardir) const;

    /// Return the end parameter in the specified direction.
    /// \param pardir the parameter direction in question. Indexing
    /// starts at 0.
    virtual double endParam(int pardir) const;

    /// Return true if the specified point lies within eps from the
    /// boundary.
    /// \param par the parameter in which to evaluate. Size of array
    /// should be equal to numParams().
    /// \param eps the tolerance defining boundary neighbourhood.
    virtual bool boundaryPoint(const double* par, double eps) const;

    /// Return the domain of the surface.
    /// \return The returned vector consists of (umin, umax, vmin,
    /// vmax).
    std::vector<double> getMima() const;

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

    /// Return the CompositeBox for the parametric object.
    /// \return The compositeBox for the parametric object.
    virtual CompositeBox compositeBox() const;

    /// Return the boundary objects of this object.
    /// \param bd_objs the boundary objects of this object.
    virtual void 
    getBoundaryObjects(std::vector<shared_ptr<BoundaryFunctionInt> >& bd_objs);

    /// Return the dimension of the geometric space.
    int dimension()
    { return dim_; }
    
    /// Return info on parameter domain which needs special treatment
    /// near a degenerated edge.
    /// \param dir the parameter direction in question. Indexing
    /// starts at 0.
    /// \param deg_edge the degenerate edge
    /// \param threshold a tolerance
    /// \return a safe parameter value that is isolated from the
    /// degenerate value
    double isolateDegPar(int dir, int deg_edge, double threshold);

    // Functions used from coincidence checking

    /// Make sure that the input parameter lies inside the range of
    /// the parametric surface.  Set t equal to umin/vmin if it lies
    /// below umin/vmin, or umax/vmax if it lies above umax/vmax.
    /// \param pardir the parameter direction in question. Indexing
    /// starts at 0.
    /// \param t the input parameter
    void assureInRange(int pardir, double& t);
 
    // Find the interval in which t lies, if t is within tol of an
    // existing knot, t will be changed to that knot.

    /// Return the knot interval for which t lies inside, moving the
    /// value t if it lies close to a knot.
    /// \param pardir the parameter direction in question. Indexing
    /// starts at 0.
    /// \param t the parameter value.
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

    /// Return true if the object is degenerate in the specified
    /// direction.
    /// \param epsge the geometric tolerance defining degeneracy.
    /// \param pardir the parameter direction in question.
    bool isDegenerate(double epsge, int pardir);

    /// Return true if the object is degenerate in the specified
    /// direction and parameter.
    /// \param epsge the geometric tolerance defining degeneracy.
    /// \param pardir the parameter direction in question.
    /// \param par the parameter in which to evaluate. Size of array
    /// should be equal to numParams().
    virtual bool isDegenerate(double epsge, int pardir, double* par);

    /// Return pointer to the parametric surface.
    /// \return The pointer to the parametric surface.
    shared_ptr<ParamSurface> getSurface()
    { return surf_; }

    /// Return pointer to the parametric surface.
    /// \return The pointer to the parametric surface.
    shared_ptr<const ParamSurface> getSurface() const
    { return surf_;}

    /// Return the partial derivatives in the input parameter point.
    /// \param u the u-parameter
    /// \param v the v-parameter
    /// \param deriv_u the partial derivative in the u-direction.
    /// \param deriv_v the partial derivative in the u-direction.
    void derivs(double u, double v,
		Point& deriv_u, Point& deriv_v) const
    {
	surf_->point(temp_point_array_, u, v, 1);
	deriv_u = temp_point_array_[1];
	deriv_v = temp_point_array_[2];
    }

protected:
    // Data members

    shared_ptr<ParamSurface> surf_; // Our 2D-parametric 1D
					   // function (R^2->R).

    int dim_; // Space dimension. 1.

    Param2FunctionInt *parentfunction_; // Probably undivided sf, to
					// avoid num error.

    std::vector<std::pair<double, int> > segment_[2];

    // pardir (0 || 1) is that opposite bdidx (i.e. dir of constant
    // parameter tpar along edge given by bdidx).
    Go::SplineCurve* extractBdCurve(const Go::SplineSurface&, int bdidx,
				    int& pardir, double& tpar) const;

    // Generated data
    // @@sbr Possibly include more of the computed datas. For instance
    // structure conc monotonicity.
    bool bd_deg_[4];  // Degeneracy, sequence left, right, bottom, top
    double deg_tol_;
    RectDomain domain_;
    // Approximating mesh
    mutable int nmesh_[2];
    mutable std::vector<double> mesh_;

    mutable std::vector<Point> temp_point_array_;

};


} // namespace Go


#endif // _PARAM2FUNCTIONINT_H

