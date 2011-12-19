//===========================================================================
//                                                                           
// File: ParamSurfaceInt.h 
//                                                                           
// Created: 
//                                                                           
// Author: oan
//                                                                           
// Revision: $Id: ParamSurfaceInt.h,v 1.58 2007-01-15 10:12:38 vsk Exp $
//                                                                           
// Description:
//                                                                           
//===========================================================================

#ifndef _PARAMSURFACEINT_H
#define _PARAMSURFACEINT_H


#include "GoTools/intersections/ParamGeomInt.h"
#include "GoTools/geometry/ParamSurface.h"


namespace Go {


class AlgObj3DInt;
class ImplicitizeSurfaceAlgo;
class RotatedBox;


/// Class that represents the "intersection object" of a parametric
/// surface.

class ParamSurfaceInt : public ParamGeomInt {
public:
    /// Constructor.
    /// \param surf the parametric surface defining the intersection
    /// object.
    /// \param parent the parent object to this object.
    explicit ParamSurfaceInt(shared_ptr<ParamSurface> surf,
			     ParamGeomInt* parent = 0);
    
    /// Destructor.
    virtual ~ParamSurfaceInt();
    
    /// Evaluate the object in the input parameter.
    /// \param pt the Point to be returned.
    /// \param tpar the parameter in which to evaluate. The size of
    /// the array should be equal to numParams().
    virtual void point(Point& pt, const double* tpar) const 
    { surf_->point(pt, tpar[0], tpar[1]); }

    /// Evaluate the object in the input parameter, with the specified
    /// number of derivatives.
    /// \param pt the vecotr of points to be returned.
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
	if (from_right) {
	    surf_->point(pt, tpar[0], tpar[1], derivs, from_right[0], 
			 from_right[1], resolution);
	} else {
	    surf_->point(pt, tpar[0], tpar[1], derivs, true, true,
			 resolution);
	}
    }

    /// Return a pointer to this object.
    /// \return Pointer to this object.
    virtual ParamSurfaceInt* getParamSurfaceInt();

//     shared_ptr<ParamSurface> getSurface()
//     { return surf_; }
//     shared_ptr<const ParamSurface> getSurface() const
//     { return surf_;}

    /// Return pointer to the parametric surface defining this object.
    /// \return Pointer to the parametric surface defining this
    /// object.
    shared_ptr<ParamSurface> getParamSurface()
    { return surf_; }

    /// Return pointer to the parametric surface defining this object.
    /// \return Pointer to the parametric surface defining this
    /// object.
    shared_ptr<const ParamSurface> getParamSurface() const
    { return surf_;}

    /// Return pointer to a subsurface of the parent surface for this
    /// object.  If no such surface exist, we use the surface of this
    /// object instead.
    /// \param domain the parametric domain of the subsurface.
    /// \return Pointer to a subsurface of the parent surface for this
    /// object.
    shared_ptr<ParamSurface> getParentParamSurface(RectDomain& domain);

    /// Return pointer to the parent surface for this object.  If no
    /// such surface exist, we use the surface of this object instead.
    /// \return Pointer to a subsurface of the parent surface for this
    /// object.
    shared_ptr<ParamSurface> getParentParamSurface();

    /// Return an intersection object for the input surface, using
    /// this object as parent.
    /// \param surf the parametric surface defining the intersection
    /// object.
    virtual shared_ptr<ParamSurfaceInt>
    makeIntObject(shared_ptr<ParamSurface> surf);

    /// Return an intersection object for the input curve, using
    /// parameter parent as parent.
    /// \param crv the parametric curve defining the intersection
    /// object.
    /// \param parent the parent to the created intersection object.
    virtual shared_ptr<ParamCurveInt> 
    makeIntCurve(shared_ptr<ParamCurve> crv, ParamGeomInt* parent);

    /// The number of parameters in the object.
    virtual int numParams() const
    { return 2;}

//     int checkCoincidence(double start, double end,
// 			 ParamSurfaceInt *other,
// 			 double other_start, double other_end);

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
    /// \param dir the direction of the constant parameter
    /// curve. Indexing starts at 0.
    /// \param par the isoparameter for the curve.
    shared_ptr<ParamCurve>
    getConstantParameterCurve(int dir, double par);

    shared_ptr<ParamCurve>
    getConstantParameterCurve(int dir, double par, double tmin, double tmax);

    void 
	minimumAlongCurve(int dir, double par, double tmin, double tmax,
			  Point& pnt, double& minpar, Point& minval, 
			  double& mindist);

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

    virtual bool canDivideTinyTriang(int pardir);

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
    virtual double startParam(int pardir) const;

    /// Return the end parameter in the specified direction.
    /// \param pardir the parameter direction in question. Indexing
    /// starts at 0.
    virtual double endParam(int pardir) const;

    /// Return true if the specified point lies within eps from the
    /// boundary.
    /// \param par the parameter in which to evaluate. Size of array
    /// should be equal to numParams().
    /// \param eps the tolerance defining the boundary neighbourhood.
    virtual bool boundaryPoint(const double* par, double eps) const;

    /// Return the domain of the surface.
    /// \return The returned vector consists of (umin, umax, vmin,
    /// vmax).
    std::vector<double> getMima() const;

    /// Check if parameter point lies close to a corner of the
    /// parameter domain.
    /// \param par the input parameter point.
    /// \param epspar the parametric tolerance defining the
    /// neighbourhood.
    virtual bool inCorner(const double *par, double epspar) const;

    bool atDegenerateBd(const double *par, double epsge, double epspar);

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

    /// Return the subsurface(s) on the input domain.  For a trimmed
    /// surface there may be more than one surface, for a surface
    /// defined on a rectangular domain there will be only one.
    /// \param from_upar start parameter in the first parameter
    /// direction.
    /// \param from_vpar start parameter in the second parameter
    /// direction.
    /// \param to_upar end parameter in the first parameter direction.
    /// \param to_vpar end parameter in the second parameter
    /// direction.
    /// \param fuzzy allowed alteration of an input parameter value.
    /// Typically this applies to a spline surface, where we do not
    /// want knots to lie too close whilst not being equal.
    virtual std::vector<shared_ptr<ParamSurfaceInt> >
    subSurfaces(double from_upar, double from_vpar,
		double to_upar, double to_vpar,
		double fuzzy);

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

    /// Check if the object is periodic in the specified direction.
    /// Analyze periodicity of surface based on number of repeating
    /// knots and control points. The return value is -1 if the
    /// surface ends are disjoint, otherwise k if cv is C^k
    /// continuous. These are sufficient but not necessary conditions
    /// for periodicity, so it is possible that a higher degree of
    /// periodicity exists.  Should not be called on this layer,
    /// should be overruled by inherited class.
    /// \param pardir the parameter direction in question.
    /// \return -1 if the curve ends are disjoint, or k if the surface
    /// is proven to be C^k continuous.
    virtual int checkPeriodicity(int pardir) const;

    /// The dimension of the geometric space.    
    virtual int dimension() const
    { return dim_; }

    /// Return the rectangular domain of the surface.
    /// \return Return the rectangular domain of the surface.
    virtual const RectDomain& getDomain() const
    { return domain_;}

    /// Return a rectangular domain of the surface which is slightly reduced
    /// in the degenerate case
    /// \return Return the rectangular domain of the surface.
    const RectDomain& getDegDomain(double epsge);

    /// Estimates if the current surface is simple enough for a
    /// singularity iteration. Checks the span of the normal cone and
    /// the size of the surface
    /// \return True if the surface is characterized as simple.
    virtual bool isSimple();

    /// Verify whether the object is a spline.
    /// \return True if the object is a spline.
    virtual bool isSpline();

    /// Check if a curve is an iso parametric curve in the current surface.
    /// NB! Only valid for splines
    virtual bool isIsoParametric(ParamCurveInt *curve, int dir, double par,
				 double ptol, double tol)
	{ return false; }  // Default behaviour

    /// We try to treat problems which will never result in a simple
    /// case by shrinking the domain slightly, resulting in smaller
    /// cones.  This is useful for scenarios where the normals are
    /// parallell in a boundary point.
    /// \param axis1 first vector defining a projection plane
    /// \param axis2 second vector defining a projection plane
    /// \return The optimized cone angle
    virtual double getOptimizedConeAngle(Point& axis1, Point& axis2);

    /// Find the knot intervals for which u and v lie inside, moving
    /// the value u or v if they lie close to a knot.
    /// \param u the \a u parameter value
    /// \param v the \a v parameter value
    /// \param utol the parametric tolerance deciding if the input
    /// parameter u should be moved.
    /// \param vtol the parametric tolerance deciding if the input
    /// parameter v should be moved.
    virtual void knotIntervalFuzzy(double& u, double&v,
				   double utol, double vtol) const; 

    /// Return the value of the knot next to the input parameter par.
    /// \param dir the parameter direction in question. Indexing
    /// starts at 0.
    /// \param par the parameter value
    /// \param forward if true we return the closest knot to the
    /// right, otherwise the closest knot to the left.
    /// \param tol the tolerance to determine if \a par is already
    /// located on the start of the next segment.
    /// \return The knot closest to the input parameter.
    virtual double nextSegmentVal(int dir, double par,
				  bool forward, double tol) const;

    // Degeneracy
    /// Verfify whether the object is degenerate.
    /// \param epsge the geometric tolerance defining degeneracy.
    bool isDegenerate(double epsge);

    /// Verfify whether the object is degenerate in the specified
    /// direction.
    /// \param epsge the geometric tolerance defining degeneracy.
    /// \param dir the parameter direction in question.
    virtual int isDegenerate(double epsge, int dir);

    /// Verfify whether the object is degenerate in the specified
    /// direction and parameter.
    /// \param epsge the geometric tolerance defining degeneracy.
    /// \param dir the parameter direction in question.
    /// \param par the parameter in which to evaluate. Size of array
    /// should be equal to numParams().
    virtual bool isDegenerate(double epsge, int dir, double *par);

    /// Return info on parameter domain which needs special treatment
    /// near a degenerated edge.
    /// \param dir the parameter direction in question. Indexing
    /// starts at 0.
    /// \param deg_edge the degenerate edge
    /// \param threshold a tolerance
    /// \return a safe parameter value that is isolated from the
    /// degenerate value
    double isolateDegPar(int dir, int deg_edge, double threshold, double *deg_factor=NULL);

    /// Set dege_triang info.
    /// Variable useful for surfaces with degenerated edge(s).
    void setDegTriang()
    { deg_triang_ = true; }

    /// Get deg_triang info.
    /// Variable useful for surfaces with degenerated edge(s).
    bool getDegTriang()
    { return deg_triang_; }
    
    /// Estimate the parameter value of a surface a specified distance
    /// from a given edge.
    /// \param dir the parameter direction in question.
    /// \param atstart true if we are the start of the parameter
    /// interval, false otherwise.
    /// \param tol influencing how far from the edge we should move.
    /// \return the computed parameter value.
    double getParOffBd(int dir, bool atstart, double tol) const;

    /// Return the coefficients necessary to calculate the first
    /// fundamental form of the surface in the parameter point.
    /// \param u the first parameter value.
    /// \param v the second parameter value.
    /// \param u_from_right if true then evaluate from the right,
    /// otherwise from the left.
    /// \param v_from_right if true then evaluate from the right,
    /// otherwise from the left.
    /// \param E the first coefficient for the first fundamental form.
    /// \param F the second coefficient for the first fundamental
    /// form.
    /// \param G the third coefficient for the first fundamental form.
    void  first_fund_form(double u, 
			  double v, 
			  bool u_from_right, 
			  bool v_from_right, 
			  double& E, 
			  double& F, 
			  double& G) const
    {
	// Getting first derivatives
	surf_->point(temp_point_array_, u, v, 1,
		     u_from_right, v_from_right); 
	
	E = temp_point_array_[1] * temp_point_array_[1];
	F = temp_point_array_[1] * temp_point_array_[2];
	G = temp_point_array_[2] * temp_point_array_[2];
    }

    /// Return the coefficients necessary to calculate the second
    /// fundamental form of the surface in the parameter point.
    /// \param u the first parameter value.
    /// \param v the second parameter value.
    /// \param u_from_right if true then evaluate from the right,
    /// otherwise from the left.
    /// \param v_from_right if true then evaluate from the right,
    /// otherwise from the left.
    /// \param L the first coefficient for the second fundamental
    /// form.
    /// \param M the second coefficient for the second fundamental
    /// form.
    /// \param N the third coefficient for the second fundamental
    /// form.
    void second_fund_form(double u, 
			  double v, 
			  bool u_from_right,
			  bool v_from_right,
			  double& L, 
			  double& M, 
			  double& N) const
    {
	// Getting second derivatives
	surf_->point(temp_point_array_, u, v, 2,
		     u_from_right, v_from_right); 
	// Getting normal (we suppose that the surface contains a
	// continuous first derivative...)
	surf_->normal(temp_point_array_[0], u, v); // Setting first
					           // point to normal
					           // vec

	L = temp_point_array_[3] * temp_point_array_[0];  // d2r/du2 . N
	M = temp_point_array_[4] * temp_point_array_[0];  // d2r/dudv . N
	N = temp_point_array_[5] * temp_point_array_[0];  // d2r/dv2 . N
    }

    /// Return the partial derivatives in the input parameter point.
    /// \param u the first parameter value.
    /// \param v the second parameter value.
    /// \param deriv_u the partial derivative in the u-direction.
    /// \param deriv_v the partial derivative in the u-direction.
    /// \param from_right_1 if true then calculate from the left int
    /// the first parameter parameter direction, otherwise from the
    /// right.
    /// \param from_right_2 if true then calculate from the left int
    /// the second parameter parameter direction, otherwise from the
    /// right.
    void derivs(double u, 
		double v, 
		Point& deriv_u, 
		Point& deriv_v,
		bool from_right_1 = true,
		bool from_right_2 = true) const
    {
	surf_->point(temp_point_array_, u, v, 1,
		     from_right_1, from_right_2);
	deriv_u = temp_point_array_[1];
	deriv_v = temp_point_array_[2];
    }

    /// Calculate the normal in the specified parameter point.
    /// \param u the first parameter value.
    /// \param v the second parameter value.
    /// \param normal the calculated normal in the parameter point.
    /// \param from_right_1 if true then calculate from the left int
    /// the first parameter parameter direction, otherwise from the
    /// right.
    /// \param from_right_2 if true then calculate from the left int
    /// the second parameter parameter direction, otherwise from the
    /// right.
    void normal(double u, 
		double v, 
		Point& normal,
		bool from_right_1 = true,
		bool from_right_2 = true) const
    {
	surf_->point(temp_point_array_, u, v, 1,
		     from_right_1, from_right_2);
	normal = temp_point_array_[1] % temp_point_array_[2];
    }

    /// Create a box containing the geometric sample mesh in the input
    /// coordinate system.
    /// \param axis the axis defining the coordinate system. The size
    /// of vector axis is 2, the last point is given as the cross
    /// product axis[0]%axis[1].
    /// \return the rotated box
    virtual RotatedBox getRotatedBox(std::vector<Point>& axis) const;

//     virtual void dumpToStream(std::ostream& os) const
//     {
// 	surf_->writeStandardHeader(os);
// 	surf_->write(os);
//     }

    /// Use the corner points of the parameter domain to create axes.
    /// \param axis1 the first axis.
    /// \param axis2 the second axis.
    void axisFromCorners(Point& axis1, Point& axis2) const;

    /// Requests from selfintersection computation.  Split at G1
    /// discontinuities
    /// \param angtol angular tolerance defining G1 discontinuity
    /// \param subG1 vector of subdivided patches that are G1
    virtual void splitAtG0(double angtol,
			   std::vector<shared_ptr<ParamSurfaceInt> >& subG1);

    /// Iterate for a surface singularity.
    /// \param eps numerical tolerance
    /// \param sing_par parameters of the singularity
    /// \param sing_pt the singular point
    /// \param sing_val a number that measures the quality of the
    /// singularity
    /// \param seed initial guess for singularity in parameter plane
    void getSingularity(double eps, double sing_par[], Point& sing_pt,
			double& sing_val, double *seed);

    /// Check if the current surface can selfintersect.  Uses normal
    /// cone and tangent cones.
    /// \param epsge a geometric tolerance
    bool canSelfIntersect(double epsge) const;

    /// Returns the normal surface corresponding to this surface.
    /// \return Pointer to the SplineSurface defining the normal
    /// surface.  Not implemented for a general parametric surface.
    virtual shared_ptr<ParamSurfaceInt> getNormalSurface() const
    {
	shared_ptr<ParamSurfaceInt> dummy;
	return dummy;
    }

    /// Check whether the objects fulfills the requirements to
    /// implicitize.
    /// \return Return true if we may implicitize the object.
    virtual bool canImplicitize()
    { return false; }

    /// Implicitize the object.
    /// \param tol the geoemtric tolerance for the implicitization.
    /// \return True if the implicitization was a success.
    virtual bool implicitize(double tol)
    { return false; }

    /// Get the implicit representation of the object.
    /// Garbage is returned if we are not able to implicitize.
    /// \param tol geometric tolerance for the implicitization
    /// procedure.
    /// \param tol2 geometric estimate for the accuracy of the
    /// implicitized object.  Not yet in use!!!
    /// \param alg_obj_3d_int the algebraic object containing the
    /// implicitized surface.
    /// \return True if the implicitization was a success.
    // @@sbr Todo (tol2)!!!
    virtual bool getImplicit(double tol, double& tol2,
			     AlgObj3DInt& alg_obj_3d_int)
    { return false; }

//     /// Find the maximum curvature in a number of isocurves.
//     /// \param dir_is_u the direction of the isocurves.
//     /// \param nmb_params the number of isocurves to check.
//     /// \param iso_from start parameter of domain of interest in the
//     /// isocurves.
//     /// \param iso_to end parameter of domain of interest in the
//     /// isocurves.
//     /// \param guess_param initial guess for the parameter value with
//     /// the highest curvature.
//     /// \param max_curv_params the parameter values for the maximum
//     /// curvature in each isocurve.
//     /// \param max_curvatures the maximum curvature in each isocurve,
//     /// corresponding the parameters in max_curv_params.
     void maxCurvatures(bool dir_is_u, int nmb_params, double iso_from,
		       double iso_to, double guess_param,
 		       std::vector<double>& max_curv_params,
 		       std::vector<double>& max_curvatures);


protected:
    // Data members
    shared_ptr<ParamSurface> surf_;
    int dim_; // Space dimension.
    double deg_tol_;
    bool deg_triang_;

    std::vector<std::pair<double, int> > segment_[2];

   // Generated data
    bool bd_deg_[4];  // Degeneracy, sequence left, right, bottom, top
    RectDomain domain_;
    RectDomain deg_domain_;   // Reduced domain for degenerate surfaces used in closest
                              // point computations
    // Approximating mesh
    mutable int nmesh_[2];
    mutable std::vector<double> mesh_;

    mutable DirectionCone cone_;

    mutable bool lw_set_;
    mutable double length_[2], wiggle_[2];

//     const static int nder_ = 2; // Order of derivatives to be
// 				// calulated.
//     double delta_;
//     std::vector<Point> ft_;     // Partial derivatives up to order
// 				// nder
//     std::vector<Point> gs_;     // Partial derivatives up to order
// 				// nder
    mutable std::vector<Point> temp_point_array_;

    // Implicit objects
    shared_ptr<ImplicitizeSurfaceAlgo> impl_sf_algo_;
    double implicit_tol_; // Error tolerance used when creating
			  // impl_sf_algo_.
    int impl_deg_;
    shared_ptr<AlgObj3DInt> implicit_obj_;
    double implicit_err_;
 
    
 private:
    void makeMesh(int size1, int size2) const;

    void computeDegDomain(double aepsge);

};


} // namespace Go


#endif // _PARAMSURFACEINT_H
