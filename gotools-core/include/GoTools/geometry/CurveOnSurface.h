/*
 * Copyright (C) 1998, 2000-2007, 2010, 2011, 2012, 2013 SINTEF ICT,
 * Applied Mathematics, Norway.
 *
 * Contact information: E-mail: tor.dokken@sintef.no                      
 * SINTEF ICT, Department of Applied Mathematics,                         
 * P.O. Box 124 Blindern,                                                 
 * 0314 Oslo, Norway.                                                     
 *
 * This file is part of GoTools.
 *
 * GoTools is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Affero General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version. 
 *
 * GoTools is distributed in the hope that it will be useful,        
 * but WITHOUT ANY WARRANTY; without even the implied warranty of         
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the          
 * GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public
 * License along with GoTools. If not, see
 * <http://www.gnu.org/licenses/>.
 *
 * In accordance with Section 7(b) of the GNU Affero General Public
 * License, a covered work must retain the producer line in every data
 * file that is created or manipulated using GoTools.
 *
 * Other Usage
 * You can be released from the requirements of the license by purchasing
 * a commercial license. Buying such a license is mandatory as soon as you
 * develop commercial activities involving the GoTools library without
 * disclosing the source code of your own applications.
 *
 * This file may be used in accordance with the terms contained in a
 * written agreement between you and SINTEF ICT. 
 */

#ifndef _CURVEONSURFACE_H
#define _CURVEONSURFACE_H

#include "GoTools/geometry/ParamSurface.h"
#include <memory>
#include "GoTools/geometry/ParamCurve.h"
#include "GoTools/geometry/SplineCurve.h"
#include "GoTools/utils/config.h"

namespace Go
{

// class SplineCurve;

    /** \brief A curve living on a parametric surface. It either has got 
	information about the curve in geometry space and in the parameter
	domain of the surface or the ability to compute the other representation 
	given one. The curve may have information on whether it is a constant
	parameter or boundary curve on the surface.
     *
     */

class GO_API CurveOnSurface : public ParamCurve
{
public:
    /// Define an empty CurveOnSurface that can be assigned or \ref read()
    /// into
    CurveOnSurface();

    /// Construct a CurveOnSurface by specifying the surface and either a space curve
    /// or a curve in the parameter plane.  Does not clone any of the input, just 
    /// sets (smart) pointers.
    /// \param surf pointer to the underlying surface
    /// \param curve pointer to the curve specifying the CurveOnSurface.  This
    ///              curve may either be in the parametric domain of the surface
    ///              (2D curve), or a space curve that the user has assured 
    ///              to be coincident with the surface.  'preferparameter' specifies
    ///              which kind of curve this is.
    /// \param preferparameter if this is set to 'true', then 'curve' is assumed
    ///                        to be a curve in the parametric domain of the surface.
    ///                        Otherwise, it is assumed to be a space (3D) curve.
  CurveOnSurface(shared_ptr<ParamSurface> surf,
		     shared_ptr<ParamCurve> curve,
                     bool preferparameter);
    
    /// Constructor for constant parameter curves
  CurveOnSurface(shared_ptr<ParamSurface> surf,
		 shared_ptr<ParamCurve> curve,
		 int constdir, double constpar, int boundary);

  /// Constructor for constant parameter curves. The curve is to
  /// be restricted to the parameter interval [par1, par2].
  CurveOnSurface(shared_ptr<ParamSurface> surf,
		 int constdir, double constpar, 
		 double par1, double par2, int boundary);
    
    /// Constructor with all information
  CurveOnSurface(shared_ptr<ParamSurface> surf,
		 shared_ptr<ParamCurve> pcurve,
		 shared_ptr<ParamCurve> spacecurve,
		 bool preferparameter, int ccm,
		 int constdir, double constpar, int boundary,
		 bool same_orientation);

    /// Construct a CurveOnSurface by specifying the surface, the space curve and the
    /// curve in the parameter plane.  The arguments are checked for consistency, and if
    /// they appear incoherent, an exception will be thrown.
    /// Does not clone any of the input, just sets (smart) pointers.
    /// \param surf pointer to the underlying surface
    /// \param pcurve pointer to the curve representing the CurveOnSurface in the
    ///               parametric domain of the surface.
    /// \param spacecurve pointer to the curve that is the spatial (3D) representation 
    ///                   of the CurveOnSurface.
    /// \param preferparameter specify whether the parametric curve or the space
    ///                        curve are preferred for internal computations.

    /// \param ccm curve creation method. Specify how curve was created.
    ///            0 = undefined, 1 = projection onto surface,
    ///            2 = intersection of two surfaces,
    ///            3 = isoparametric curve (i.e. either a u-curve or a v-curve).
    CurveOnSurface(shared_ptr<ParamSurface> surf,
		   shared_ptr<ParamCurve> pcurve,
		   shared_ptr<ParamCurve> spacecurve,
		   bool preferparameter, int ccm = 0);

    /// Copy constructor.  The copy constructor will not clone() the underlying
    /// surface, but it will clone() both the parametric and the spatial curve.
    /// \param surface_curve the CurveOnSurface to copy into 'this' CurveOnSurface.
    CurveOnSurface(const CurveOnSurface& surface_curve);

    /// Assignment operator.  Like the copy constructor, the assignment operator
    /// clone()s the curves, and not the surface.
    /// \param other the CurveOnSurface to copy into 'this' CurveOnSurface.
    CurveOnSurface& operator= (const CurveOnSurface& other);
    
    /// Destructor.
    /// Trivial because memory is managed by shared_ptr.
    virtual ~CurveOnSurface();


    // inherited from Streamable
    virtual void read (std::istream& is);
    virtual void write (std::ostream& os) const;

    // inherited from GeomObject
    /// Axis align box surrounding this object
    /// Computed with respect to the space curve if this one exists,
    /// otherwise the underlying surface
    virtual BoundingBox boundingBox() const;

    /// Cone surrounding the set of tangent directions corresponding to
    /// the surface. Only computed if the space curve exists
    virtual DirectionCone directionCone() const;

    /// Dimension of geometry space
    virtual int dimension() const;

    /// Type of current object. 
    virtual ClassType instanceType() const;
    static ClassType classType()
    { return Class_CurveOnSurface; }
    
    // The function clone() calls the copy constructor,
    // so clone() also makes deep copies of the curves,
    // but not the underlying surface.
// #ifdef _MSC_VER
// #if _MSC_VER < 1300
//     virtual GeomObject* clone() const
//       { return new CurveOnSurface(*this); }
// #else
//     virtual CurveOnSurface* clone() const
//       { return new CurveOnSurface(*this); }
// #endif // _MSC_VER < 1300
// #else
    virtual CurveOnSurface* clone() const
    { return new CurveOnSurface(*this); }
// #endif

    // inherited from ParamCurve
    virtual void point(Point& pt, double tpar) const;

    /// Inherited from \ref ParamCurve. Only works for 'derivs' = 0 or 1.
    /// \see ParamCurve::point()
    virtual void point(std::vector<Point>& pts, 
		       double tpar,
		       int derivs, bool from_right = true) const;
    
    // inherited from ParamCurve
    virtual double startparam() const;

    // inherited from ParamCurve
    virtual double endparam() const;

    // inherited from ParamCurve
    virtual void reverseParameterDirection(bool switchparam = false);

    // inherited from ParamCurve
    virtual void setParameterInterval(double t1, double t2);

    // inherited from ParamCurve
    virtual SplineCurve* geometryCurve();

    // inherited from ParamCurve
    virtual bool isDegenerate(double degenerate_epsilon);


// #if ((_MSC_VER > 0) && (_MSC_VER < 1300))
//     // inherited from ParamCurve
//     virtual ParamCurve* subCurve(double from_par, double to_par,
// 				   double fuzzy = DEFAULT_PARAMETER_EPSILON) const;
// #else
    // inherited from ParamCurve
    virtual CurveOnSurface* subCurve(double from_par, double to_par,
				       double fuzzy =
				       DEFAULT_PARAMETER_EPSILON) const;
// #endif

    /// Split curve in a specified parameter value
    virtual
      std::vector<shared_ptr<ParamCurve> > 
      split(double param,
	    double fuzzy = DEFAULT_PARAMETER_EPSILON) const; 

   // inherited from ParamCurve
    virtual void closestPoint(const Point& pt,
			      double         tmin,
			      double         tmax,
			      double&        clo_t,
			      Point&       clo_pt,
			      double&        clo_dist,
			      double const   *seed = 0) const;

    // Inherited from ParamCurve
    // Compute the total length of this curve
    virtual double length(double tol);

    // inherited from ParamCurve.  NB: does not check whether the resulting ParamCurve
    // stays inside parameter domain (or that the space curve stays on surface).
    virtual void appendCurve(ParamCurve* cv, bool reparam=true);

    // inherited from ParamCurve.  NB: does not check whether the resulting ParamCurve
    // stays inside parameter domain (or that the space curve stays on surface).
    virtual void appendCurve(ParamCurve* cv,
			     int continuity, double& dist, bool reparam=true);

    /// Set the underlying surface to the one pointed to by the argument
    /// \param surface the pointer to the surface we will set as underlying for this
    ///                CurveOnSurface.
    void setUnderlyingSurface(shared_ptr<ParamSurface> surface)
    {surface_ = surface;}

    /// Replace the curves describing the curve on surface. The curve preference
    /// is not changed
    void setCurves(shared_ptr<ParamCurve> spacecurve,
		      shared_ptr<ParamCurve> parametercurve)
    {
      spacecurve_ = spacecurve;
      pcurve_ = parametercurve;
    }

    /// Replace the space curve corresponding to this curve on surface curve.
    /// Use with care!
    void setSpaceCurve(shared_ptr<ParamCurve> spacecurve)
    {
      spacecurve_ = spacecurve;
    }

    /// Replace the parameter curve corresponding to this curve on surface curve.
    /// Used for instance in relation to reparameterizations of the related surface.
    /// Use with care!
    void setParameterCurve(shared_ptr<ParamCurve> parametercurve)
    {
      pcurve_ = parametercurve;
    }

    /// Remove parameter curve information
    void unsetParameterCurve()
    {
      pcurve_.reset();
      prefer_parameter_ = false;
    }

    /// Generate parameter curve if it doesn't exist
    bool ensureParCrvExistence(double epsgeo,
			       const RectDomain* domain_of_interest = NULL);

    /// Make parameter curve between given end parameters
    bool makeParameterCurve(double tol, const Point& par1, const Point& par2 );

    /// Translate an existing parameter curve in a given direction
    /// Return value: Whether or not the tranlation was possible
    bool translateParameterCurve(const Point& dir);

    /// Translate an existing parameter curve in a given direction and swap sign
    /// is specified
    /// Return value: Whether or not the tranlation was possible
    bool translateSwapParameterCurve(const Point& dir, double sgn, int pdir);

    /// Represent the parameter curve associated with this
    /// surface curve on a different domain
    bool setDomainParCrv(double umin, double umax, double vmin, double vmax,
			 double uminprev, double umaxprev,
			 double vminprev, double vmaxprev);

    /// Generate space curve if it doesn't exist
    bool ensureSpaceCrvExistence(double tol);

    /// Set curve information
    void setCurveTypeInfo(int ccm)
    {
      ccm_ = ccm;
    }

    /// Update constant parameter curves
    bool updateIsoCurves();

    /// Update constant parameter curves. New constant parameter curve 
    /// information is provided. Used in for instance change of parameter
    /// directions in the associated surface
    bool updateIsoCurves(int constdir, double constpar, int boundary);

    /// Update curves. Reapproximate the dependent curve within the tolerance epsge.
    bool updateCurves(double epsge);

    /// Update curves
    bool updateCurves(Point vx1, Point vx2, double epsge);

    /// Inherited from \ref ParamCurve.  If the parametric curve is set to be the
    /// 'prefered' one, this function will return the next segment value for the 
    /// parametric curve; otherwise it will return the next segment value for the 
    /// spatial 3D curve.
    /// See also \ref ParamCurve::nextSegmentVal()
    virtual double nextSegmentVal(double par, bool forward, double tol) const;


    /// Get a shared pointer to the underlying surface
    /// \return a shared pointer to the underlying surface
    shared_ptr<ParamSurface> underlyingSurface()
    { return surface_; }

    /// Get a shared pointer to the curve in the parameter domain.
    /// \return a shared pointer to the curve in the parameter domain
    shared_ptr<ParamCurve> parameterCurve()
      { return pcurve_; }

    /// Get a shared pointer to the space curve
    /// \return a shared pointer to the space curve.
    shared_ptr<ParamCurve> spaceCurve()
    { return spacecurve_; }

    /// Get a constant, shared pointer to the underlying surface
    /// \return a const-pointer to the underlying surface
    shared_ptr<const ParamSurface> underlyingSurface() const
    { return surface_; }

    /// Get a constant, shared pointer to the curve in the parameter domain.
    /// \return a const-pointer to the curve in the parameter domain.
    shared_ptr<const ParamCurve> parameterCurve() const
    { return pcurve_; }

    /// Get a constant, shared pointer to the space curve
    /// \return a const-shared pointer to the space curve.
    shared_ptr<const ParamCurve> spaceCurve() const
    { return spacecurve_; }
    
    /// Query whether the parameter curve or the space curve is prefered for computation
    /// in this object.
    bool parPref() const
    { return prefer_parameter_; }

    void setParPref(bool prefer) 
    { prefer_parameter_ = prefer; }

     /// Get the rectangle enclosing the underlying surface's parametric domain.
    /// \return the RectDomain for the underlying surface.
    RectDomain containingDomain() const;

    /// Fetch the surface parameter corresponding to a curve parameter
    Point faceParameter(double crv_par,
			const RectDomain* domain_of_interest = NULL) const;

    /// Check if this curve lies at a surface boundary, in that case which one
    /// -1 = none, 0 = umin, 1 = umax, 2 = vmin,  3 = vmax
    /// NB! same_orientation is set only if the curve is defined to lie
    ///  at constant parameter curve. Set to 'false' if direction is opposite
    ///  that of the surface.
    int whichBoundary(double tol, bool& same_orientation) const;

    /// Check if the curve is a constant parameter curve and marked as such
    bool isConstantCurve() const
    {
      return (constdir_ > 0);
    }

    /// Check if the curve is a constant parameter curve with regard to the 
    /// parameter tol
    bool isConstantCurve(double tol, int& pardir, double& parval) const;

    /// Legality test regarding the consistence of geometry and parameter curve.
    /// Check if the curves have got the same parameter domain
    bool sameParameterDomain() const;

    /// Legality test regarding the consistence of geometry and parameter curve.
    /// Check if the two curves have the same orientation
    bool sameOrientation() const;
    
    /// Legality test regarding the consistence of geometry and parameter curve.
    /// Check if the two curves describe the same trace with respect to the
    /// tolerance tol. Pointwise check.
    bool sameTrace(double tol, int nmb_sample = 5) const;

    /// Legality test regarding the consistence of geometry and parameter curve.
    /// Check if the two curves represent identical curves with respect to the
    /// tolerance tol. Pointwise check.
    bool sameCurve(double tol, int nmb_sample = 5) const;

    // Make sure that domain and orientation correspond. Routine does
    // not check difference between trace of parameter curve and space
    // curve.
    void makeCurvesConsistent(bool prefer_parameter);

    /// Return curve creation method
    /// 0 = undefined, 1 = projection onto surface,
    /// 2 = intersection of two surfaces,
    /// 3 = isoparametric curve (i.e. either a u-curve or a v-curve).
    int curveCreationMethod() const;

    /// Maximum distance between curve represented by parameter curve and
    /// space curve in a number of sampling points
    double maxTraceDiff(int nmb_sample = 5) const;

    /// Check if the curve is axis rotational. Only true if a connection
    /// to an axis rotational elementary curve exist
    /// The axis and rotational angle is only specified if the curve
    /// is actually rotational
    virtual bool isAxisRotational(Point& centre, Point& axis, Point& vec,
				  double& angle);

    /// Check if the curve is linear
    virtual bool isLinear(Point& dir, double tol);

   /// Check if the curve lies in a plane passing through a given axis
    virtual bool isInPlane(const Point& loc, const Point& axis,
			   double eps, Point& normal) const;

    /// Check if the curve lies in a plane with a given normal
    virtual bool isInPlane(const Point& norm,
			   double eps, Point& pos) const;

private:
    /// The underlying surface
    shared_ptr<ParamSurface> surface_;
    /// The 2D curve in the parameter domain of the surface.
    /// May point to null.
    shared_ptr<ParamCurve> pcurve_;
    /// An instance of the curve in the surface. May point to null.
    shared_ptr<ParamCurve> spacecurve_;
    /// Which representation to prefer if both exist
    bool prefer_parameter_;
    /// Specify how curve was created. Useful if mismatch between
    /// pcurve_ and spacecurve_.
    /// 0 = undefined, 1 = projection onto surface,
    /// 2 = intersection of two surfaces,
    /// 3 = isoparametric curve (i.e. either a u-curve or a v-curve).
    int ccm_;

    /// More detailed specification of a constant parameter curve
    int constdir_; //1=u_parameter constant, 2=v_parameter constant
    double constval_;  // Value of constant parameter
    int at_bd_;  // -1=no, 0=umin, 1=umax, 2=vmin, 3=vmax
    bool same_orientation_;  // Indicates of the constant parameter curve is
    // oriented the same way as the surface in this parameter direction

    // Fixes performed to get a legal curve on surface
    mutable int fix_performed_;

    // Assuming the surface is closed with multiple param points
    // corresponding to a point on the space curve. Pick the one which
    // is in correspondance with the space curve tangent. If the
    // tangent folows the surface seem nothing is done.
    void pickParamPoint(std::vector<Point>& par_candidates,
			double tpar, double epsgeo);

    // The returned point is a (2D) tangent in the parameter domain.
    Point projectSpaceCurveTangent(const Point& par_pt, double tpar);

};


} // namespace Go

#endif // _CURVEONSURFACE_H

