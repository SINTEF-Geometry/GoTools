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

#ifndef _SPLINESURFACE_H
#define _SPLINESURFACE_H


#include "GoTools/geometry/ParamSurface.h"
#include "GoTools/geometry/BsplineBasis.h"
#include "GoTools/geometry/RectDomain.h"
#include "GoTools/utils/ScratchVect.h"
#include "GoTools/utils/config.h"

namespace Go
{

class Interpolator;
class SplineCurve;
class DirectionCone;
class ElementarySurface;

/// Structure for storage of results of grid evaluation of the basis function of a spline surface.
/// Positional evaluation information in one parameter value
struct BasisPtsSf
{
  /// Parameter tripple in which the basis functions are evaluated
  double param[2]; 
  /// Index of the knot interval where the parameter value is situated for all
  /// parameter directions. The indices of the non-zero basis functions are
  /// left_idx[i]-order[i]+1, ..., left_idx[i] for i=0,1
  int left_idx[2]; 
  /// The value of all basis functions, size equal to (degree_u+1)*(degree_v+1)*(degree_w+1)
    std::vector< double > basisValues; 

    void preparePts(double u, double v, int idx_u, int idx_v, int size)
	{
	    param[0] = u;
	    param[1] = v;
	    left_idx[0] = idx_u;
	    left_idx[1] = idx_v;
	    basisValues.resize(size);
	}
};

/// Structure for storage of results of grid evaluation of the basis function of a spline surface.
/// Position and first derivatives
struct BasisDerivsSf
{
  /// Parameter double in which the basis functions are evaluated
  double param[2];
  /// Index of the knot interval where the parameter value is situated for all
  /// parameter directions. The indices of the non-zero basis functions are
  /// left_idx[i]-order[i]+1, ..., left_idx[i] for i=0,1
  int left_idx[2];
  /// The value of all basis functions, size equal to (degree_u+1)*(degree_v+1)
  std::vector< double > basisValues; 
  /// the derivative of all basis functions in u direction, same size as previous
  std::vector< double > basisDerivs_u; 
  /// the derivative of all basis functions in v direction, same size as previous
  std::vector< double > basisDerivs_v;

    void prepareDerivs(double u, double v, int idx_u, int idx_v, int size)
	{
	    param[0] = u;
	    param[1] = v;
	    left_idx[0] = idx_u;
	    left_idx[1] = idx_v;
	    basisValues.resize(size);
	    basisDerivs_u.resize(size);
	    basisDerivs_v.resize(size);
	}
};

/// Structure for storage of results of grid evaluation of the basis function of a spline surface.
/// Position, first and second derivatives
struct BasisDerivsSf2
{
  /// Parameter double in which the basis functions are evaluated
  double param[2];   
  /// Index of the knot interval where the parameter value is situated for all
  /// parameter directions. The indices of the non-zero basis functions are
  /// left_idx[i]-order[i]+1, ..., left_idx[i] for i=0,1
  int left_idx[2];   
  /// The value of all basis functions, size equal to (degree_u+1)*(degree_v+1)
  std::vector< double > basisValues; 
  
  /// the derivative of all basis functions in u direction, same size as previous
  std::vector< double > basisDerivs_u;
  /// the derivative of all basis functions in v direction, same size as previous
  std::vector< double > basisDerivs_v;

  /// the second derivative of all basis functions twice in u direction, same size as previous
  std::vector< double > basisDerivs_uu;
  /// the second derivative of all basis functions in u and v direction, same size as previous
  std::vector< double > basisDerivs_uv;
  /// the second derivative of all basis functions twice in v direction, same size as previous
    std::vector< double > basisDerivs_vv;

    void prepareDerivs(double u, double v, int idx_u, int idx_v,
		       int size)
	{
	    param[0] = u;
	    param[1] = v;
	    left_idx[0] = idx_u;
	    left_idx[1] = idx_v;
	    basisValues.resize(size);
	    basisDerivs_u.resize(size);
	    basisDerivs_v.resize(size);
	    basisDerivs_uu.resize(size);
	    basisDerivs_uv.resize(size);
	    basisDerivs_vv.resize(size);
	}
};

/// Structure for storage of results of grid evaluation of the basis function of a spline surface.
/// Position, first, second and third derivatives
struct BasisDerivsSf3
{
  /// Parameter double in which the basis functions are evaluated
  double param[2];
  /// Index of the knot interval where the parameter value is situated for all
  /// parameter directions. The indices of the non-zero basis functions are
  /// left_idx[i]-order[i]+1, ..., left_idx[i] for i=0,1
  int left_idx[2];
  /// The value of all basis functions, size equal to (degree_u+1)*(degree_v+1)
  std::vector< double > basisValues;

  /// the derivative of all basis functions in u direction, same size as previous
  std::vector< double > basisDerivs_u;
  /// the derivative of all basis functions in v direction, same size as previous
  std::vector< double > basisDerivs_v;

  /// the second derivative of all basis functions twice in u direction, same size as previous
  std::vector< double > basisDerivs_uu;
  /// the second derivative of all basis functions in u and v direction, same size as previous
  std::vector< double > basisDerivs_uv;
  /// the second derivative of all basis functions twice in v direction, same size as previous
    std::vector< double > basisDerivs_vv;

    std::vector< double > basisDerivs_uuu;
    std::vector< double > basisDerivs_uuv;
    std::vector< double > basisDerivs_uvv;
    std::vector< double > basisDerivs_vvv;

    void prepareDerivs(double u, double v, int idx_u, int idx_v, int size)
    {
      param[0] = u;
      param[1] = v;
      left_idx[0] = idx_u;
      left_idx[1] = idx_v;
      basisValues.resize(size);
      basisDerivs_u.resize(size);
      basisDerivs_v.resize(size);
      basisDerivs_uu.resize(size);
      basisDerivs_uv.resize(size);
      basisDerivs_vv.resize(size);
      basisDerivs_uuu.resize(size);
      basisDerivs_uuv.resize(size);
      basisDerivs_uvv.resize(size);
      basisDerivs_vvv.resize(size);
    }
};

/// Structure for storage of results of grid evaluation of the basis function of a spline surface.
/// Position, and a given number of uni-directed derivatives
struct BasisDerivsSfU
{
  size_t derivs;
  /// Parameter double in which the basis functions are evaluated
  double param[2];
  /// Index of the knot interval where the parameter value is situated for all
  /// parameter directions. The indices of the non-zero basis functions are
  /// left_idx[i]-order[i]+1, ..., left_idx[i] for i=0,1
  int left_idx[2];

  /// The value of all basis functions and derivatives
  std::vector< std::vector<double> > values;

  /// Resize data structures
  /// \param u u parameter
  /// \param u v parameter
  /// \param idx_u u index
  /// \param idx_u v index
  /// \param deriv Number of derivatives
  /// \param size Number of functions
  void prepareDerivs(double u, double v, int idx_u, int idx_v,
		      size_t deriv, size_t size)
  {
    derivs = deriv;
    param[0] = u;
    param[1] = v;
    left_idx[0] = idx_u;
    left_idx[1] = idx_v;
    values.resize(size);
    for (size_t i=0;i<size;++i)
      values[i].resize(2*derivs+1);
  }
};

/// \brief SplineSurface provides methodes for storing,
/// reading and manipulating rational and non-rational
/// B-spline surfaces.
///
/// Non-rational B-spline surfaces represented on the form
/// \f[ \sum_{i=1}^{n_1}\sum_{j=1}^{n_2} P_{i,j} B_{i,o_1}(u) B_{j,o_2}(v), \f]
/// where the B-spline coefficients are stored in a vector
/// called coefs. The coefs are stored as
/// \f[ P_{0,0}, P_{1,0},  \ldots , P_{n_1, n_2}, \f]
/// where \f$ P_{i,j} \f$ is represented as dim doubles.
///
/// NURBS surfaces are represented on the form
/// \f[ \frac{\sum_{i=1}^{n_1}\sum_{j=1}^{n_2} w_{i,j}P_{i,j} B_{i,o_1}(u) B_{j,o_2}(v)} {\sum_{i=1}^{n_1}\sum_{j=1}^{n_2} w_{i,j} B_{i,o_1}(u) B_{j,o_2}(v)} \f]
/// where \f$ B_{i,o} \f$ is the <em> i</em>'th non-rational B-spline of
/// order \e o. The Projected coefficients are stored in the
/// coefs vector, i.e.
/// \f[ w_{0,0} \cdot P_{0,0}, w_{1,0} \cdot P_{1,0}, \ldots , w_{n_1, n_2} \cdot P_{n_1, n_2}. \f]
/// In addition the cefficients for the surface in projective
/// space are kept,  i.e
/// \f[ P_{0,0}, w_{0,0},  P_{1,0}, w_{1,0}, \ldots , P_{n_1, n_2}, w_{n_1, n_2}. \f]
/// With this representation surface is in the convex hull
/// of the coefficients stored in coefs both for rational
/// and non rational surfaces.
class GO_API SplineSurface : public ParamSurface
{
 public:
    /// Creates an uninitialized SplineSurface, which can only be assigned to 
    /// or read(...) into.
    SplineSurface()
      : ParamSurface(), dim_(-1), rational_(false), is_elementary_surface_(false)
    {
    }

    /// Create a SplineSurface by explicitly providing all spline-related 
    /// information.
    /// \param number1 number of control points along the u-parameter
    /// \param number2 number of control points along the v-parameter
    /// \param order1 BsplineBasis order along the u-parameter
    /// \param order2 BsplineBasis order along the v-parameter
    /// \param knot1start pointer to the array describing the knotvector 
    ///                   for the u-parameter
    /// \param knot2start pointer to the array describing the knotvector
    ///                   for the v-parameter
    /// \param coefsstart pointer to the array where the control points
    ///                   are consecutively stored.  The storage order is
    ///                   such that control points along the u-parameter have
    ///                   the shortest stride (stored right after each other).
    ///                   If the surface is rational, pay attention to the
    ///                   comments below.
    /// \param dim        dimension of the space in which the surface lies 
    ///                   (usually 3).  
    /// \param rational Specify whether the surface is rational or not.
    ///                 If the surface is rational, coefficients must be in
    ///                 the following format: 
    ///                 wP1 wP2 .... wPdim w.   Ie. a (dim+1)-dimensional form.
    ///                 (This is the same form that is used within SISL).
    template <typename RandomIterator1,
	      typename RandomIterator2,
	      typename RandomIterator3>
    SplineSurface(int number1,
		    int number2,
		    int order1,
		    int order2,
		    RandomIterator1 knot1start,
		    RandomIterator2 knot2start,
		    RandomIterator3 coefsstart,
		    int dim,
		    bool rational = false)
	: ParamSurface(), dim_(dim), rational_(rational),
        basis_u_(number1, order1, knot1start),
        basis_v_(number2, order2, knot2start), 
        is_elementary_surface_(false)
    {
	if (rational) {
	    int n = (dim+1)*number1*number2;
	    rcoefs_.resize(n);
	    std::copy(coefsstart, coefsstart + n, rcoefs_.begin());
	    coefs_.resize(dim*number1*number2);
	    updateCoefsFromRcoefs();
	} else {
	    int n = dim*number1*number2;
	    coefs_.resize(n);
	    std::copy(coefsstart, coefsstart + n, coefs_.begin());
	}
    }

    /// Create a SplineSurface by explicitly providing all spline-related 
    /// information.
    /// \param basis_u the BsplineBasis to be used along the u-direction
    /// \param basis_v the BsplineBasis to be used along the v-direction
    /// \param coefsstart pointer to the array where the control points
    ///                   are consecutively stored.  The storage order is
    ///                   such that control points along the u-parameter have
    ///                   the shortest stride (stored right after each other).
    ///                   If the surface is rational, pay attention to the
    ///                   comments below.
    /// \param dim dimension of the space in which the surface lies (usually 3).
    /// \param rational Specify whether the surface is rational or not.
    ///                 If the surface is rational, coefficients must be in the
    ///                 following format:
    ///                 wP1 wP2 .... wPdim w.   Ie. a (dim+1)-dimensional form.
    ///                 (This is the same form that is used within SISL).
    template <typename RandomIterator>
    SplineSurface(const BsplineBasis& basis_u,
		    const BsplineBasis& basis_v,
		    RandomIterator coefsstart,
		    int dim,
		    bool rational = false)
	: ParamSurface(), dim_(dim), rational_(rational),
        basis_u_(basis_u),
        basis_v_(basis_v),
        is_elementary_surface_(false)
    {
	int number1 = basis_u.numCoefs();
	int number2 = basis_v.numCoefs();
	if (rational) {
	    int n = (dim+1)*number1*number2;
	    rcoefs_.resize(n);
	    std::copy(coefsstart, coefsstart + n, rcoefs_.begin());
	    coefs_.resize(dim*number1*number2);
	    updateCoefsFromRcoefs();
	} else {
	    int n = dim*number1*number2;
	    coefs_.resize(n);
	    std::copy(coefsstart, coefsstart + n, coefs_.begin());
	}
    }

    /// Virtual destructor, enables safe inheritance.
    virtual ~SplineSurface();

    // inherited from Streamable
    virtual void read (std::istream& is);

    // inherited from Streamable
    virtual void write (std::ostream& os) const;


    // inherited from GeomObject
    virtual BoundingBox boundingBox() const;

    // inherited from GeomObject
    virtual int dimension() const;
    
    /// quick swap of two SplineSurface objects with each other
    void swap(SplineSurface& other);

    // inherited from GeomObject
    virtual ClassType instanceType() const;

    // inherited from GeomObject
    static ClassType classType()
    { return Class_SplineSurface; }
// #if ((_MSC_VER > 0) && (_MSC_VER < 1300))
//     virtual GeomObject* clone() const
//     { return new SplineSurface(*this); }
// #else
    virtual SplineSurface* clone() const;
// #endif

    /// Return a copy of the spline surface represented by this surface
    virtual SplineSurface* asSplineSurface() 
    {
        // We return a copy of this object, to avoid differing memory handling depending on surface type.
        return clone();
    }

    /// Return this surface
    virtual SplineSurface* getSplineSurface() 
    {
      return this;
    }

    // inherited from ParamSurface
// #if ((_MSC_VER > 0) && (_MSC_VER < 1300))
//     virtual const Domain& parameterDomain() const;
// #else
    virtual const RectDomain& parameterDomain() const;
// #endif

    // inherited from ParamSurface
    virtual RectDomain containingDomain() const;

    /// Check if a parameter pair lies inside the domain of this surface
    virtual bool inDomain(double u, double v, double eps=1.0e-4) const;

    /// Check if a parameter pair lies inside the domain of this surface
    /// return value = 0: outside
    ///              = 1: internal
    ///              = 2: at the boundary
    virtual int inDomain2(double u, double v, double eps=1.0e-4) const;

    /// Check if a parameter pair lies at the boundary of this surface
    virtual bool onBoundary(double u, double v, double eps=1.0e-4) const;

    virtual Point closestInDomain(double u, double v) const;

    // inherited from ParamSurface
    virtual CurveLoop outerBoundaryLoop(double degenerate_epsilon
					  = DEFAULT_SPACE_EPSILON) const;
    // inherited from ParamSurface
    virtual std::vector<CurveLoop> allBoundaryLoops(double degenerate_epsilon
						      = DEFAULT_SPACE_EPSILON) const;

    // inherited from ParamSurface
    virtual void point(Point& pt, double upar, double vpar) const;

    // Output: Partial derivatives up to order derivs (pts[0]=S(u,v),
    // pts[1]=dS/du=S_u, pts[2]=S_v, pts[3]=S_uu, pts[4]=S_uv, pts[5]=S_vv, ...)
    // inherited from ParamSurface
    virtual void point(std::vector<Point>& pts, 
		       double upar, double vpar,
		       int derivs,
		       bool u_from_right = true,
		       bool v_from_right = true,
		       double resolution = 1.0e-12) const;

    /// Get the start value for the u-parameter
    /// \return the start value for the u-parameter
    virtual double startparam_u() const;

    /// Get the start value for the v-parameter
    /// \return the start value for the v-parameter
    virtual double startparam_v() const;

    /// Get the start parameter value for the parameter direction idx
    double startparam(int idx) const
    {
      if (idx == 0)
	return startparam_u();
      else
	return startparam_v();
    }

    /// Get the end value for the u-parameter
    /// \return the end value for the u-parameter
    virtual double endparam_u() const;

    /// Get the end value for the v-parameter
    /// \return the end value for the v-parameter
    virtual double endparam_v() const;

    /// Get the end parameter value for the parameter direction idx
    double endparam(int idx) const
    {
      if (idx == 0)
	return endparam_u();
      else
	return endparam_v();
    }

     // inherited from ParamSurface
    virtual void normal(Point& n, double upar, double vpar) const;

    /// Enumerates the method for computing the normal cone
    enum NormalConeMethod { 
	SederbergMeyers = 0,
	SMCornersFirst = 1,
	sislBased = 2
    };
    /// Returns a cone that contains the convex hull of all normalized
    /// tangents of the surface. Note: It is an overestimate.
    /// \return the normal cone of the surface
    DirectionCone normalCone(NormalConeMethod method) const;

    /// Function that calls normalCone(NormalConeMethod) with method =
    /// sislbased. Needed because normalCone() is virtual! 
    /// (Inherited from ParamSurface).
    /// \return a DirectionCone (not necessarily the smallest) containing all normals 
    ///         to this surface.
    virtual DirectionCone normalCone() const;

    // inherited from ParamSurface
    virtual DirectionCone tangentCone(bool pardir_is_u) const;

    // Creates a composite box enclosing the surface. The composite box
    // consists of an inner and an edge box. The inner box is
    // made from the interior control points of the surface, while the
    // edge box is made from the boundary control points.
    // Inherited from ParamSurface
    virtual CompositeBox compositeBox() const;

    /// Not yet implemented
    SplineSurface* normal() const;

    /// Returns the normal surface corresponding to this surface, as 
    /// described in: \n
    /// <tt> Computing normal vector Bezier patches, Przemyslaw Kiciak,</tt>
    /// <tt> Computer Aided Design 18 (2001), 699-710</tt>.
    /// \return pointer to a newly created SplineSurface which is the 
    ///         normal surface.  User assumes ownership of this object,
    ///         and is responsible for destroying it.
    SplineSurface* normalSurface() const;

    /// Get the derivative surface.
    /// Expresses the i,j-th derivative of a spline surface as
    /// a spline surface. Ported from sisl routine s1386.
    /// \param ider1 what derivative to use along the u parameter
    /// \param ider2 what derivative to use along the v parameter
    /// \return pointer to a newly constructed SplineSurface which is 
    ///         the derivative surface.  User assumes ownership of this
    ///         object, and is responsible for destroying it.
    SplineSurface* derivSurface(int ider1, int ider2) const;

    /// Get a SplineSurface which represent a part of 'this' SplineSurface
    /// \param from_upar start value for u-parameter in the sub-surface to be 
    ///                  generated 
    /// \param from_vpar start value for v-parameter in the sub-surface to be
    ///                  generated
    /// \param to_upar end value for u-parameter in the sub-surface to be generated
    /// \param to_vpar end value for v-parameter in the sub-surface to be generated
    /// \param fuzzy tolerance used to determine whether given parameter values
    ///              are located directly \em knot values.
    SplineSurface* subSurface(double from_upar,
			      double from_vpar,
			      double to_upar,
			      double to_vpar,
			      double fuzzy =
			      DEFAULT_PARAMETER_EPSILON) const;

    // inherited from ParamSurface
    virtual std::vector<shared_ptr<ParamSurface> >
    subSurfaces(double from_upar, double from_vpar,
		double to_upar, double to_vpar,
		double fuzzy = DEFAULT_PARAMETER_EPSILON) const;

    /// Mirror a surface around a specified plane
    virtual SplineSurface* mirrorSurface(const Point& pos, const Point& norm) const;

    // inherited from ParamSurface
    virtual void closestPoint(const Point& pt,
			      double&        clo_u,
			      double&        clo_v, 
			      Point&         clo_pt,
			      double&        clo_dist,
			      double         epsilon,
			      const RectDomain* domain_of_interest = NULL,
			      double   *seed = 0) const;

    // inherited from ParamSurface
    virtual void closestBoundaryPoint(const Point& pt,
				      double&        clo_u,
				      double&        clo_v, 
				      Point&         clo_pt,
				      double&        clo_dist,
				      double         epsilon,
				      const RectDomain* rd = NULL,
				      double *seed = 0) const;

    /// IF POSSIBLE create the surface defined by appending the specified
    /// suface to 'this' surface along a specified parameter direction.
    /// Requires consistent surface types and surface types that support
    /// append
    virtual shared_ptr<ParamSurface> getAppendSurface(ParamSurface *sf,
						     int join_dir,
						     int cont,
						     double& dist,
						     bool repar=true);

    /// Appends a surface to the end of 'this' surface along a specified parameter
    /// direction.  The knotvector and control point grid of 'this' surface will be
    /// extended, and the degree will be raised if necessary.  Note that there \em
    /// will be side-effects on the other surface; its order might be raised, its
    /// knotvector will become k-regular and reparametrized, one of its edges will be
    /// moved to coincide with an edge of 'this' surface, etc.
    /// \param sf the surface to append to 'this' surface
    /// \param join_dir the parameter direction along which to extend 'this' surface
    ///                 by appending the 'sf' surface.  If it is equal to 1, we will
    ///                 extend the surface along the u-parameter; if it equals 2, we
    ///                 will extend the surface along the v-parameter.
    /// \param continuity the desired level of continuity at the transition between
    ///                   the two surfaces (can be from -1 to order()).  The higher 
    ///                   the value, the more the surfaces will have to be locally 
    ///                   modified around the seam.
    /// \param dist upon function return, this variable will hold the estimated
    ///             maximum distorsion after the 'smoothing' of the seam in order 
    ///             to achieve the desired continuity.   
    /// \param repar The parametrization of the concerned parameter in the other 
    ///              surface will \em always be shifted so that it starts where the
    ///              parametrization of 'this' surface ends.  However, if 'repar' is 
    ///              set to 'true', it will also be \em scaled as a function of position
    ///              of control points close to the transition.
    /// \param return Set in the rational case. Maximum adjustment of input weight
    double appendSurface(ParamSurface* sf, int join_dir,
		       int continuity, double& dist, bool repar=true);

    //void appendSurface(SplineSurface* sf, int join_dir, int continuity, double& dist);

    /// Short hand function to call \ref appendSurface with C^1 continuity.
    /// \param sf the surface to append to 'this' surface
    /// \param join_dir the parameter direction along which to extend 'this' surface
    ///                 by appending the 'sf' surface.  If it is equal to 1, we will
    ///                 extend the surface along the u-parameter; if it equals 2, we
    ///                 will extend the surface along the v-parameter.
    /// \param repar The parametrization of the concerned parameter in the other 
    ///              surface will \em always be shifted so that it starts where the
    ///              parametrization of 'this' surface ends.  However, if 'repar' is 
    ///              set to 'true', it will also be \em scaled as a function of position
    ///              of control points close to the transition.
    void appendSurface(ParamSurface* sf, int join_dir, bool repar=true);

    // inherited from ParamSurface
    virtual void getBoundaryInfo(Point& pt1, Point& pt2, 
				 double epsilon, SplineCurve*& cv,
				 SplineCurve*& crosscv, double knot_tol = 1e-05) const;

    /// Get the boundary curve and the outward cross tangent curve between
    /// two parameter values on a given boundary.
    /// \param par1 first parameter along the boundary
    /// \param par2 second parameter along the boundary
    /// \param bdindex index of the boundary in question.  The boundaries are
    ///                indexed as in \ref getBoundaryIdx().
    /// \param cv returns a pointer to a generated SplineCurve representing the
    ///           boundary between the given parameters.
    ///           The user assumes ownership and is responsible for deletion.
    /// \param crosscv returns a pointer to a generated SplineCurve representing
    ///                the cross tangent curve between the given parameters on the
    ///                boundary. The direction is outwards from the surface.
    ///                The user assumes ownership and is responsible for
    ///                deletion.
    /// \param knot_tol tolerance used when working with the knot-vector, to specify how
    ///                 close a parameter value must be to a knot in order to be considered
    ///                 'on' the knot.
    void getBoundaryInfo(double par1, double par2,
			 int bdindex, SplineCurve*& cv,
			 SplineCurve*& crosscv, double knot_tol = 1e-05) const;

    /// Given two points on the surface boundary, find the number of the
    /// corresponding boundary and the curve parameter of the closest points
    /// on this surface boundary.
    /// 
    /// Ordering of boundaries: 
    /// \verbatim
    ///                      1 
    ///           ---------------------- 
    ///           |                    |
    ///         2 |                    | 3
    ///      v    |                    |
    ///      ^    ----------------------
    ///      |-> u            0
    /// \endverbatim
    /// \param pt1 point in geometry space
    /// \param pt2 point in geometry space
    /// \param epsilon geometric tolerance
    /// \param bdindex if the two points are on a common boundary, 
    ///                the index of the boundary, otherwise -1.
    /// \param par1 if bdindex is 0 or 1, the u parameter of pt1, 
    ///             otherwise the v parameter.
    /// \param par2 if bdindex is 0 or 1, the u parameter of pt2, 
    ///             otherwise the v parameter.
    /// \param knot_tol tolerance used when working with the knot-vector, to specify how
    ///                 close a parameter value must be to a knot in order to be considered
    ///                 'on' the knot.
    void getBoundaryIdx(Point& pt1, Point& pt2, 
			double epsilon, int& bdindex,
			double& par1, double& par2, double knot_tol = 1e-05) const;

    /// Given one point on the surface boundary, find the index number of the
    // corresponding boundary.  The Boundaries are indexed as in \getBoundaryIdx().
    /// \param pt1 the point we want to find the boundary for
    /// \param epsilon the geometric distance the point can have from a boundary and 
    ///                still be considered as laying on the boundary.
    /// \param bdindex the index of the boundary on which the point lies.  If the 
    ///                point was found not to lie on any boundary, the value will 
    ///                be -1.
    /// \param knot_tol tolerance used when working with the knot-vector, to specify how
    ///                 close a parameter value must be to a knot in order to be considered
    ///                 'on' the knot.
    void getBoundaryIdx(Point& pt1, double epsilon, int& bdindex,
			double knot_tol = 1e-05) const;

    // inherited from ParamSurface
    virtual void turnOrientation();

    // inherited from ParamSurface
    virtual void swapParameterDirection();

    // inherited from ParamSurface
    virtual void reverseParameterDirection(bool direction_is_u);

    /** If the two points pt1 and pt2 are on the same edge the
     edge index is returned, otherwise -1 is returned.
     Edges are counted anticlockwise, and starts with the
     edge (upar, 0). */

    /// find whether two points are lying on the same surface boundary; if
    /// this is the case, return the index of the boundary (as specified
    /// in \ref getBoundaryIdx()), else return 1
    /// \return The boundary index on which the points lie, if they are found
    ///         to both lie on the same boundary.  In the opposite case, -1 
    ///         is returned.
    int boundaryIndex(Point& param_pt1, Point& param_pt2) const;

    /// Compute the total area of this surface up to some tolerance
    /// \param tol the relative tolerance when approximating the area, i.e.
    ///            stop iteration when error becomes smaller than
    ///            tol/(surface area)
    /// \return the area calculated
    virtual double area(double tol) const;

    /// get a const reference to the BsplineBasis for the first parameter
    /// \return const reference to the BsplineBasis for the first parameter
    const BsplineBasis& basis_u() const
    { return basis_u_; }

    /// get a const reference to the BsplineBasis for the second parameter
    /// \return const reference to the BsplineBasis for the second parameter
    const BsplineBasis& basis_v() const
    { return basis_v_; }

    /// get a reference to the BsplineBasis for the first parameter
    /// \return reference to the BsplineBasis for the first parameter
    BsplineBasis& basis_u()
    { return basis_u_; }

    /// get a reference to the BsplineBasis for the second parameter
    /// \return reference to the BsplineBasis for the second parameter
    BsplineBasis& basis_v()
    { return basis_v_; }

    /// get one of the BsplineBasises of the surface
    /// \param i specify whether to return the BsplineBasis for the first 
    ///          parameter (0), or for the second parameter (1).
    /// \return const reference to the requested BsplineBasis.
    const BsplineBasis& basis(int i) const
    { return (i==0) ? basis_u_ : basis_v_; }

    /// Query the number of control points along the first parameter direction
    /// \return the number of control points along the first parameter direction
    int numCoefs_u() const
    { return basis_u_.numCoefs(); }

    /// Query the number of control points along the second parameter direction
    /// \return the number of control points along the second parameter direction
    int numCoefs_v() const
    { return basis_v_.numCoefs(); }

    /// Query the number of elements in the SplineSurface
    int numElem() const
    {
      return basis_u_.numElem()*basis_v_.numElem();
    }

    /// Query the number of elements in one parameter direction of 
    /// the SplineSurface
    int numElem(int pardir) const
    {
      return (pardir == 0) ? basis_u_.numElem() :
	basis_v_.numElem();
    }

    /// Query the order of the BsplineBasis for the first parameter
    /// \return  the order of the BsplineBasis for the first parameter
    int order_u() const
    { return basis_u_.order(); }

    /// Query the order of the BsplineBasis for the second parameter
    /// \return the order of the BsplineBasis for the second parameter
    int order_v() const
    { return basis_v_.order(); }

    /// Query whether the surface is rational
    /// \return 'true' if the surface is rational, 'false' otherwise
    bool rational() const
    { return rational_; }
    
    /// Get an iterator to the start of the internal array of non-rational control 
    /// points.
    /// \return an (nonconst) iterator to the start of the internal array of non-
    ///         rational control points
    std::vector<double>::iterator coefs_begin()
    { return coefs_.begin(); }

    /// Get an iterator to the one-past-end position of the internal array of non-
    /// rational control points
    /// \return an (nonconst) iterator to the one-past-end position of the internal
    ///         array of non-rational control points
    std::vector<double>::iterator coefs_end()
    { return coefs_.end(); }

    /// Get a const iterator to the start of the internal array of non-rational
    /// control points.
    /// \return a const iterator to the start of the internal array of non-rational
    ///         control points.
    std::vector<double>::const_iterator coefs_begin() const
    { return coefs_.begin(); }

    /// Get a const iterator to the one-past-end position of the internal array of
    /// non-rational control points.
    /// \return a const iterator to the one-past-end position of the internal array of
    ///         non-rational control points.
    std::vector<double>::const_iterator coefs_end() const
    { return coefs_.end(); }

    /// Get an iterator to the start of the internal array of \em rational control 
    /// points.
    /// \return an (nonconst) iterator ro the start of the internal array of rational
    ///         control points.
    std::vector<double>::iterator rcoefs_begin()
    { return rcoefs_.begin(); }

    /// Get an iterator to the one-past-end position of the internal array of 
    /// \em rational control points.
    /// \return an (nonconst) iterator to the start of the internal array of rational
    ///         control points.
    std::vector<double>::iterator rcoefs_end()
    { return rcoefs_.end(); }

    /// Get a const iterator to the start of the internal array of \em rational
    /// control points.
    /// \return a const iterator to the start of the internal array of rational
    ///         control points
    std::vector<double>::const_iterator rcoefs_begin() const
    { return rcoefs_.begin(); }

    /// Get a const iterator to the one-past-end position of the internal array
    /// of \em rational control points.
    /// \return a const iterator to the one-past-end position of the internal array 
    ///         of rational control points.
    std::vector<double>::const_iterator rcoefs_end() const
    { return rcoefs_.end(); }

    /// Get an iterator to the start of the internal array of active control 
    /// points.
    /// \return an (nonconst) iterator to the start of the internal array of 
    ///         rational or non-rational control points
    std::vector<double>::iterator ctrl_begin()
    { return rational_ ? rcoefs_.begin() : coefs_.begin(); }

    /// Get an iterator to the one-past-end position of the internal array of 
    /// active control points
    /// \return an (nonconst) iterator to the one-past-end position of the internal
    ///         array of rational or non-rational control points
    std::vector<double>::iterator ctrl_end()
    { return rational_ ? rcoefs_.end() : coefs_.end(); }

    /// Get a const iterator to the start of the internal array of active
    /// control points.
    /// \return a const iterator to the start of the internal array of rational or 
    ///         non-rational control points.
    std::vector<double>::const_iterator ctrl_begin() const
    { return rational_ ? rcoefs_.begin() : coefs_.begin(); }

    /// Get a const iterator to the one-past-end position of the internal array of
    /// active control points.
    /// \return a const iterator to the one-past-end position of the internal array of
    ///         rational or non-rational control points.
    std::vector<double>::const_iterator ctrl_end() const
    { return rational_ ? rcoefs_.end() : coefs_.end(); }

    /// Replace one specified coefficient (local enumeration)
    void replaceCoefficient(int ix, Point coef);

    /// Return all weights corresponding to this surface, a non-rational volume
    /// has all weights equal to one
    void getWeights(std::vector<double>& weights) const;

    // inherited from ParamSurface
    virtual bool isDegenerate(bool& b, bool& r,
		      bool& t, bool& l, double tolerance) const;

    /// Check for paralell and anti paralell partial derivatives in surface corners
    virtual void getDegenerateCorners(std::vector<Point>& deg_corners, double tol) const;

    /// Return surface corners, geometric and parametric points
    /// in that sequence
    virtual void 
      getCornerPoints(std::vector<std::pair<Point,Point> >& corners) const;

    /// set the parameter domain to a given rectangle
    /// \param u1 new min. value of first parameter span
    /// \param u2 new max. value of first parameter span
    /// \param v1 new min. value of second parameter span
    /// \param v2 new max. value of second parameter span
    virtual void setParameterDomain(double u1, double u2, double v1, double v2);

    /// Insert a new knot in the knotvector of the first parameter
    /// \param apar the parameter value at which a new knot will be inserted
    void insertKnot_u(double apar);
    
    /// Insert new knots in the knotvector of the first parameter
    /// \param new_knots a vector containing the parameter values of the
    ///                  new knots to insert.
    void insertKnot_u(const std::vector<double>& new_knots);

    /// Insert a new knot in the knotvector of the second parameter
    /// \param apar the parameter value at which a new knot will be inserted
    void insertKnot_v(double apar);
    
    /// Insert new knots in the knotvector of the second parameter
    /// \param new_knots a vector containing the parameter values of the 
    ///                  new knots to insert.
    void insertKnot_v(const std::vector<double>& new_knots);

    /// Remove a knot from the knotvector of the first parameter.
    /// \param tpar the parameter value of the knot to be removed
    void removeKnot_u(double tpar);

    /// Remove a knot from the knotvector of the second parameter.
    /// \param tpar the parameter value of the knot to be removed.
    void removeKnot_v(double tpar);

    /// Inserts knots in the u knot vector, such that all knots
    /// have multiplicity order
    void makeBernsteinKnotsU();

    /// Inserts knots in the v knot vector, such that all knots
    /// have multiplicity order
    void makeBernsteinKnotsV();

    /// Ensure k-regularity of this surface in both parameter directions
    void makeSurfaceKRegular();

    /// Returns the number of knot intervals in u knot vector.
    /// \return the number of knot intervals in the knotvector for the first 
    ///         parameter
    int numberOfPatches_u() const;

    /// Returns the number of knot intervals in v knot vector.
    /// \return the number of knot intervals in the knotvector for the second
    ///         parameter
    int numberOfPatches_v() const;

     /// Returns the size of the knot interval (knot[iknot],knot[iknot+1])
    /// in the specified direction. An index outside the legal range will 
    /// result in a zero knot span
    double knotSpan(int pardir, int iknot) const;

   /// Raise the order of the spline surface as indicated by parameters.
    /// \param raise_u the order of the BsplineBasis associated with the first
    ///                parameter will be raised this many times.
    /// \param raise_v the order of the BsplineBasis associated with the second
    ///                parameter will be raised this many times.
    void raiseOrder(int raise_u, int raise_v);

    /// Generate and return a SplineCurve that represents a constant parameter 
    /// curve on the surface
    /// \param parameter value of the fixed parameter
    /// \param pardir_is_u 'true' if the first parameter is the running parameter,
    ///                    'false' otherwise.
    /// \return pointer to a newly constructed SplineCurve representing the 
    ///         specified constant parameter curve.  It is the user's reponsibility
    ///         to delete it when it is no longer needed.
    SplineCurve* constParamCurve(double parameter,
				 bool pardir_is_u) const;

    /// Generate and return two SplineCurves, representing a constant parameter 
    /// curve on the surface as well as it cross-tangent curve.
    /// \param parameter value of the fixed parameter
    /// \param pardir_is_u 'true' if the first parameter is the running parameter,
    ///                    'false' otherwise.
    /// \param cv upon function completion, 'cv' will point to a newly constructed
    ///           SplineCurve representing the specified constant parameter curve.
    ///           It is the user's responsibility to delete it when it is no longer
    ///           needed.
    /// \param crosscv upon function completion, 'crosscv' will point to a newly
    ///                constructed SplineCurve representing the cross-tangent curve
    ///                of the specified constant parameter curve.  It is the user's
    ///                reponsibility to delete it when it is no longer needed.
    void constParamCurve(double parameter, 
			 bool pardir_is_u, 
			 SplineCurve*& cv, 
			 SplineCurve*& crosscv) const;

    /// Generate SplineCurves that represents constant parameter 
    /// curves on the surface. Do this for both directions
    /// \param parameter params_u the parameter values in u-directions for curves
    ///        running along the v-direction
    /// \param parameter params_v the parameter values in v-directions for curves
    ///        running along the u-direction
    /// \param curves_u upon function return, holds constant parameter curves
    ///        for the u-direction parameters, the curves are running along
    ///        the v-direction
    /// \param curves_v upon function return, holds constant parameter curves
    ///        for the v-direction parameters, the curves are running along
    ///        the u-direction
    void getConstParamCurves(const std::vector<double>& params_u,
			     const std::vector<double>& params_v,
			     std::vector<shared_ptr<SplineCurve> >& curves_u,
			     std::vector<shared_ptr<SplineCurve> >& curves_v);

    // inherited from ParamSurface
    virtual std::vector< shared_ptr<ParamCurve> >
    constParamCurves(double parameter, bool pardir_is_u) const;

    /// Here edge_number means:
    ///  0 -> bottom edge,
    ///  1 -> right edge,
    ///  2 -> top edge and
    ///  3 -> left edge.
    /// The returned pointer is new'ed and the responsibility
    /// of the caller. You should either delete it, or
    /// manage it with a smart pointer.


    /// Generate and return a SplineCurve representing one of the surface's four
    /// edges.
    /// \param ccw_edge_number indicates which edge the user requests. 
    /// \verbatim 
    /// 0 -> bottom edge
    /// 1 -> right edge
    /// 2 -> top edge
    /// 3 -> left edge \endverbatim
    /// \return A pointer to a newly constructed SplineCurve representing the
    ///         requested edge.  It is the user's responsibility to delete it when
    ///         it is no longer needed.
  // @@sbr What about the orientation? Seems it follows that of the sf.
    SplineCurve* edgeCurve(int ccw_edge_number) const;

    /// Remake 'this'surface to interpolate (or approximate)  a given point grid.
    /// \param interpolator1 Interpolator to apply to the first parameter direction
    /// \param interpolator2 Interpolator to apply to the second parameter direction
    /// \param num_points1 number of points to interpolate along the first parameter 
    ///                    direction
    /// \param num_points2 number of points to interpolate along the second parameter
    ///                    direction
    /// \param dim dimension of the points to interpolate (usually 3)
    /// \param param1_start pointer to the start of the input grid's parametrization
    ///                     along the first parameter
    /// \param param2_start pointer to the start of the input grid's parametrization
    ///                     along the second parameter
    /// \param data_start pointer to the start of the storage array for the data point
    ///                   grid to interpolate
    void interpolate(Interpolator& interpolator1,
		     Interpolator& interpolator2,
		     int num_points1,
		     int num_points2,
		     int dim,
		     const double* param1_start,
		     const double* param2_start,
		     const double* data_start);

    /// Evaluate points in a grid
    /// The nodata value is applicable for bounded surfaces
    /// and will not be used in this context
    virtual void evalGrid(int num_u, int num_v, 
			  double umin, double umax, 
			  double vmin, double vmax,
			  std::vector<double>& points,
			  double nodata_val = -9999) const;

    /// Evaluate points and normals on an entire grid, taking computational advantage
    /// over calculating all these values simultaneously rather than one-by-one.
    /// \param num_u number of values to evaluate along first parameter direction
    /// \param num_v number of values to evaluate along second parameter direction
    /// \param points upon function return, this vector holds all the evaluated points
    /// \param normals upon function return, this vector holds all the evaluated normals
    /// \param param_u upon function return, this vector holds all the numerical values
    ///                for the first parameter where evaluation has taken place
    /// \param param_v upon function return, this vector holds all the numerical values
    ///                for the second parameter where evaluation has taken place.
    /// \param normalize tells whether the normal vectors should be normalized
    void gridEvaluator(int num_u, int num_v,
		       std::vector<double>& points,
		       std::vector<double>& normals,
		       std::vector<double>& param_u,
		       std::vector<double>& param_v,
		       bool normalize = true) const;

    /// Evaluate points on an entire grid, taking computational advantage
    /// over calculating all these values simultaneously rather than one-by-one.
    /// \param num_u number of values to evaluate along first parameter direction
    /// \param num_v number of values to evaluate along second parameter direction
    /// \param points upon function return, this vector holds all the evaluated points
    /// \param param_u upon function return, this vector holds all the numerical values
    ///                for the first parameter where evaluation has taken place
    /// \param param_v upon function return, this vector holds all the numerical values
    ///                for the second parameter where evaluation has taken place.
    void gridEvaluator(int num_u, int num_v,
		       std::vector<double>& points,
		       std::vector<double>& param_u,
		       std::vector<double>& param_v) const
    {
      gridEvaluator(num_u, num_v,
		    points, param_u, param_v,
		    startparam_u(),
		    endparam_u(),
		    startparam_v(),
		    endparam_v());
    }

    /// Evaluate points on an entire grid, taking computational advantage
    /// over calculating all these values simultaneously rather than one-by-one.
    /// \param num_u number of values to evaluate along first parameter direction
    /// \param num_v number of values to evaluate along second parameter direction
    /// \param points upon function return, this vector holds all the evaluated points
    /// \param param_u upon function return, this vector holds all the numerical values
    ///                for the first parameter where evaluation has taken place
    /// \param param_v upon function return, this vector holds all the numerical values
    ///                for the second parameter where evaluation has taken place.
    /// \param start_u start value of first parameter range
    /// \param end_u end value of first parameter range
    /// \param start_v start value of second parameter range
    /// \param end_v end value of second parameter range
    void gridEvaluator(int num_u, int num_v,
		       std::vector<double>& points,
		       std::vector<double>& param_u,
		       std::vector<double>& param_v,
		       double start_u,
		       double end_u,
		       double start_v,
		       double end_v) const;

    /// Evaluate points on an entire grid, taking computational advantage
    /// over calculating all these values simultaneously rather than one-by-one.
    /// \param points upon function return, this vector holds all the evaluated points
    /// \param param_u this vector holds all the numerical values
    ///                for the first parameter where evaluation should take place
    /// \param param_v this vector holds all the numerical values
    ///                for the second parameter where evaluation should take place
    void gridEvaluator(std::vector<double>& points,
		       const std::vector<double>& param_u,
		       const std::vector<double>& param_v) const;

    /// Evaluate points and derivatives on an entire grid, taking computational advantage
    /// over calculating all these values simultaneously rather than one-by-one.
    /// \param params_u the values for the first parameter where evaluation takes place
    /// \param params_v the values for the second parameter where evaluation takes place
    /// \param points upon function return, this vector holds all the evaluated points
    /// \param derivs_u upon function return, this vector holds all derivations w.r.t. u
    /// \param derivs_u upon function return, this vector holds all derivations w.r.t. v
    /// \param evaluate_from_right specifies directional derivatives, true=right, false=left, 
    ///                            same behaviour for both parameter directions
    void gridEvaluator(const std::vector<double>& params_u,
		       const std::vector<double>& params_v,
		       std::vector<double>& points,
		       std::vector<double>& derivs_u,
		       std::vector<double>& derivs_v,
		       bool evaluate_from_right = true) const;

    /// Evaluate positions and first derivatives of all basis values in a given parameter pair
    /// For non-rationals this is an interface to BsplineBasis::computeBasisValues 
    /// where the basis values in each parameter direction are multiplied to 
    /// compute B_i(u)*B_j(v), for rationals the routine evaluates the rational
    /// basis functions, i.e. the basis functions are divided by the denominator of the
    /// surface in the given parameter direction
    /// \param param the parameter pair in which to compute
    /// \param basisValues the value of all basis functions, size equal to 
    ///                    (degree_u+1)*(degree_v+1)
    /// \param basisDerivs the derivative of all basis functions, same size as previous
    /// \param evaluate_from_right specifies directional derivatives, true=right, false=left, 
    ///                            same behavious for both parameter directions
    void computeBasis(double param[], 
		      std::vector< double > &basisValues,
		      std::vector< double > &basisDerivs_u,
		      std::vector< double > &basisDerivs_v,
		      bool evaluate_from_right = true) const;

    /// Evaluate positions and first derivatives of all basis values, when the B-splines
    /// and their derivatives have all been precalculated and are given as input
    /// \param bas_vals_u the basis values and derivatives in first parameter direction
    /// \param bas_vals_v the basis values and derivatives in second parameter direction
    /// \param left_u index of first non-zero interval in first parameter direction
    /// \param left_v index of first non-zero interval in second parameter direction
    /// \param basisValues the value of all basis functions, size equal to 
    ///                    (degree_u+1)*(degree_v+1)
    /// \param basisDerivs the derivative of all basis functions, same size as previous
    void computeBasis(const std::vector<double>::const_iterator& bas_vals_u,
		      const std::vector<double>::const_iterator& bas_vals_v,
		      int left_u,
		      int left_v,
		      std::vector<double>& basisValues,
		      std::vector<double>& basisDerivs_u,
		      std::vector<double>& basisDerivs_v) const;

    /// Convenience to be used in computations of basis grids
    typedef std::vector<double>  Dvector; 
    /// Convenience to be used in computations of basis grids
    typedef std::vector<Dvector> Dmatrix; 

    /// Evaluate positions of all basis values in a specified
    /// grid.  For non-rationals this is an interface to BsplineBasis::computeBasisValues 
    /// where the basis values in each parameter direction are multiplied to 
    /// compute B_i(u)*B_j(v), for rationals the routine evaluates the rational 
    /// basis functions, i.e. the basis functions are divided by the denominator of the surface in 
    /// the given parameter direction
    /// \param param_u this vector holds all the numerical values for the first parameter 
    ///                where evaluation will tak place
    /// \param param_v this vector holds all the numerical values for the second parameter 
    ///                where evaluation will tak place
    /// \param basisValues the value of all basis functions. The first matrix dimension 
    ///                will correspond to all evaluation parameters, i.e. the size is 
    ///                param_u.size()*param_v*size(). The second 
    ///                dimension correspond to all coefficients, i.e. 
    ///                nmb_coef_u*nmb_coef_v. Most entries in the matrix
    ///                will be zero
    void computeBasisGrid(const Dvector& param_u,
			  const Dvector& param_v,
			  Dmatrix& basisValues) const; 

    /// Compute basis values (position) in the parameter (param_u,param_v).
    /// Store result in a BasisPtsSf entity
    void computeBasis(double param_u,
		      double param_v,
		      BasisPtsSf& result) const;

    /// Compute basis values (position and 1. derivatives) in the parameter 
    /// (param_u,param_v). Store result in a BasisDerivSf entity
    void computeBasis(double param_u,
		      double param_v,
		      BasisDerivsSf& result,
		      bool evaluate_from_right = true) const;

    /// Compute basis values (position and uni-directed derivatives) in the parameter
    /// (param_u,param_v). Store result in a BasisDerivsSfU entity
    void computeBasis(double param_u,
		      double param_v,
                      int derivs,
		      BasisDerivsSfU& result,
		      bool evaluate_from_right = true) const;

    /// Compute basis values (position and 1. and 2. derivatives) in the parameter 
    /// (param_u,param_v). Store result in a BasisDerivSf2 entity
     void computeBasis(double param_u,
		      double param_v,
		      BasisDerivsSf2& result,
		      bool evaluate_from_right = true) const;

    /// Compute basis values (position and 1., 2. and 3. derivatives) in the parameter
    /// (param_u,param_v). Store result in a BasisDerivSf3 entity
     void computeBasis(double param_u,
                       double param_v,
                       BasisDerivsSf3& result,
                       bool evaluate_from_right = true) const;

    /// Compute basis grid (position) in the parameter pairs combined from param_u
    /// and param_v. Store result in a vector of BasisPtsSf.
    void computeBasisGrid(const Dvector& param_u,
			  const Dvector& param_v,
			  std::vector<BasisPtsSf>& result) const; 

    /// Evaluate positions and first derivatives of all basis values in a specified
    /// grid.  For non-rationals this is an interface to BsplineBasis::computeBasisValues 
    /// where the basis values in each parameter direction are multiplied to 
    /// compute B_i(u)*B_j(v), for rationals the routine evaluates the rational 
    /// basis functions, i.e. the basis functions are divided by the denominator of the surface in 
    /// the given parameter direction
    /// \param param_u this vector holds all the numerical values for the first parameter 
    ///                where evaluation will tak place
    /// \param param_v this vector holds all the numerical values for the second parameter 
    ///                where evaluation will tak place
    /// \param basisValues the value of all basis functions. The first matrix dimension 
    ///                will correspond to all evaluation parameters, i.e. the size is 
    ///                param_u.size()*param_v*size(). The second 
    ///                dimension correspond to all coefficients, i.e. 
    ///                nmb_coef_u*nmb_coef_v. Most entries in the matrix
    ///                will be zero
    /// \param basisDerivs_u the derivative of all basis functions, size as above
    /// \param basisDerivs_v the derivative of all basis functions
    /// \param evaluate_from_right specifies directional derivatives, true=right, false=left
    void computeBasisGrid(const Dvector& param_u,
			  const Dvector& param_v,
			  Dmatrix& basisValues,
			  Dmatrix& basisDerivs_u,
			  Dmatrix& basisDerivs_v,
			  bool evaluate_from_right = true) const; 


    /// Compute basis grid (position and 1. derivatives) in the parameter pairs 
    /// combined from param_u and param_v. Store result in a vector of BasisDerivSf.
    void computeBasisGrid(const Dvector& param_u,
			  const Dvector& param_v,
			  std::vector<BasisDerivsSf>& result,
			  bool evaluate_from_right = true) const; 


    /// Compute basis grid (position and 1. and 2. derivatives) in the parameter pairs 
    /// combined from param_u and param_v. Store result in a vector of BasisDerivSf2.
    void computeBasisGrid(const Dvector& param_u,
			  const Dvector& param_v,
			  std::vector<BasisDerivsSf2>& result,
			  bool evaluate_from_right = true) const;

    /// Compute basis grid (position and 1., 2. and 3. derivatives) in the parameter pairs
    /// combined from param_u and param_v. Store result in a vector of BasisDerivSf3.
    void computeBasisGrid(const Dvector& param_u,
                          const Dvector& param_v,
                          std::vector<BasisDerivsSf3>& result,
                          bool evaluate_from_right = true) const;


        // inherited from ParamSurface
    virtual double nextSegmentVal(int dir, double par, bool forward, double tol) const;

    /// Replace one boundary curve of this surface
    /// Update spline spaces to enable replacement
    /// NB! Requires both the surface and the new boundary curve to be
    /// either rational or non-rational
    /// bd_nmb = 0: umin
    ///        = 1: umax
    ///        = 2: vmin
    ///        = 3: vmax
    /// Return value: true if a replacement is performe
    bool replaceBoundaryCurve(int bd_nmb, shared_ptr<SplineCurve> bd_crv,
			      bool unify=true);

    /// Check if the surface is of type spline
    virtual bool isSpline() const
    {
      return true;
    }

    /// Adds the given deformation vector to the coefficients.
    void deform(const std::vector<double>& vec, int vdim = 0);

    /// Add coefficients from another surface. Weights are not summed for rational cases
    /// Nothing is done and exception is raised if
    /// - Spline spaces are different in any paramter direction (order or knot vectors are not identical)
    /// - The geometry spaces of the surfaces have different dimension
    /// - One surface is rational while the other is not
    /// - The weights are not equal within a given tolerance (only if surfaces are rational)
    /// \param other The other spline surface with coefficients to be added into this surface
    /// \param tol tolerance used to test if weights are considered equal in rational case
    void add(const SplineSurface* other, double tol = 1.0e-10);

    /// Create function by multiplying the coefficients of the this 
    /// surface with a given vector. 
    SplineSurface* multCoefs(const Point& vec) const;

    /// Ensure that the current surface is represented as a rational surface
    void representAsRational();

    /// Set the average weight at one boundary to a given value (if rational)
    /// pardir = 0 : 1. parameter direction
    /// pardir = 1 : 2. parameter direction
    /// output : Variance between existing weights
    double setAvBdWeight(double wgt, int pardir, bool at_start);

    /// Check if the surface is axis rotational. Only true if a connection
    /// to an axis rotational elementary surface exist
    virtual bool isAxisRotational(Point& centre, Point& axis, Point& vec,
				  double& angle);

    /// This surface is planar if it represents a plane or the
    /// spline surface is linear in both parameter direction and planar
    virtual bool isPlanar(Point& normal, double tol);

   /// Check if a polynomial element (for spline surfaces) intersects the
    /// (trimming) boundaries of this surface
    /// \param elem_ix: Element index counted according to distinct knot
    /// values. Sequence of coordinates: x runs fastest, then y
    /// \param eps: Intersection tolerance
    /// \return -1: Not a spline surface or element index out of range
    ///          0: Not on boundary or touching a boundary curve
    ///          1: On boundary (intersection with boundary found)
    /// Note that a touch with the boundaries of the underlying surfaces
    /// is not consdered a boundary intersection while touching a trimming
    /// curve is seen as an intersection
    virtual int ElementOnBoundary(int elem_ix, double eps)
    {
      return 0;  
    }

   /// Check if a polynomial element (for spline surfaces) intersects the
    /// (trimming) boundaries of this ftSurface, is inside or outside
    /// \param elem_ix: Element index counted according to distinct knot
    /// values. Sequence of coordinates: x runs fastest, then y
    /// \param eps: Intersection tolerance
    /// \return -1: Not a spline surface or element index out of range
    ///          0: Outside trimmed volume
    ///          1: On boundary (intersection with boundary found)
    ///          2: Internal to trimmed surfaces
    /// Note that a touch with the boundaries of the underlying surface
    /// is not consdered a boundary intersection while touching a trimming
    /// curve is seen as an intersection
    virtual int ElementBoundaryStatus(int elem_ix, double eps)
    {
      return 2;
    }

    /// Fetch curves in the parameter domain surrounding the specified element
    /// elem_par - parameter values of element boundaries, sequence umin,
    /// umax, vmin, vmax
    std::vector<shared_ptr<SplineCurve> > getElementBdParCvs(int elem_ix,
							     double elem_par[]);

    /// Check if the surface has stored information about an original
    /// surface
    virtual bool hasParentSurface() const
    {
      return (elementary_surface_.get() != NULL);
    }
      
    /// Return an eventual original surface
    virtual shared_ptr<ParamSurface> getParentSurface();
      
    /// Query if the surface was generated from an ElementarySurface
    bool isElementarySurface() 
    {
      return is_elementary_surface_;
    }

    /// Get shared pointer to ElementarySurface, if it exists. If not, return
    /// empty shared pointer.
    ///
    /// NOTE: The ElementarySurface returned by this function should
    /// ideally be the one corresponing to the current
    /// SplineSurface. However, there is no guarantee for this. One
    /// may check to see if this is the case by using
    /// checkElementarySurface().
    shared_ptr<ElementarySurface> getElementarySurface();

    /// Return associated elementary surface, if any
    virtual ElementarySurface* elementarySurface()
    {
      return elementary_surface_.get();
    }
      
    /// Set shared pointer to the ElementarySurface that is
    /// represented by \c this.
    ///
    /// NOTE: The current SplineSurface (i.e. \c this), must be
    /// created from the argument ElementarySurface by
    /// ElementarySurface::createSplineSurface(), otherwise undefined
    /// behaviour may occur. One may check to see if this is the case
    /// by using checkElementarySurface().
    void setElementarySurface(shared_ptr<ElementarySurface> elsurf);

    /// Check to see if \c this corresponds to the ElementarySurface
    /// set by setElementarySurface().
    /// NOTE: Not yet implemented!
    bool checkElementarySurface();

    /// Linearly extend surface a given length along a given parameter direction
    /// (u if 'in_u' is true, v otherwise), before parameter start value
    /// ('at_end' = false) or after the parameter end value ('at_end' = true).
    void enlarge(double len, bool in_u, bool at_end);

    // Linearly extend surface a given length along each parameter direction,
    // before the parameter start value and after the parameter end value.
    void enlarge(double l_umin, double l_umax, double l_vmin, double l_vmax);
    
 private:

    // Canonical data
    int dim_;
    bool rational_;
    BsplineBasis basis_u_;
    BsplineBasis basis_v_;
    std::vector<double> coefs_;   // Like ecoef in SISL
    std::vector<double> rcoefs_;  // Like rcoef in SISL, only used if rational

    // Generated data
    mutable RectDomain domain_;
    mutable CurveLoop spatial_boundary_;

    // Data about origin or history
    bool is_elementary_surface_;
    shared_ptr<ElementarySurface> elementary_surface_;

    // Helper functions
    void updateCoefsFromRcoefs();
    std::vector<double>& activeCoefs() { return rational_ ? rcoefs_ : coefs_; }
    bool normal_not_failsafe(Point& n, double upar, double vpar) const;
    bool search_for_normal(bool interval_in_u,
			   double fixed_parameter,
			   double interval_start, // normal is not defined here
			   double interval_end,
			   Point& normal) const;

    // Members new to this class
 public:
    /// Evaluate points and possibly normals in a regular grid.
    /// The B-spline basis functions must be pre evaluated.
    /// This is a convenience funtion lifted to the user interface to allow
    /// for better performance if several spline surfaces with the same
    /// spline spaces are to be evaluated in the same points. Must be used with care.
    void pointsGrid(int m1, int m2, int derivs,
		    const double* basisvals1,
		    const double* basisvals2,
		    const int* knotint1,
		    const int* knotint2,
		    double* result,
		    double* normals = 0) const;
 private:

    // Rewritten pointsGrid, to avoid reformatting results.
    // Does not return the derivatives, works only for a 3D non-rational
    // spline.
    void pointsGridNoDerivs(int m1, int m2,
			    const double* basisvals1,
			    const double* basisvals2,
			    const int* knotint1,
			    const int* knotint2,
			    double* result,
			    double* normals = 0,
			    bool normalize = true) const;

    void accumulateBasis(const std::vector<double>::const_iterator& basisvals_u,
			 const std::vector<double>::const_iterator& basisvals_v,
			 const std::vector<double>& weights,
			 std::vector<double>& basisValues) const;

    void accumulateBasis(const std::vector<double>::const_iterator& basisvals_u,
			 const std::vector<double>::const_iterator& basisvals_v,
			 const std::vector<double>& weights,
			 std::vector<double>& basisValues,
			 std::vector<double>& basisDerivs_u,
			 std::vector<double>& basisDerivs_v) const;

    void accumulateBasis(const std::vector<double>::const_iterator& basisvals_u,
			 const std::vector<double>::const_iterator& basisvals_v,
			 const std::vector<double>& weights,
			 std::vector<double>& basisValues,
			 std::vector<double>& basisDerivs_u,
			 std::vector<double>& basisDerivs_v,
			 std::vector<double>& basisDerivs_uu,
			 std::vector<double>& basisDerivs_uv,
			 std::vector<double>& basisDerivs_vv) const;

    void accumulateBasis(const std::vector<double>::const_iterator& basisvals_u,
                         const std::vector<double>::const_iterator& basisvals_v,
                         const std::vector<double>& weights,
                         BasisDerivsSf3& result) const;

    void accumulateBasis(const std::vector<double>& basisvals_u,
			 const std::vector<double>& basisvals_v,
			 const std::vector<double>& weights,
                         BasisDerivsSfU& output) const;

    // Actually computes the closest point. Only difference is the explicit
    // robust_seedfind-parameter, which is always true in the virtual
    // closestPoint() function, when that function calls closestPointImpl().
    void closestPointImpl(const Point& pt,
			  double&        clo_u,
			  double&        clo_v, 
			  Point&         clo_pt,
			  double&        clo_dist,
			  double         epsilon,
			  const RectDomain* domain_of_interest = NULL,
			  bool robust_seedfind = true,
			  double   *seed = 0) const;

    // Compute the area of the subsurface determined by a subdomain of the parametric space
    // Uses cross product of diagonal vectors as an estimation of the area.
    double areaDiagonalCross(double tol,
			     double start_u, double end_u,
			     double start_v, double end_v) const;

    void s1773(const double ppoint[],double aepsge, double estart[],double eend[],double enext[],
	       double gpos[],int *jstat) const;

    void s1773_s9corr(double gd[],double acoef1,double acoef2,
		      double astart1,double aend1,double astart2,double aend2) const;

    void s1773_s9dir(double *cdist,double *cdiff1,double *cdiff2,
		     double PS[],const double *eval1,std::vector<Point> eval2,
		     double aepsge, int idim,int *jstat) const;


};


} // namespace Go




#endif // _SPLINESURFACE_H

