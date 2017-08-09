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

#ifndef _LRSPLINESURFACE_H
#define _LRSPLINESURFACE_H

#include <array>
#include <functional>
#include <set>
#include <map>
#include <vector>
#include <unordered_map>
#include <iostream> // @@ debug
#include <memory>

#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/geometry/LineCloud.h"
#include "GoTools/lrsplines2D/Mesh2D.h"
#include "GoTools/lrsplines2D/LRBSpline2D.h"
#include "GoTools/lrsplines2D/Element2D.h"

namespace Go
{
  // =============================================================================
  class LRSplineSurface : public ParamSurface
// =============================================================================
{
   public:
  //
  // SEARCH STRUCTURES
  //
  // Structure representing a refinement to carry out.  This structure is used as an
  // argument to the LRSplineSurface::refine() methods below.  It is particularly useful
  // when passing along a whole batch of refinements. 
  struct Refinement2D {
    double kval;      // value of the fixed parameter of the meshrectangle to insert 
    double start;     // start value of the meshrectangle's non-fixed parameter
    double end;       // end value of the meshrectangle's non-fixed parameter
    Direction2D d;    // direction of the meshrectangle (XFIXED or YFIXED) YCONSTANT & XCONSTANT?
    int multiplicity; // multiplicity of the meshrectangle 

    void setVal(double val, double st, double e, Direction2D dir, int mult)
    {
      kval = val;
      start = st;
      end = e;
      d = dir;
      multiplicity = mult;
    }
  };

  // 'BSKey' defines the key for storing/looking-up individual B-spline basis 
  // functions. 'u_min','v_min', 'u_max' and 'v_max' designate the support corners
  // in the parametric domain. 'u_mult' and 'v_mult' designate the knot multiplicities
  // at the lower left corner.  Taken together, these six values uniquely determine
  // a particular B-spline function of a given bi-degree in a given mesh.
  //
  // NB: Note the oddity of using a structure with 'double' values as a key to 
  // a STL map.  This works as long as the values are _never_ separately computed, 
  // but  _always_ copied from the _same_ source before insertion or lookup.  
  // This can safely be assumed here though, as the doubles will all be taken
  // (copied) directly from the mesh (which never does any arithmetic on knot 
  // values, only stores them).  The reason for choosing to store 'double' values
  // within the key rather than just indices to an external knotvector is that
  // a later mesh refinement will change indexing of the knot vectors, whereas 
  // keys in an STL map are not allowed to change (for good reason).  The parameter
  // values themselves do not change, though, which make them suitable for use
  // in the key.
  struct BSKey 
  {
    double u_min, v_min, u_max, v_max;
    int u_mult1, v_mult1, u_mult2, v_mult2;
    bool operator<(const BSKey& rhs) const; // needed for sorting when used in an STL map
  };

  // these maps could be redefined as hash tables later, as this is likely to improve
  // performance (at the expense of having to specify hash functions for these types of keys).
  typedef std::map<BSKey, std::unique_ptr<LRBSpline2D> > BSplineMap; // storage of basis functions


  struct double_pair_hash {
    size_t operator()(const std::pair<double, double>& dp) const {
      // Two hashes for independent variables, combined using XOR.  It is expected
      // that the resulting hash is probably as good as the input hashes.
      return std::hash<double>()(dp.first) ^ std::hash<double>()(dp.second);
    }
  };


  // Function for generating the key to use when storing B-spline function 'b'.  (This is an 
  // implementation detail that should not worry users).
  static BSKey generate_key(const LRBSpline2D& b, const Mesh2D& m);
  static BSKey generate_key(const LRBSpline2D& b);

  // The ElementMap will be used to keep track over which BasisFunctions are covering each
  // element.  An element is represented by its lower-left coordinates (doubles). 
  // Note: the same comment regarding the use of 'double' as a key in a map applies here 
  // (c.f. comment above on BSKey).
   struct ElemKey
  {
    double u_min, v_min;
    bool operator<(const ElemKey& rhs) const; // needed for sorting when used in an STL map
  };
    
  // these maps could be redefined as hash tables later, as this is likely to improve
  // performance (at the expense of having to specify hash functions for these types of keys).
  typedef std::map<ElemKey, std::unique_ptr<Element2D> > ElementMap; // storage of basis functions
  // Function for generating the key to use when storing elemen 'elem'.  (This is an 
  // implementation detail that should not worry users).

  static ElemKey generate_key(const double&, const double&);

  // ----------------------------------------------------
  // ---- CONSTRUCTORS, COPY, SWAP AND ASSIGNMENT -------
  // ----------------------------------------------------
  // Construct a LR-spline tensor-product surface with k-multiple 
  // knots at endpoints.
  // The coefficients are assumed to be stored sequentially, 
  // with the shortest stride in the u-direction.
  template<typename KnotIterator, typename CoefIterator>
  LRSplineSurface(int deg_u,
	   int deg_v,
	   int coefs_u,
	   int coefs_v,
	   int dimension,
	   KnotIterator knotvals_u_start,
	   KnotIterator knotvals_v_start,
	   CoefIterator coefs_start,
	   double knot_tol = 1.0e-8);

  // Construct a LR-spline tensor-product surface with k-multiple
  // knots at endpoints, and with all coefficients set to the origin
  // in the appropriate dimension.
  template<typename KnotIterator>
  LRSplineSurface(int deg_u,
	   int deg_v,
	   int coefs_u,
	   int coefs_v,
	   int dimension,
	   KnotIterator knotvals_u_start,
	   KnotIterator knotvals_v_start,
	   double knot_tol = 1.0e-8);

  // Construct a LRSplineSurface based on a spline surface
  LRSplineSurface(SplineSurface *surf, double knot_tol);

  // construct empty, invalid spline
  LRSplineSurface() 
    {
      curr_element_ = NULL;
    } 

  // Copy constructor
  LRSplineSurface(const LRSplineSurface& rhs);

    /// Virtual destructor, enables safe inheritance.
  virtual ~LRSplineSurface();

  // Assignment operator.
  const LRSplineSurface& operator= (const LRSplineSurface& other);


  // Constructor reading from an input stream
  LRSplineSurface(std::istream& is) { read(is);}

  // Swap operator (swap contents of two LRSplineSurfaces).
  void swap(LRSplineSurface& rhs);

  // ----------------------------------------------------
  // ------- READ AND WRITE FUNCTIONALITY ---------------
  // ----------------------------------------------------
  virtual void  read(std::istream& is);       
  virtual void write(std::ostream& os) const; 

  // ----------------------------------------------------
  // Inherited from GeomObject
  // ----------------------------------------------------
    virtual ClassType instanceType() const;

    static ClassType classType()
    { return Class_LRSplineSurface; }

    virtual LRSplineSurface* clone() const
    { return new LRSplineSurface(*this); }

    /// Return the spline surface represented by this surface, if any
    /// The user may get the spline surface lying in the (refined)
    /// regular grid by calling the function
    /// expandToFullTensorProduct().
    virtual SplineSurface* asSplineSurface(); 

  // inherited from GeomObject
  virtual BoundingBox boundingBox() const;

  // inherited from GeomObject
  // Dimension of function codomain (e.g. dimension of geometry space)
  virtual int dimension() const
  {
    // The weight is no longer part of the dimension.
    return bsplines_.size() > 0 ? 
      bsplines_.begin()->second->dimension() : 0;
      // bsplines_.begin()->second->dimension() - (int)rational_ : 0;
  }
    
   // ----------------------------------------------------
  // Inherited from ParamSurface
  // ----------------------------------------------------
    virtual const RectDomain& parameterDomain() const;

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

    // Fetch boundary curve
    // edge_num: 0 = umin, 1 = umax, 2 = vmin, 3 = vmax
    // The orientation of the curve is always from the smallest parameter
    // value to the largest
    SplineCurve* edgeCurve(int edge_num) const;

    // inherited from ParamSurface
    virtual void point(Point& pt, double upar, double vpar) const;

    void point(Point& pt, double upar, double vpar, Element2D* elem) const;

    // Output: Partial derivatives up to order derivs (pts[0]=S(u,v),
    // pts[1]=dS/du=S_u, pts[2]=S_v, pts[3]=S_uu, pts[4]=S_uv, pts[5]=S_vv, ...)
    // inherited from ParamSurface
    virtual void point(std::vector<Point>& pts, 
		       double upar, double vpar,
		       int derivs,
		       bool u_from_right = true,
		       bool v_from_right = true,
		       double resolution = 1.0e-12) const;

    void point(std::vector<Point>& pts, 
	       double upar, double vpar,
	       int derivs,
	       Element2D* elem,
	       bool u_from_right = true,
	       bool v_from_right = true,
	       double resolution = 1.0e-12) const;

    /// Closest point iteration taking benifit from information about
    /// an element in which to start searching
    void closestPoint(const Point& pt,
		      double& clo_u,
		      double& clo_v, 
		      Point& clo_pt,
		      double& clo_dist,
		      double epsilon,
		      int maxiter,
		      Element2D* elem = NULL,
		      const RectDomain* rd = NULL,
		      double *seed = NULL) const;

    /// Evaluate points in a grid
    /// The nodata value is applicable for bounded surfaces
    /// and will not be used in this context
    virtual void evalGrid(int num_u, int num_v, 
			  double umin, double umax, 
			  double vmin, double vmax,
			  std::vector<double>& points,
			  double nodata_val = -9999) const;

    /// Get the start value for the u-parameter
    /// \return the start value for the u-parameter
    virtual double startparam_u() const;

    /// Get the start value for the v-parameter
    /// \return the start value for the v-parameter
    virtual double startparam_v() const;

    /// Get the end value for the u-parameter
    /// \return the end value for the u-parameter
    virtual double endparam_u() const;

    /// Get the end value for the v-parameter
    /// \return the end value for the v-parameter
    virtual double endparam_v() const;

     // inherited from ParamSurface
    virtual void normal(Point& n, double upar, double vpar) const;

    void normal(Point& n, double upar, double vpar, Element2D* elem) const;

    /// Function that calls normalCone(NormalConeMethod) with method =
    /// SederbergMeyers. Needed because normalCone() is virtual! 
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

    // Fetch a part of a LRSplineSurface
    LRSplineSurface*
      subSurface(double from_upar, double from_vpar,
		 double to_upar, double to_vpar,
		 double fuzzy) const;

   // inherited from ParamSurface
    virtual std::vector<shared_ptr<ParamSurface> >
    subSurfaces(double from_upar, double from_vpar,
		double to_upar, double to_vpar,
		double fuzzy = DEFAULT_PARAMETER_EPSILON) const;

    /// Mirror a surface around a specified plane
    virtual LRSplineSurface* mirrorSurface(const Point& pos, const Point& norm) const;
    // inherited from ParamSurface
    virtual void closestBoundaryPoint(const Point& pt,
				      double&        clo_u,
				      double&        clo_v, 
				      Point&         clo_pt,
				      double&        clo_dist,
				      double         epsilon,
				      const RectDomain* rd = NULL,
				      double *seed = 0) const;
    // inherited from ParamSurface
    virtual void getBoundaryInfo(Point& pt1, Point& pt2, 
				 double epsilon, SplineCurve*& cv,
				 SplineCurve*& crosscv, double knot_tol = 1e-05) const;

    void
    getBoundaryInfo(double par1, double par2,
		    int bdindex, SplineCurve*& cv,
		    SplineCurve*& crosscv, double knot_tol) const;

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

    // inherited from ParamSurface
    virtual void turnOrientation();

    // inherited from ParamSurface
    virtual void swapParameterDirection();

    // inherited from ParamSurface
    virtual void reverseParameterDirection(bool direction_is_u);

    virtual void setParameterDomain(double u1, double u2, double v1, double v2);

    /// Compute the total area of this surface up to some tolerance
    /// \param tol the relative tolerance when approximating the area, i.e.
    ///            stop iteration when error becomes smaller than
    ///            tol/(surface area)
    /// \return the area calculated
    virtual double area(double tol) const;

    /// Return surface corners, geometric and parametric points
    /// in that sequence
    virtual void 
      getCornerPoints(std::vector<std::pair<Point,Point> >& corners) const;

    SplineCurve*
      constParamCurve(double parameter, bool pardir_is_u) const;


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

    // inherited from ParamSurface
    virtual std::vector< shared_ptr<ParamCurve> >
    constParamCurves(double parameter, bool pardir_is_u) const;

    bool isDegenerate(bool& b, bool& r,
		      bool& t, bool& l, double tolerance) const;

    /// Check for parallel and anti parallel partial derivatives in surface corners
    virtual void getDegenerateCorners(std::vector<Point>& deg_corners, double tol) const;

    // inherited from ParamSurface
    // This does not really have any meaning in s locally refined
    // surface unless we know the iso-parameter. The closest
    // possible approach is to use the global knot vector.
    virtual double nextSegmentVal(int dir, double par, bool forward, double tol) const;

    // ----------------------------------------------------
  // --------------- QUERY FUNCTIONS --------------------
  // ----------------------------------------------------
  // point evaluation
  Point operator()(double u, double v, int u_deriv = 0, int v_deriv = 0) const; // evaluation

  /* virtual void point(Point &pt, double u, double v, int iEl=-1) const; */
  /* virtual void point(Point &pt, double u, double v, int iEl, bool u_from_right, bool v_from_right) const; */
  /* virtual void point(std::vector<Point> &pts, double upar, double vpar,  */
  /* 		     int derivs, int iEl=-1) const; */
  Point operator()(double u, double v, int u_deriv, int v_deriv, 
		   Element2D* elem) const; // evaluation


  // Query parametric domain (along first (x) parameter: d = XFIXED; along second (y) parameter: YFIXED)
  double paramMin(Direction2D d) const;
  double paramMax(Direction2D d) const;

  // Query spline polynomial degree
  // @@sbr201302 Does the theory allow LRSplines with basis functions with different degrees.
  int degree(Direction2D d) const;

  // get a reference to the box partition (the underlying mesh)
  const Mesh2D& mesh() const {return mesh_;} 
  
  // Total number of separate basis functions defined over the box partition
  int numBasisFunctions() const {return (int)bsplines_.size();}

  int numElements() const {return (int)emap_.size();}

  // @@@ VSK. This functionality interface is fetched from the Trondheim code
  // We need a storage for last element evaluated. Index or reference?
  // Should the element be identified by index or reference?
  // How should the set of elements be traversed? Iterator?
  void computeBasis (double param_u, double param_v, BasisPtsSf     & result, int iEl=-1 ) const;
  void computeBasis (double param_u, double param_v, BasisDerivsSf  & result, int iEl=-1 ) const;
  void computeBasis (double param_u, double param_v, BasisDerivsSf2 & result, int iEl=-1 ) const;
  void computeBasis (double param_u,
		     double param_v,
		     std::vector<std::vector<double> >& result,
		     int derivs=0,
		     int iEl=-1 ) const;
  int getElementContaining(double u, double v) const;

  // Returns a pair.  
  // The first element of this pair is a pair of doubles representing the lower-left
  // corner of the element in which the point at (u, v) is located.  
  // The second element of this pair is a vector of pointers to the LRBSpline2Ds that cover
  // this element. (Ownership of the pointed-to LRBSpline2Ds is retained by the LRSplineSurface).
//  const ElementMap::value_type&
  Element2D*  coveringElement(double u, double v) const;

  // Construct a mesh of pointers to elements. The mesh has one entry for
  // each possible knot domain. If a knot has multiplicity zero in an area
  // several entries will point to the same element.
  // The construction can speed up evaluation in many points by making
  // it possible to avoid searching of the correct element
  void constructElementMesh(std::vector<Element2D*>& elements) const;
 
  // Returns pointers to all basis functions whose support covers the parametric point (u, v). 
  // (NB: ownership of the pointed-to LRBSpline2Ds is retained by the LRSplineSurface.)
  std::vector<LRBSpline2D*> basisFunctionsWithSupportAt(double u, double v) const;

  // Return an iterator to the beginning/end of the map storing the LRBSpline2Ds
  BSplineMap::const_iterator basisFunctionsBegin() const {return bsplines_.begin();}
  BSplineMap::const_iterator basisFunctionsEnd()   const {return bsplines_.end();}

#if 1
  BSplineMap::iterator basisFunctionsBeginNonconst() {return bsplines_.begin();}
  BSplineMap::iterator basisFunctionsEndNonconst() {return bsplines_.end();}
#endif

  BSplineMap::iterator bsplineFromDomain(double start_u, double start_v, double end_u,
					 double end_v, int startmult_u, int startmult_v,
					 int endmult_u, int endmult_v);

  std::vector<LRBSpline2D*> getBoundaryBsplines(Direction2D d, bool atstart);

  // The following function returns 'true' if the underlying mesh is a regular grid, i.e. 
  // the surface is a tensor product spline surface.
  bool isFullTensorProduct() const;

  /// Tolerance for equality of knots
  double getKnotTol()
  {
    return knot_tol_;
  }

  // ----------------------------------------------------
  // --------------- EDIT FUNCTIONS ---------------------
  // ----------------------------------------------------
  
  // Insert a single refinement in the mesh.  
  // If 'absolute' is set to 'false' (default), the refinement will _increment_ multiplicty
  // of the involved meshrectangles by 'mult'.  If set to 'true', the refinement 
  // will _set_ the multiplicity of the involved meshrectangles to 'mult' (however, if this results in
  // a _decrease_ of multiplicity for any involved meshrectangle, the method will throw an error instead).
  // The method will also throw an error if the resulting multiplicity for any meshrectangle would
  // end up being higher than degree+1.
  void refine(Direction2D d, double fixed_val, double start, double end, int mult = 1, bool absolute=false);

  // Same function as previous, but information about the refinement is passed along in a 'Refinement2D' structure
  // (defined above).  The 'absolute' argument works as in the previous refine() method.
  void refine(const Refinement2D& ref, bool absolute=false);

  // Insert a batch of refinements simultaneously.  The 'absolute' argument works as in the two 
  // preceding refine() methods.
  void refine(const std::vector<Refinement2D>& refs, bool absolute=false);

  // @@@ VSK. Index or iterator? Must define how the elements or bsplines 
  // are refined and call one of the other functions (refine one or refine
  // many). Is there a limit where one should be chosen before the other?
  // Refinement of one element really doesn't make sense, but it could b
  // a part of a larger strategy
  void refineBasisFunction(int index);
  void refineBasisFunction(const std::vector<int> &indices);
  void refineElement(int index);
  void refineElement(const std::vector<int> &indices);
  
  // Set the coefficient of the LRBSpline2D pointed to by 'target' to 'value'.
  // Conditions for calling this function are:
  // 1) 'target' should be a LRBSpline2D that belongs to the LRSplineSurface.
  // 2) The dimension of 'value' should be equal to the dimension of the LRSplineSurface image (e.g.
  //    the value returned by LRSplineSurface::dimension().
  void setCoef(const Point& value, const LRBSpline2D* target);

  void setCoefTimesGamma(const Point& value, const LRBSpline2D* target);

  // Set the coefficient of the LRBSpline2D with support as specified by the knots
  // with indices 'umin_ix', 'vmin_ix', 'umax_ix' and 'vmax_ix' in the mesh, and whose
  // knot multiplicities at the lower-left corner are indicated by u_mult and v_mult respectively.  
  // Throws if no such LRBSpline2D exists, or if the dimension of 'value' is not 
  // equal to LRSplineSurface::dimension().
  void setCoef(const Point& value, int umin_ix, int vmin_ix, int umax_ix, int vmax_ix, 
	       int u_mult = 1, int v_mult = 1);
  
  // Convert the LRSplineSurface to its full tensor product spline representation (NB: not reversible!)
  void expandToFullTensorProduct(); 

  // Add another LR B-spline surface to the current one.
  // NB! The surfaces must be defined on the same mesh and all scaling factors must
  // correspond. The function will throw if the requirements are not satisfied
  void addSurface(const LRSplineSurface& other_sf, double fac=1.0);

  // Convert a 1-D LRSplineSurface ("function") to a 3-D spline, by using the Greville points as x-
  // and y-coordinates, and the LRSplineSurface function value as z-coordinate.  
  // Requires that the LRSplineSurface is 1-D, and that the degree is > 0.
  void to3D();

  //Go::LineCloud getElementBds(int num_pts) const;

  bool rational() const;

  // Translate the surface along a given vector.
  void translate(const Point& vec);

  ElementMap::const_iterator elementsBegin() const { return emap_.begin();}
  ElementMap::const_iterator elementsEnd()   const { return emap_.end();}

  // Set element to be used as first try in point evaluation
  void setCurrentElement(Element2D* curr_el) 
  {
    curr_element_ = curr_el;
  }

  // ----------------------------------------------------
  // --------------- DEBUG FUNCTIONS --------------------
  // ----------------------------------------------------

  /* // Produce an ASCII image of the underlying mesh.  */
  /* // NB: requires a wide character stream (wostream). */
  /* void plotMesh(std::wostream& os) const; */

  /* // Produce a series of ASCII images showing the support of */
  /* // each basis function and its relation to the mesh. */
  /* // NB: requires a wide character stream (wostream). */
  /* void plotBasisFunctionSupports(std::wostream& os) const; */

  LineCloud getElementBds(int num_pts = 5) const;

 private:

  // ----------------------------------------------------
  // ----------------- PRIVATE DATA ---------------------
  // ----------------------------------------------------

  double knot_tol_;       // Tolerance for when to consider two knot values distinct rather than 
                          // a single one of higher multiplicity.

  bool rational_;

  Mesh2D mesh_;           // Represents mesh topology, multiplicites, as well as knot values.

  BSplineMap bsplines_;   // Map of individual b-spline basis functions.  

  ElementMap emap_;       // Map of individual elements

  // Generated data
  mutable RectDomain domain_;
  mutable Element2D* curr_element_;

   // Private constructor given mesh and LR B-splines
  LRSplineSurface(double knot_tol, bool rational,
		  Mesh2D& mesh, std::vector<std::unique_ptr<LRBSpline2D> >& b_splines);

#if 0
  // @@sbr Remove this when LRSplineEvalGrid does not need them any longer!
  ElementMap::iterator elementsBeginNonconst() { return emap_.begin();}
  ElementMap::iterator elementsEndNonconst() { return emap_.end();}
#endif

  // Locate all elements in a mesh
  static ElementMap construct_element_map_(const Mesh2D&, const BSplineMap&);

  // Collect all LR B-splines overlapping a specified area
//    std::vector<std::unique_ptr<LRBSpline2D> > 
    std::vector<LRBSpline2D*> 
    collect_basis(int from_u, int to_u, 
		  int from_v, int to_v) const;

    void 
      s1773(const double ppoint[],double aepsge, 
	    double estart[],double eend[],double enext[],
	    double gpos[], int maxiter,
	    Element2D* elem, int *jstat) const;

    //DEBUG
    void checkSupport(LRBSpline2D* basis) const;

}; 

// end class LRSplineSurface



}; // end namespace Go

#include "LRSplineSurface_impl.h"


#endif
