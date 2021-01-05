//===========================================================================
//                                                                           
// File: LRSplineVolume.h                                                    
//                                                                           
// Created: Mon Feb 25 11:08:53 2013                                         
//                                                                           
// Author: Sverre Briseid <Sverre.Briseid@sintef.no>
//                                                                           
// Revision: $Id$
//                                                                           
// Description:
//                                                                           
//===========================================================================

#ifndef _LRSPLINEVOLUME_H
#define _LRSPLINEVOLUME_H


#include "GoTools/trivariate/SplineVolume.h"
#include "GoTools/lrsplines2D/LRSplineSurface.h"
#include "GoTools/lrsplines3D/Direction3D.h"
#include "GoTools/lrsplines3D/LRBSpline3D.h"
#include "GoTools/lrsplines2D/BSplineUniLR.h"

#include <array>
#include <functional>
#include <set>
#include <map>
#include <vector>
#include <unordered_map>
#include <iostream> // @@ debug


namespace Go
{
// =============================================================================
class LRSplineVolume : public ParamVolume
// =============================================================================
{
   public:
  //
  // SEARCH STRUCTURES
  //
  // Structure representing a refinement to carry out.  This structure is used as an
  // argument to the LRSplineVolume::refine() methods below.  It is particularly useful
  // when passing along a whole batch of refinements. 
  struct Refinement3D
  {
    double kval;       // value of the fixed parameter of the meshrectangle to insert
    double start1;     // start value of the meshrectangle's non-fixed 1st parameter
    double end1;       // end value of the meshrectangle's non-fixed 1st parameter
    double start2;     // start value of the meshrectangle's non-fixed 2nd parameter
    double end2;       // end value of the meshrectangle's non-fixed 2nd parameter
    Direction3D d;     // direction of the meshrectangle (XDIR, YDIR or ZDIR)
    int multiplicity;  // multiplicity of the meshrectangle

    Refinement3D()
    {
      kval = 0.0;
      start1 = 0.0;
      end1 = 0.0;
      start2 = 0.0;
      end2 = 0.0;
      d = XDIR;
      multiplicity = -1;
    }
    
   Refinement3D(double val,
		double st1, double e1,
		double st2, double e2,
		Direction3D dir, int mult)
    {
      kval = val;
      start1 = st1;
      end1 = e1;
      start2 = st2;
      end2 = e2;
      d = dir;
      multiplicity = mult;
    }
    
    void setVal(double val,
		double st1, double e1,
		double st2, double e2,
		Direction3D dir, int mult)
    {
      kval = val;
      start1 = st1;
      end1 = e1;
      start2 = st2;
      end2 = e2;
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
      double u_min, v_min, w_min, u_max, v_max, w_max;
      int u_mult1, v_mult1, w_mult1, u_mult2, v_mult2, w_mult2;
    bool operator<(const BSKey& rhs) const; // needed for sorting when used in an STL map
  };

  // these maps could be redefined as hash tables later, as this is likely to improve
  // performance (at the expense of having to specify hash functions for these types of keys).
  typedef std::map<BSKey, std::unique_ptr<LRBSpline3D> > BSplineMap; // storage of basis functions


  struct double_pair_hash {
    size_t operator()(const std::pair<double, double>& dp) const {
      // Two hashes for independent variables, combined using XOR.  It is expected
      // that the resulting hash is probably as good as the input hashes.
      return std::hash<double>()(dp.first) ^ std::hash<double>()(dp.second);
    }
  };


  // Function for generating the key to use when storing B-spline function 'b'.  (This is an 
  // implementation detail that should not worry users).
  static BSKey generate_key(const LRBSpline3D& b, const Mesh3D& m);
  static BSKey generate_key(const LRBSpline3D& b);

  // The ElementMap will be used to keep track over which BasisFunctions are covering each
  // element.  An element is represented by its lower-left coordinates (doubles). 
  // Note: the same comment regarding the use of 'double' as a key in a map applies here 
  // (c.f. comment above on BSKey).
   struct ElemKey
  {
    double u_min, v_min, w_min;
    //inline 
    bool operator<(const ElemKey& rhs) const; // needed for sorting when used in an STL map
  };
    
  // these maps could be redefined as hash tables later, as this is likely to improve
  // performance (at the expense of having to specify hash functions for these types of keys).
  typedef std::map<ElemKey, std::unique_ptr<Element3D> > ElementMap; // storage of basis functions
  // Function for generating the key to use when storing elemen 'elem'.  (This is an 
  // implementation detail that should not worry users).

    static ElemKey generate_key(const double&, const double&, const double&);


  // ----------------------------------------------------
  // ---- CONSTRUCTORS, COPY, SWAP AND ASSIGNMENT -------
  // ----------------------------------------------------
  // Construct a LR-spline tensor-product surface with k-multiple 
  // knots at endpoints.
  // The coefficients are assumed to be stored sequentially, 
  // with the shortest stride in the u-direction.
  template<typename KnotIterator, typename CoefIterator>
  LRSplineVolume(int deg_u,
		 int deg_v,
		 int deg_w,
		 int coefs_u,
		 int coefs_v,
		 int coefs_w,
		 int dimension,
		 KnotIterator knotvals_u_start,
		 KnotIterator knotvals_v_start,
		 KnotIterator knotvals_w_start,
		 CoefIterator coefs_start,
		 double knot_tol = 1.0e-8);

  // Construct a LR-spline tensor-product surface with k-multiple
  // knots at endpoints, and with all coefficients set to the origin
  // in the appropriate dimension.
  template<typename KnotIterator>
  LRSplineVolume(int deg_u,
		 int deg_v,
		 int deg_w,
		 int coefs_u,
		 int coefs_v,
		 int coefs_w,
		 int dimension,
		 KnotIterator knotvals_u_start,
		 KnotIterator knotvals_v_start,
		 KnotIterator knotvals_w_start,
		 double knot_tol = 1.0e-8);

  // Construct a LRSplineVolume based on a spline surface
  LRSplineVolume(SplineVolume *surf, double knot_tol);

  // construct empty, invalid spline
  LRSplineVolume() {} 

  // Copy constructor
  LRSplineVolume(const LRSplineVolume& rhs);

  // Assignment operator.
  const LRSplineVolume& operator= (const LRSplineVolume& other);

  // Constructor reading from an input stream
  LRSplineVolume(std::istream& is) { read(is);}

  // ----------------------------------------------------
  // ------- READ AND WRITE FUNCTIONALITY ---------------
  // ----------------------------------------------------
  virtual void  read(std::istream& is);       
  virtual void write(std::ostream& os) const; 


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

  // Swap operator (swap contents of two LRSplineVolumes).
  void swap(LRSplineVolume& rhs);

  // ----------------------------------------------------
  // Inherited from GeomObject
  // ----------------------------------------------------
  virtual ClassType instanceType() const;

  static ClassType classType()
    { return Class_LRSplineVolume; }

  virtual LRSplineVolume* clone() const
    { return new LRSplineVolume(*this); }

  // ----------------------------------------------------
  // Inherited from ParamVolume
  // ----------------------------------------------------
  // umin, umax, vmin, vmax, wmin, wmax
  virtual const Array<double,6> parameterSpan() const;

  // inherited from ParamVolume
  virtual void point(Point& pt, double upar, double vpar, double wpar) const;

  void point(Point& pt, double upar, double vpar, double wpar, Element3D* elem) const;

  // Output: Partial derivatives up to order derivs (pts[0]=V(u,v,w),
  // pts[1]=dV/du=V_u, pts[2]=V_v, pts[3]=V_w, pts[4]=V_uu, pts[5]=V_uv,
  // pts[6]=V_uw, pts[7]=V_vv, pts[8]=V_vw, ...)
  // inherited from ParamVolume
  virtual void point(std::vector<Point>& pts, 
		     double upar, double vpar, double wpar,
		     int derivs,
		     bool u_from_right = true,
		     bool v_from_right = true,
		     bool w_from_right = true,
		     double resolution = 1.0e-12) const;

  void point(std::vector<Point>& pts, 
	     double upar, double vpar, double wpar,
	     int derivs,
	     Element3D* elem,
	     bool u_from_right = true,
	     bool v_from_right = true,
	     bool w_from_right = true,
	     double resolution = 1.0e-12) const;

  /// Grid evaluation in given element
  /// The sequence of points is (u1,v1,w1), (u2,v1,w1), ... (u1,v2,w1),
  /// (u2,v2,w1), ..., (u1,v1,w2), ...
  /// The function throws if any parameter is outside the element boundaries
  void elementGridEvaluate(Element3D *element,
			   std::vector<double>& upar,
			   std::vector<double>& vpar,
			   std::vector<double>& wpar,
			   std::vector<double>& points) const;

  /// Grid evaluation in given element. deriv number of derivatives are computed
  /// The sequence of points is (u1,v1,w1), (u2,v1,w1), ... (u1,v2,w1),
  /// (u2,v2,w1), ..., (u1,v1,w2), ...
  /// Derivatives are placed directly after the associated point. The sequence
  /// is F, Fu, Fv, Fw, Fuu, Fuv, Fuw, Fvv, Fvw, Fww, Fuuu, ...
  /// The function throws if any parameter is outside the element boundaries
  void elementGridEvaluate(Element3D *element,
			   double* upar, int usize,
			   double* vpar, int vsize,
			   double *wpar, int wsize,
			   int deriv,
			   std::vector<double>& points) const;

  /// Get the start value for the u-parameter
  /// \return the start value for the u-parameter
  virtual double startparam_u() const;

  /// Get the start value for the v-parameter
  /// \return the start value for the v-parameter
  virtual double startparam_v() const;

  /// Get the start value for the w-parameter
  /// \return the start value for the w-parameter
  virtual double startparam_w() const;

  /// Get the end value for the u-parameter
  /// \return the end value for the u-parameter
  virtual double endparam_u() const;

  /// Get the end value for the v-parameter
  /// \return the end value for the v-parameter
  virtual double endparam_v() const;

  /// Get the end value for the w-parameter
  /// \return the end value for the w-parameter
  virtual double endparam_w() const;

  // inherited from ParamVolume
  virtual DirectionCone tangentCone(int pardir) const;

  /// Get the derivative volume. Expresses the i,j,k-th derivative of
  /// a spline volume as a spline volume.
  /// \param ider1 what derivative to use along the u parameter
  /// \param ider2 what derivative to use along the v parameter
  /// \param ider3 what derivative to use along the w parameter
  /// \return pointer to a newly constructed SplineVolume which is 
  ///         the derivative volume.  User assumes ownership of this
  ///         object, and is responsible for destroying it.
  LRSplineVolume* derivVolume(int ider1, int ider2, int ider3) const;

  // Fetch a part of a LRSplineVolume
  LRSplineVolume*
  subVolume(double from_upar, double from_vpar, double from_wpar,
	    double to_upar, double to_vpar, double to_wpar,
	    double fuzzy) const;

    // Split volume in specified parameter values in a given direction
    std::vector<shared_ptr<SplineVolume> > 
      split(std::vector<double>& param,
	    int pardir,
	    double fuzzy = DEFAULT_PARAMETER_EPSILON) const; 

    // inherited from ParamVolume
    virtual void closestPoint(const Point& pt,
			      double&        clo_u,
			      double&        clo_v, 
			      double&        clo_w, 
			      Point&         clo_pt,
			      double&        clo_dist,
			      double         epsilon,
			      double   *seed = 0) const;

    /// Returns the corner closest to a given point together with
    /// the associated enumeration of the corner coefficient.
    /// In degenerate cases, the enumeration will reflect an arbitrary 
    /// corner
    int closestCorner(const Point& pt,
		      double& upar, double& vpar, double& wpar,
		      Point& corner, double& dist) const;

    /// Appends a volume to the end of 'this' volume along a specified parameter
    /// direction.
    /// \param vol the volume to append to 'this' volume
    /// \param join_dir the parameter direction along which to extend 'this' volume
    ///                 by appending the 'vol' volume. If it is equal to 0, we will
    ///                 extend the volume along the u-parameter; if it equals 1, we
    ///                 will extend the volume along the v-parameter; if it equals
    ///                 2, we will extend the volume along the w-parameter.
    /// \param cont continuity at jont, legal values are 0 and 1 (i,.e. C^0 or
    ///             C^1-continuity).
    /// \param repar The parametrization of the concerned parameter in the other 
    ///              volume will \em always be shifted so that it starts where the
    ///              parametrization of 'this' volume ends.  However, if 'repar' is 
    ///              set to 'true', it will also be \em scaled as a function of position
    ///              of control points close to the transition.
    void appendVolume(SplineVolume* vol, int join_dir,
		       int cont, bool repar=true);

    // inherited from ParamVolume
    virtual void reverseParameterDirection(int pardir);

    // inherited from ParamVolume
    virtual void swapParameterDirection(int pardir1, int pardir2);

    void setParameterDomain(double u1, double u2, double v1, double v2, double w1, double w2);

  /// Compute the total area of this surface up to some tolerance
  /// \param tol the relative tolerance when approximating the area, i.e.
  ///            stop iteration when error becomes smaller than
  ///            tol/(surface area)
  /// \return the area calculated
  virtual double volume(double tol) const;

  /// Return surface corners, geometric and parametric points
  /// in that sequence
  virtual void 
  getCornerPoints(std::vector<std::pair<Point,Point> >& corners) const;


  /// Generate and return a SplineSurface that represents a constant parameter 
  /// surface on the volume
  /// \param parameter value of the fixed parameter
  /// \param pardir 0 if the surface is constant in the u-parameter,
  ///               1 if the surface is constant in the v-parameter,
  ///               2 if the surface is constant in the w-parameter.
  /// \return pointer to a newly constructed SplineSurface representing the 
  ///         specified constant parameter surface.  It is the user's reponsibility
  ///         to delete it when it is no longer needed.
  virtual LRSplineSurface* constParamSurface(double parameter,
					     int pardir) const;

  void constParamSurfaces(std::vector<double>& param, int pardir,
			  std::vector<shared_ptr<LRSplineSurface> >& sfs) const;
  
  /// Fetch boundary curve
  /// bd_num: 0 = umin, 1 = umax, 2 = vmin, 3 = vmax, 4 = wmin, 5 = wmax
  LRSplineSurface* boundarySurface(int bd_num) const;

  /// Fetch all boundary surfaces corresponding to the volume.
  virtual std::vector<shared_ptr<ParamSurface> > 
  getAllBoundarySurfaces() const;

  /// Specify degeneracy
  /// \param which_sf 0=u_min, 1=u_max, 2=v_min, 3=v_max, 4=w_min, 5=w_max
  /// \param type 0 - Not degenerate, 1 - surface degenerate in specified boundary,
  ///             2 - surface degenerate to line, 3 - surface degenerate to point    
  /// \param b true if surface is degenerate along bottom curve (v=min for surface in parameterized 
  //                                                             u and v)
  /// \param r true if surface is degenerate along right curve
  /// \param t true if surface is degenerate along top curve
  /// \param l true if surface is degenerate along left curve
  /// \param tol Tolerance used in test
  bool isDegenerate(int which_sf, int& type, bool& b, bool& r,
		    bool& t, bool& l, double tol) const;

  /// Check for parallel and anti parallel partial derivatives in surface corners
  virtual void getDegenerateCorners(std::vector<Point>& deg_corners, double tol) const;

  // inherited from ParamVolume
  // This does not really have any meaning in s locally refined
  // surface unless we know the iso-parameter. The closest
  // possible approach is to use the global knot vector.
  virtual double nextSegmentVal(int dir, double par, bool forward, double tol) const;

  // ----------------------------------------------------
  // --------------- QUERY FUNCTIONS --------------------
  // ----------------------------------------------------
  // point evaluation
  Point operator()(double u, double v, double w, int u_deriv = 0, int v_deriv = 0, int w_deriv = 0) const; // evaluation

  /* virtual void point(Point &pt, double u, double v, int iEl=-1) const; */
  /* virtual void point(Point &pt, double u, double v, int iEl, bool u_from_right, bool v_from_right) const; */
  /* virtual void point(std::vector<Point> &pts, double upar, double vpar,  */
  /* 		     int derivs, int iEl=-1) const; */

  Point operator()(double u, double v, double w, int u_deriv, int v_deriv, int w_deriv, Element3D* elem) const;

  // Query parametric domain (along first (x) parameter: d = XDIR; along second (y) parameter: YDIR; third parameter ZDIR)
  double paramMin(Direction3D d) const;
  double paramMax(Direction3D d) const;

  // Query spline polynomial degree
  // @@sbr201302 Does the theory allow LRSplines with basis functions with different degrees.
  int degree(Direction3D d) const;

  // get a reference to the box partition (the underlying mesh)
  const Mesh3D& mesh() const {return mesh_;} 
  
  // Total number of separate basis functions defined over the box partition
  int numBasisFunctions() const {return (int)bsplines_.size();}

  int numElements() const {return (int)emap_.size();}

  // @@@ VSK. This functionality interface is fetched from the Trondheim code
  // We need a storage for last element evaluated. Index or reference?
  // Should the element be identified by index or reference?
  // How should the set of elements be traversed? Iterator?
  void computeBasis (double param_u, double param_v, double param_w,
		     BasisPtsSf     & result, int iEl=-1 ) const;
  void computeBasis (double param_u, double param_v, double param_w,
		     BasisDerivsSf  & result, int iEl=-1 ) const;
  void computeBasis (double param_u, double param_v, double param_w,
		     BasisDerivsSf2 & result, int iEl=-1 ) const;
  void computeBasis (double param_u,
		     double param_v,
		     double param_w,
		     std::vector<std::vector<double> >& result,
		     int derivs=0,
		     int iEl=-1 ) const;
  int getElementContaining(double u, double v, double w) const;

  // Returns a pair.  
  // The first element of this pair is a pair of doubles representing the lower-left
  // corner of the element in which the point at (u, v) is located.  
  // The second element of this pair is a vector of pointers to the LRBSpline3Ds that cover
  // this element. (Ownership of the pointed-to LRBSpline3Ds is retained by the LRSplineVolume).
  Element3D* coveringElement(double u, double v, double w) const;

  // Construct a mesh of pointers to elements. The mesh has one entry for
  // each possible knot domain. If a knot has multiplicity zero in an area
  // several entries will point to the same element.
  // The construction can speed up evaluation in many points by making
  // it possible to avoid searching of the correct element
  void constructElementMesh(std::vector<Element3D*>& elements) const;
 
  // Returns pointers to all basis functions whose support covers the parametric point (u, v). 
  // (NB: ownership of the pointed-to LRBSpline3Ds is retained by the LRSplineVolume.)
    std::vector<LRBSpline3D*> basisFunctionsWithSupportAt(double u, double v, double w) const;

  // Return an iterator to the beginning/end of the map storing the LRBSpline3Ds
  BSplineMap::const_iterator basisFunctionsBegin() const {return bsplines_.begin();}
  BSplineMap::const_iterator basisFunctionsEnd()   const {return bsplines_.end();}

#if 1
  BSplineMap::iterator basisFunctionsBeginNonconst() {return bsplines_.begin();}
  BSplineMap::iterator basisFunctionsEndNonconst() {return bsplines_.end();}
#endif

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
  void refine(Direction3D d, double fixed_val,
	      double start1, double end1,
	      double start2, double end2,
	      int mult = 1, bool absolute=false);

  // Same function as previous, but information about the refinement is passed along in a 'Refinement3D' structure
  // (defined above).  The 'absolute' argument works as in the previous refine() method.
  void refine(const Refinement3D& ref, bool absolute=false);

  // Insert a batch of refinements simultaneously.  The 'absolute' argument works as in the two 
  // preceding refine() methods.
  void refine(const std::vector<Refinement3D>& refs, bool absolute=false);
#if 0
  // Insert a zero multiplicity knot 
  void zero_knot(Direction3D d, double knotval);
#endif  
  // @@@ VSK. Index or iterator? Must define how the elements or bsplines 
  // are refined and call one of the other functions (refine one or refine
  // many). Is there a limit where one should be chosen before the other?
  // Refinement of one element really doesn't make sense, but it could b
  // a part of a larger strategy
  void refineBasisFunction(int index);
  void refineBasisFunction(const std::vector<int> &indices);
  void refineElement(int index);
  void refineElement(const std::vector<int> &indices);
  
  // Set the coefficient of the LRBSpline3D pointed to by 'target' to 'value'.
  // Conditions for calling this function are:
  // 1) 'target' should be a LRBSpline3D that belongs to the LRSplineVolume.
  // 2) The dimension of 'value' should be equal to the dimension of the LRSplineVolume image (e.g.
  //    the value returned by LRSplineVolume::dimension().
  void setCoef(const Point& value, const LRBSpline3D* target);

  // Set the coefficient of the LRBSpline3D with support as specified by the knots
  // with indices 'umin_ix', 'vmin_ix', 'umax_ix' and 'vmax_ix' in the mesh, and whose
  // knot multiplicities at the lower-left corner are indicated by u_mult and v_mult respectively.  
  // Throws if no such LRBSpline3D exists, or if the dimension of 'value' is not 
  // equal to LRSplineVolume::dimension().
  void setCoef(const Point& value,
	       int umin_ix, int vmin_ix, int wmin_ix, int umax_ix, int vmax_ix, int wmax_ix,
	       int u_mult = 1, int v_mult = 1, int w_mult = 1);
  
  // Convert the LRSplineVolume to its full tensor product spline representation (NB: not reversible!)
  void expandToFullTensorProduct(); 

  // Add another LR B-spline volume to the current one.
  // NB! The volumes must be defined on the same mesh and all scaling factors must
  // correspond. The function will throw if the requirements are not satisfied
  void addVolume(const LRSplineVolume& other_vol, double fac=1.0);

  // Convert a 1-D LRSplineVolume ("function") to a 3-D spline, by using the Greville points as x-
  // and y-coordinates, and the LRSplineVolume function value as z-coordinate. This is incorrect,
  // except for locally linearly independent LR-splines, but generally gives a good approximation 
  // if the gamma values do not vary too much.
  // Requires that the LRSplineVolume is 1-D, and that the degree is > 0.
  void to3D();

  bool rational() const;

  /// Return the spline surface represented by this surface, if any
  /// The user may get the spline surface lying in the (refined)
  /// regular grid by calling the function
  /// expandToFullTensorProduct().
  virtual SplineVolume* asSplineVolume(); 

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

  // Set element to be used as first try in point evaluation
  void setCurrentElement(Element3D* curr_el)
  {
    curr_element_ = curr_el;
  }

  ElementMap::const_iterator elementsBegin() const { return emap_.begin();}
  ElementMap::const_iterator elementsEnd()   const { return emap_.end();}

  /// Check if the volume is of type spline
  virtual bool isSpline() const
    {
      return false;  // This is a spline
    }

  // inherited from ParamVolume
  virtual void translate(const Point& vec);

  /// Test if the volume is lefthanded. Returns false if volume does
  /// not lie in 3-dimensional space. The test is done by taking the
  /// XXX of the derivatives at the centre point. Assumes volume is
  /// not self intersecting, and has linearly independent derivatives
  /// at the centre.
  bool isLeftHanded();


private:

  // ----------------------------------------------------
  // ----------------- PRIVATE DATA ---------------------
  // ----------------------------------------------------

  double knot_tol_;       // Tolerance for when to consider two knot values distinct rather than 
                          // a single one of higher multiplicity.

  bool rational_;

  Mesh3D mesh_;           // Represents mesh topology, multiplicites, as well as knot values.

  // Map of individual univariate b-spline basis functions, 1. par. dir.  
  std::vector<std::unique_ptr<BSplineUniLR> > bsplinesuni1_;  // To be kept sorted   

  // Map of individual univariate b-spline basis functions, 2. par. dir.  
  std::vector<std::unique_ptr<BSplineUniLR> > bsplinesuni2_;   

  // Map of individual univariate b-spline basis functions, 3. par. dir.  
  std::vector<std::unique_ptr<BSplineUniLR> > bsplinesuni3_;   

  BSplineMap bsplines_;   // Map of individual b-spline basis functions.  

  ElementMap emap_;       // Map of individual elements

  // Generated data
  mutable Array<double,6> domain_;
  mutable Element3D* curr_element_;

  // Private constructor given mesh and LR B-splines
  LRSplineVolume(double knot_tol, bool rational,
                 Mesh3D& mesh, std::vector<LRBSpline3D*> b_splines,
		 int first_ixu, int first_ixv, int first_ixw);

#if 0
  // @@sbr Remove this when LRSplineEvalGrid does not need them any longer!
  ElementMap::iterator elementsBeginNonconst() { return emap_.begin();}
  ElementMap::iterator elementsEndNonconst() { return emap_.end();}
#endif

  // Refine volume without reconstructing the element map
  void doRefine(const std::vector<Refinement3D>& refs, bool absolute);

  // Define constant parameter surface given a degree+1 multiple knot at parameter
  LRSplineSurface* defineBivariate(int pardir, double parval, bool at_start) const;

  // Locate all elements in a mesh
  static ElementMap construct_element_map_(const Mesh3D&, const BSplineMap&);



  // Collect all LR B-splines overlapping a specified area
  std::vector<LRBSpline3D* >
  collect_basis(int from_u, int to_u, 
		int from_v, int to_v,
		int from_w, int to_w) const;

}; 

// end class LRSplineVolume






}; // end namespace Go

#include "LRSplineVolume_impl.h"

#endif // _LRSPLINEVOLUME_H

