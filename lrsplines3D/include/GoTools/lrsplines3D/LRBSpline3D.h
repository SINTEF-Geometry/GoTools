//===========================================================================
//                                                                           
// File: LRBSpline3D.h                                                       
//                                                                           
// Created: Mon Feb 25 11:07:24 2013                                         
//                                                                           
// Author: Sverre Briseid <Sverre.Briseid@sintef.no>
//                                                                           
// Revision: $Id$
//                                                                           
// Description:
//                                                                           
//===========================================================================

#ifndef _LRBSPLINE3D_H
#define _LRBSPLINE3D_H

#include <vector>
#include <assert.h>

#include "GoTools/utils/Point.h"
#include "GoTools/lrsplines3D/Direction3D.h"
#include "GoTools/lrsplines3D/Element3D.h"
#include "GoTools/lrsplines3D/Mesh3D.h"
#include "GoTools/lrsplines2D/BSplineUniLR.h"
#include "GoTools/utils/StreamUtils.h"
#include "GoTools/geometry/Streamable.h"

namespace Go
{



//==============================================================================
/// This class represents a single LR B-spline basis function, intended for 
/// use in a LRSplineVolume.  It contains the following information
/// Its coefficient
/// Its gmama multplier (scaling factor to ensure partition of unity (c.f. LR-spline paper)
/// Its three knot vectors represented as instances of  BSplineUniLR

class LRBSpline3D : public Streamable
//==============================================================================
{
 public:

  // ---------------------------------------------------------
  // --- CONSTRUCTORS, READING, WRITING AND SWAP FUNCTIONS --- 
  // ---------------------------------------------------------

  /// Constructor to create an empty (invalid) LRBSpline3D
  LRBSpline3D() 
    {}

  //template<typename Iterator>
  /// Constructor
  /// \param c_g coefficient
  /// \param weight rational weight (1 if non-rational volume)
  /// \param bspline_u univariate B-spline in first parameter direction
  /// \param bspline_v univariate B-spline in second second direction
  /// \param bspline_v univariate B-spline in third second direction
  /// \param gamma scaling factor
  /// \param rational indicates rational surface
  LRBSpline3D(const Point& c_g, double weight,
	      BSplineUniLR *bspline_u,
	      BSplineUniLR *bspline_v,
	      BSplineUniLR *bspline_w,
	      double gamma, 
	      bool rational = false)
    : coef_times_gamma_(c_g),
    weight_(weight),
    rational_(rational),
    gamma_(gamma),
    bspline_u_(bspline_u),
    bspline_v_(bspline_v),
    bspline_w_(bspline_w),
    coef_fixed_(0)
    {
      bspline_u_->incrCount();
      bspline_v_->incrCount();
      bspline_w_->incrCount();
    }

  /// Copy constructor
  LRBSpline3D(const LRBSpline3D& rhs);

  /// Swap the contents of two LRBSpline3Ds
  /// @@sbr201303 What about the degrees? Mesh? The rest?
  void swap(LRBSpline3D& rhs) 
  {
    coef_times_gamma_.swap(rhs.coef_times_gamma_);
    std::swap(weight_, rhs.weight_);
    std::swap(rational_, rhs.rational_);
    std::swap(gamma_, rhs.gamma_);
    std::swap(bspline_u_, rhs.bspline_u_);
    std::swap(bspline_v_, rhs.bspline_v_);
    std::swap(bspline_w_, rhs.bspline_w_);
    std::swap(coef_fixed_,rhs.coef_fixed_);
 }

  /// Destructor
  ~LRBSpline3D() 
    {  
      //std::cout << "Delete LRBSpline " << this << std::endl;
      bspline_u_->decrCount();
      bspline_v_->decrCount();
      bspline_w_->decrCount();
   }

  /// Write the LRBSpline3D to a stream
  virtual void write(std::ostream& os) const;
  
  /// Read the LRBSpline3D from a stream
  /// Do not use this function. Will create memory loss
  virtual void read(std::istream& is);

  /// Read the LRBSpline3D from a stream, and collect univariate B-splines
  void read(std::istream& is, 
	    std::vector<std::unique_ptr<BSplineUniLR> >& bsplineuni_u,
	    int& left1, 
	    std::vector<std::unique_ptr<BSplineUniLR> >& bsplineuni_v,
	    int& left2,
	    std::vector<std::unique_ptr<BSplineUniLR> >& bsplineuni_w,
	    int& left3);

  // ---------------------------
  // --- EVALUATION FUNCTION ---
  // ---------------------------

  /// Read the LRBSpline2D from a stream
  /// Do not use this function. Will create memory loss
  double evalBasisFunc(double u, double v, double w) const;

  /// Similar to 'eval' below, but returns the value of the
  /// LRBSpline3D's underlying _basis_ function, rather than the
  /// function value itself. (In other words, the basis function's
  /// value is not multiplied by the spline coefficient, nor gamma,
  /// before it is returned.) If u_deriv = n1, v_deriv = n2 and
  /// w_deriv = n3, the function returns the partial derivative
  /// d^(n1+n2+n3) B / du^n1 dv^n2 dw^n3
  // @@@ VSK. Appropriate for rationals?
  /* double evalBasisFunction(double u, double v, double w, */
  /*                          const double* const kvals_u, */
  /*                          const double* const kvals_v, */
  /*                          const double* const kvals_w, */
  /* 			   int u_deriv = 0, int v_deriv = 0, int w_deriv = 0, */
  /* 			   bool u_at_end = false, bool v_at_end = false,  */
  /* 			   bool w_at_end = false) const; */

  double evalBasisFunction(double u, double v, double w,
                           int u_deriv = 0, int v_deriv = 0, int w_deriv = 0,
                           bool u_at_end = false, bool v_at_end = false, 
			   bool w_at_end = false) const;

  /// Evaluate the LRBSpline3D or its derivative in (u, v, w), looking
  /// up the knot values from the arrays pointed to by 'kvals_u', 'kvals_v' and
  /// 'kvals_w' (the actual indices to the relevant knots are already
  /// stored in the BSplineUniLR). It returns the partial derivative
  /// computed by evalBasisFuntion times the coefficient times the scaling factor.
  /// If the basis function is positioned at the upper boundary of the
  /// domain in either the u, v or w direction, the consideration of
  /// half-open intervals is reversed in order to cover the closure of
  /// the domain. The basis function itself has no knowledge of the
  /// underlying mesh, so this information has to be given explicitly
  /// using the 'u_at_end', 'v_at_end' and 'w_at_end' parameters.
  // @@@ VSK. I am not sure if I like the way of distributing position and
  // derivatives, but it is not top level. What about mixed derivatives?
  // What about rationals? Should maybe look at SplineSurface for interface.


  Point eval(double u, double v, double w,
             int u_deriv = 0, int v_deriv = 0, int w_deriv = 0,
             bool u_at_end = false, bool v_at_end = false, bool w_at_end = false) const
  {
    if (rational_)
      {
        //MESSAGE("Rational case under construction!");
        if (u_deriv > 0 || v_deriv > 0 || w_deriv > 0)
          MESSAGE("Rational case with derivs not yet supported!");
      }

    Point geom_pos =
      evalBasisFunction(u, v, w, u_deriv, v_deriv, w_deriv, u_at_end, v_at_end, w_at_end) *
      coefTimesGamma();

    return geom_pos;

  }

  /* Point eval(double u, double v, double w, */
  /* 	     const double* const kvals_u, const double* const kvals_v, const double* const kvals_w, */
  /* 	     int u_deriv = 0, int v_deriv = 0, int w_deriv = 0, */
  /* 	     bool u_at_end = false, bool v_at_end = false, bool w_at_end = false) const */
  /* {  */
  /*   if (rational_) */
  /*     { */
  /* 	MESSAGE("Rational case under construction!"); */
  /* 	if (u_deriv > 0 || v_deriv > 0 || w_deriv > 0) */
  /* 	  MESSAGE("Rational case with derivs not yet supported!"); */
  /*     } */


  /*   Point geom_pos = */
  /*     evalBasisFunction(u, v, w, kvals_u, kvals_v, kvals_w, u_deriv, v_deriv, w_deriv, u_at_end, v_at_end, w_at_end) *  */
  /*     coefTimesGamma(); */

  /*   return geom_pos; */
      
  /* } */

  /// Evaluate position only, array of correct size is expected as
  /// input. No size checking
  void evalpos(double u, double v, double w, double pos[])
  {
    double bb = evalBasisFunc(u, v, w);
    int dim = coef_times_gamma_.dimension();
    for (int ki=0; ki<dim; ++ki)
      pos[ki] = bb*coef_times_gamma_[ki];
  }

  /// Evaluate positon and requested derivatives and add the result to the
  /// existing content in the output array
  /// \param u parameter value in the first direction
  /// \param v parameter value in the second direction
  /// \param w parameter value in the third direction
  /// \param deriv number of derivative to evaluate (0 = only position)
  /// \param der array of length (deriv+1)*(deriv+2)*(deriv+3)/6 for accumulated 
  /// storage of derivatives. Sequence: p, du, dv, dw, duu, duv, duw, dvv, dvw, dww, ...
  void evalder_add(double u, double v, double w,
		   int deriv,
		   Point der[],
		   bool u_at_end, bool v_at_end, bool w_at_end) const;
 
  // -----------------------
  // --- QUERY FUNCTIONS ---
  // -----------------------

  /// Access the coefficient multiplied by the gamma factor
  // (to get the pure coefficient,
  // divide by the gamma factor, which can be obtained by the gamma() member function below).
        Point& coefTimesGamma()       { return coef_times_gamma_;}
  /// Access the coefficient multiplied by the gamma factor
  const Point& coefTimesGamma() const { return coef_times_gamma_;}

 /// Access the coefficient 
        Point Coef()       { return coef_times_gamma_/gamma_;}
 /// Access the coefficient 
  const Point Coef() const { return coef_times_gamma_/gamma_;}
 
  /// Access the gamma multiplier of this LRBSpline3D (should be set to ensure partition 
  /// of one, c.f. LR-spline paper).
        double& gamma()       {return gamma_;}
  const double& gamma() const {return gamma_;}

  /// Access the rational weight of this LRBSpline3D.
        double& weight()       {return weight_;}
  const double& weight() const {return weight_;}

  /// Check if this B-spline is rational
  const bool rational() const {return rational_;}

  /// Get the dimension of the LRBSpline3Ds codomain.
  /// For rational cases the dimension is the same, i.e. interpreted as geometric dimension.
  const int dimension() const {return coef_times_gamma_.dimension();}

  /// Access the LRBSpline3D's knot vector in the given direction.  (The knot vectors
  /// only contain incices to an external, shared vector of knot values).
  /// To get the knotvector in the first direction (x-direction), 'd' should be XDIR.
  /// To get the knotvector in the second direction (y-direction), 'd' should be YDIR.
  /// To get the knotvector in the third direction (z-direction), 'd' should be ZDIR.
  const std::vector<int>& kvec(Direction3D d) const
  {return (d==XDIR) ? bspline_u_->kvec() : ((d==YDIR) ? bspline_v_->kvec() : 
					    bspline_w_->kvec());}
  std::vector<int>& kvec(Direction3D d)
    {return (d==XDIR) ? bspline_u_->kvec() : ((d==YDIR) ? bspline_v_->kvec() : 
					      bspline_w_->kvec());}

  /// Get the polynomial degree of the spline.
  const int degree(Direction3D d) const {return (int)kvec(d).size() - 2;}  

  /// Get the index to the knot that defines the start of the LRBSpline3D's support.
  // (The vector of the actual knot values is stored outside of the LRBSpline3D, as it 
  // is shared among many LRBSpline3Ds).
  const int suppMin(Direction3D d) const {return kvec(d).front();}
  /// Get the index to the knot that defines the end of the LRBSpline2D's support.
  const int suppMax(Direction3D d) const {return kvec(d).back();}

  /// Information about the domain covered by this B-spline
  double umin() const 
  { 
    return bspline_u_->min();
  }
  double umax() const 
  { 
    return bspline_u_->max();
  }
  double vmin() const 
  { 
    return bspline_v_->min();
  }
  double vmax() const 
  {
    return bspline_v_->max();
  }
  double wmin() const 
  { 
    return bspline_w_->min();
  }
  double wmax() const 
  {
    /* std::vector<int> kv = bspline_w_->kvec(); */
    /* if (kv[kv.size()-1] > 18) */
    /*   std::cout << "wmax: " << kv[kv.size()-1] << std::endl; */
    return bspline_w_->max();
  }

  /// Query if the coefficient is fixed in volume approximation
  int coefFixed() const
  {
    return coef_fixed_;
  }

  /// Specify if the coefficient should be fixed in volume approximation
  void setFixCoef(int coef_fixed)
  {
    coef_fixed_ = coef_fixed;
  }

  /// Count multiplicity in the ends of the B-spline, first parameter direction
  int endmult_u(bool atstart) const;
  /// Count multiplicity in the ends of the B-spline, second parameter direction
  int endmult_v(bool atstart) const;
  /// Count multiplicity in the ends of the B-spline, third parameter direction
  int endmult_w(bool atstart) const;
  /// Count multiplicity in the start/end of the B-spline, specified parameter direction
  int endmult(Direction3D dir, bool atstart) const;

  /// Query whether the parameter point speficied by the knots indexed by 'u_ix', 'v_ix' 
  /// and 'w_ix' is covered by the support of this LRBSpline3D.
  // (NB: The vector of the actual knot 
  // values is stored outsde of the LRBSpline3D, since it is shared among many 
  // LRBSpline3Ds.
  bool coversCorner(int u_ix, int v_ix, int w_ix) const { 
    return 
      bspline_u_->coversPar(u_ix) && bspline_v_->coversPar(v_ix) &&
      bspline_w_->coversPar(w_ix);
  }

  /// Return greville parameter in the three parameter directions
  Point getGrevilleParameter() const;
  /// Return greville parameter in the specified parameter direction
  double getGrevilleParameter(Direction3D d) const
  {
    return (d == XDIR) ? bspline_u_->getGrevilleParameter() :
      ((d == YDIR) ? bspline_v_->getGrevilleParameter() :
       bspline_w_->getGrevilleParameter());
  }



  /// Fetch univariate B-spline
  BSplineUniLR* getUnivariate(Direction3D d) const
  {
    return (d == XDIR) ? bspline_u_ : (d == YDIR) ? bspline_v_ : bspline_w_;
  }

  /// Update univariate B-spline pointer
  void setUnivariate(Direction3D d, BSplineUniLR *uni)
  {
    if (d == XDIR)
      {
	if (bspline_u_ != NULL)
	  bspline_u_->decrCount();
	bspline_u_ = uni;
	bspline_u_->incrCount();
      }
    else if (d == YDIR)
      {
	if (bspline_v_ != NULL)
	  bspline_v_->decrCount();
	bspline_v_ = uni;
	bspline_v_->incrCount();
      }
    else
      {
	if (bspline_w_ != NULL)
	  bspline_w_->decrCount();
	bspline_w_ = uni;
	bspline_w_->incrCount();
      }
  }

  // Operations related to the support of this B-spline
  /// Check if the support of this B-spline overlaps the given element
  bool overlaps(Element3D *el) const;
  /// Check if the support of this B-spline overlaps the given domain: umin, umax, vmin, vmax, wmin, wmax.
  bool overlaps(double domain[]) const; // domain: umin, umax, vmin, vmax, wmin, wmax.

  /// Add element to vector of elements in the support
  bool addSupport(Element3D *el) ;
  /// Remove element from vector of elements in the support
  void removeSupport(Element3D *el) ;
  /// Remove all elements from vector of elements in the support
  void removeSupportedElements()
  {
    support_.clear();
  }
  
  /// Iterator to start of elements in the support
  std::vector<Element3D*>::iterator supportedElementBegin() ;
  /// Iterator to end of elements in the support
  std::vector<Element3D*>::iterator supportedElementEnd() ;
#if 0
  std::vector<Element3D*> getExtendedSupport() ;
  std::vector<Element3D*> getMinimalExtendedSupport();
#endif
  /// All elements in the support
  std::vector<Element3D*> supportedElements()
    {
      return support_;
    }
  /// Set all elements in the support
  void setSupport(std::vector<Element3D*> elements)
  {
    support_ = elements;
  }

  /// Number of elements in the support
  int nmbSupportedElements() { return (int)support_.size(); }

  /// Associate LR mesh to this B-spline
  void setMesh(const Mesh3D* mesh)
  {
    bspline_u_->setMesh(mesh);
    bspline_v_->setMesh(mesh);
    bspline_w_->setMesh(mesh);
  }

  /// Fetch associated LR mesh
  const Mesh3D* getMesh()
  {
    return dynamic_cast<const Mesh3D*>(bspline_u_->getMesh());
  }

  // -----------------
  // --- OPERATORS ---
  // -----------------

  /// Operator defining a partial ordering of LRBSpline3Ds.
  bool operator<(const LRBSpline3D& rhs) const;

  /// Equality operator
  bool operator==(const LRBSpline3D &rhs) const;

 /// Set coefficient and scaling factor (gamma)
  void setCoefAndGamma(Point& coef, double gamma)
    {
      gamma_ = gamma;
      coef_times_gamma_ = coef*gamma;;
    }

  void reverseParameterDirection(int pardir);

  void swapParameterDirection(int pardir1, int pardir2);

 private:

  Point coef_times_gamma_;
  double weight_; // For rational case.
  bool rational_; // @@sbr201301 Should this also be part of the LRSplineVolume? It seems
                  // best suited for this class since it is here we use the rational part.
  double gamma_; // normalizing weight to ensure partition of unity, c.f. Section 7 of paper
  BSplineUniLR *bspline_u_;
  BSplineUniLR *bspline_v_;
  BSplineUniLR *bspline_w_;
  std::vector<Element3D*> support_;  // Elements lying in the support of this LRB-spline
  
  // Used in approximation algorithms
  int coef_fixed_;  // 0=free coefficients, 1=fixed, 2=not affected

}; // end class LRBSpline3D

 inline std::ostream& operator<<(std::ostream& os, const LRBSpline3D& b) {b.write(os); return os;}
 inline std::istream& operator>>(std::istream& is,       LRBSpline3D& b) {b.read (is); return is;}

}; // end namespace Go


#endif // _LRBSPLINE3D_H

