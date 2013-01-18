#ifndef LRBSPLINE2D_H
#define LRBSPLINE2D_H

#include <vector>
#include <assert.h>

#include "GoTools/utils/Point.h"
#include "GoTools/lrsplines2D/Direction2D.h"
#include "GoTools/lrsplines2D/Element2D.h"
#include "GoTools/lrsplines2D/Mesh2D.h"
#include "GoTools/utils/StreamUtils.h"
#include "GoTools/geometry/Streamable.h"

namespace Go
{


//==============================================================================
// This class represents a single LR B-spline basis function, intended for 
// use in a LRSplineSurface.  It contains the following information
// Its coefficient
// Its gmama multplier (scaling factor to ensure partition of unity (c.f. LR-spline paper)
// Its two knot vectors (only the indices to the knots, not the values themselves, as these
// are usually shared among many LRBSpline2Ds.  The knotvalues are therefore 
// stored separately in the LRSpline, for collective lookup.)

class LRBSpline2D : public Streamable
//==============================================================================
{
 public:

  // ---------------------------------------------------------
  // --- CONSTRUCTORS, READING, WRITING AND SWAP FUNCTIONS --- 
  // ---------------------------------------------------------

  /// Constructor to create an empty (invalid) LRBSpline2D
  LRBSpline2D() 
    { }; 

  template<typename Iterator>
  LRBSpline2D(const Point& c_g, double weight,
	      int deg_u, int deg_v, 
	      Iterator kvec_u_start, Iterator kvec_v_start,
	      double gamma, const Mesh2D* mesh,
	      bool rational = false) 
    : coef_times_gamma_(c_g),
      weight_(weight),
      rational_(rational),
      gamma_(gamma),
      kvec_u_(kvec_u_start, kvec_u_start + deg_u + 2),
      kvec_v_(kvec_v_start, kvec_v_start + deg_v + 2),
      mesh_(mesh)
    {}

  /// Swap the contents of two LRBSpline2Ds
  void swap(LRBSpline2D& rhs) 
  {
    coef_times_gamma_.swap(rhs.coef_times_gamma_);
    std::swap(weight_, rhs.weight_);
    std::swap(rational_, rhs.rational_);
    std::swap(gamma_, rhs.gamma_);
    kvec_u_.swap(rhs.kvec_u_);
    kvec_v_.swap(rhs.kvec_v_);
    //    mesh_.swap(rhs.mesh_);
  }

  ~LRBSpline2D() 
    {  
      //std::cout << "Delete LRBSpline " << this << std::endl;
    }; 

  /// Write the LRBSpline2D to a stream
  virtual void write(std::ostream& os) const;
  
  /// Read the LRBSpline2D from a stream
  virtual void read(std::istream& is);

  // ---------------------------
  // --- EVALUATION FUNCTION ---
  // ---------------------------

  /// Similar to 'eval' below, but returns the value of the
  /// LRBSpline2D's underlying _basis_ function, rather than the
  /// function value itself. (In other words, the basis function's
  /// value is not multiplied by the spline coefficient, nor gamma,
  /// before it is returned.)
  // @@@ VSK. Appropriate for rationals?
  double evalBasisFunction(double u, double v, 
			   const double* const kvals_u, 
			   const double* const kvals_v,
			   int u_deriv = 0, int v_deriv = 0,
			   bool u_at_end = false, bool v_at_end = false) const;

  // Evaluate the LRBSpline2D or its derivative in (u, v), looking
  // up the knot values from the arrays pointed to by 'kvals_u' and
  // 'kvals_v' (the actual indices to the relevant knots are already
  // stored in the LRBSpline2D). If u_deriv = n and v_deriv = m,
  // the partial derivative d^(n+m) B / du^n dv^m will be computed. If
  // the basis function is positioned at the upper boundary of the
  // domain in either the u or v direction, the consideration of
  // half-open intervals is reversed in order to cover the closure of
  // the domain. The basis function itself has no knowledge of the
  // underlying mesh, so this information has to be given explicitly
  // using the 'u_at_end' and 'v_at_end' parameters.
  // @@@ VSK. I am not sure if I like the way of distributing position and
  // derivatives, but it is not top level. What about mixed derivatives?
  // What about rationals? Should maybe look at SplineSurface for interface.
  Point eval(double u, double v, 
	     const double* const kvals_u, const double* const kvals_v,
	     int u_deriv = 0, int v_deriv = 0,
	     bool u_at_end = false, bool v_at_end = false) const
  { 
    if (rational_)
      {
	MESSAGE("Rational case under construction!");
	if (u_deriv > 0 || v_deriv > 0)
	  MESSAGE("Rational case with derivs not yet supported!");
      }


    Point geom_pos =
      evalBasisFunction(u, v, kvals_u, kvals_v, u_deriv, v_deriv, u_at_end, v_at_end) * 
      coefTimesGamma();

    return geom_pos;
      
  }

  // -----------------------
  // --- QUERY FUNCTIONS ---
  // -----------------------

  // Access the coefficient multiplied by the gamma factor (to get the pure coefficient,
  // divide by the gamma factor, which can be obtained by the gamma() member function below).
        Point& coefTimesGamma()       { return coef_times_gamma_;}
  const Point& coefTimesGamma() const { return coef_times_gamma_;}

        Point Coef()       { return coef_times_gamma_/gamma_;}
  const Point Coef() const { return coef_times_gamma_/gamma_;}
 
  // Access the gamma multiplier of this LRBSpline2D (should be set to ensure partition 
  // of one, c.f. LR-spline paper).
        double& gamma()       {return gamma_;}
  const double& gamma() const {return gamma_;}

  // Access the rational weight of this LRBSpline2D.
        double& weight()       {return weight_;}
  const double& weight() const {return weight_;}

  const bool rational() const {return rational_;}

  // Get the dimension of the LRBSpline2Ds codomain.
  // For rational cases the geometry dimension is 1 less.
  const int dimension() const {return coef_times_gamma_.dimension();}

  // Access the LRBSpline2D's knot vector in the given direction.  (The knot vectors
  // only contain incices to an external, shared vector of knot values).
  // To get the knotvector in the first direction (x-direction), 'd' should be XFIXED.
  // To get the knotvector in the second direction (y-direction), 'd' should be YFIXED.
  const std::vector<int>& kvec(Direction2D d) const {return (d==XFIXED) ? kvec_u_ : kvec_v_;}
        std::vector<int>& kvec(Direction2D d)       {return (d==XFIXED) ? kvec_u_ : kvec_v_;}

  // Get the polynomial degree of the spline.
  const int degree(Direction2D d) const {return (int)kvec(d).size() - 2;}  

  /// Get the index to the knot that defines the start (end) of the LRBSpline2D's support.
  // (The vector of the actual knot values is stored outside of the LRBSpline2D, as it 
  // is shared among many LRBSpline2Ds).
  const int suppMin(Direction2D d) const {return kvec(d).front();}
  const int suppMax(Direction2D d) const {return kvec(d).back();}

  /// Information about the domain covered by this B-spline
  double umin() const 
  { 
    return mesh_->kval(XFIXED, kvec_u_[0]);
  };
  double umax() const 
  { 
    return mesh_->kval(XFIXED, kvec_u_[kvec_u_.size()-1]);
  };
  double vmin() const 
  { 
    return mesh_->kval(YFIXED, kvec_v_[0]);    
  };
  double vmax() const 
  {
    return mesh_->kval(YFIXED, kvec_v_[kvec_v_.size()-1]);    
  };
  // Query whether the parameter point speficied by the knots indexed by 'u_ix' and 'v_ix' 
  // is covered by the support of this LRBSpline2D.  (NB: The vector of the actual knot 
  // values is stored outsde of the LRBSpline2D, since it is shared among many 
  // LRBSpline2Ds.
  bool coversCorner(int u_ix, int v_ix) const { 
    return 
      u_ix >= suppMin(XFIXED) &&  u_ix < suppMax(XFIXED) && 
      v_ix >= suppMin(YFIXED) &&  v_ix < suppMax(YFIXED);
  }

  Point getGrevilleParameter() const;

  // Operations related to the support of this B-spline
  bool overlaps(Element2D *el) const;
  bool addSupport(Element2D *el) ;
  void removeSupport(Element2D *el) ;
  std::vector<Element2D*>::iterator supportedElementBegin() ;
  std::vector<Element2D*>::iterator supportedElementEnd() ;
  std::vector<Element2D*> getExtendedSupport() ;
  std::vector<Element2D*> getMinimalExtendedSupport();
  std::vector<Element2D*> supportedElements()
    {
      return support_;
    }
  void setSupport(std::vector<Element2D*> elements)
  {
    support_ = elements;
  }

  int nmbSupportedElements() { return (int)support_.size(); };

  void setMesh(const Mesh2D* mesh)
  {
    mesh_ = mesh;
  }

  const Mesh2D* getMesh()
  {
    return mesh_;
  }

  void subtractKnotIdx(int u_del, int v_del);

  // -----------------
  // --- OPERATORS ---
  // -----------------

  // Operator defining a partial ordering of LRBSpline2Ds.
  bool operator<(const LRBSpline2D& rhs) const;

  // Equality operator
  bool operator==(const LRBSpline2D &rhs) const;

 private:

  Point coef_times_gamma_;
  double weight_; // For rational case.
  bool rational_; // @@sbr201301 Should this also be part of the LRSplineSurface? It seems
                  // best suited for this class since it is here we use the rational part.
  double gamma_; // normalizing weight to ensure partition of unity, c.f. Section 7 of paper
  std::vector<int> kvec_u_;
  std::vector<int> kvec_v_;
  std::vector<Element2D*> support_;  // Elements lying in the support of this LRB-spline
  const Mesh2D *mesh_; // Information about global knot vectors and multiplicities

}; // end class LRBSpline2D

 inline std::ostream& operator<<(std::ostream& os, const LRBSpline2D& b) {b.write(os); return os;}
 inline std::istream& operator>>(std::istream& is,       LRBSpline2D& b) {b.read (is); return is;}

}; // end namespace Go

#endif
