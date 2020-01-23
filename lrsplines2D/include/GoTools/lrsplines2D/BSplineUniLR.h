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

#ifndef BSPLINEUNILR_H
#define BSPLINEUNILR_H

#include <vector>
#include <assert.h>

#include "GoTools/lrsplines2D/MeshLR.h"
#include "GoTools/utils/checks.h"
#include "GoTools/utils/StreamUtils.h"
#include "GoTools/geometry/Streamable.h"

namespace Go
{


class BSplineUniLR : public Streamable
//==============================================================================
{
 public:

  // ---------------------------------------------------------
  // --- CONSTRUCTORS, READING, WRITING AND SWAP FUNCTIONS --- 
  // ---------------------------------------------------------

  /// Constructor to create an empty (invalid) BSplineUniLR
  BSplineUniLR() 
    { }; 

  template<typename Iterator>
    BSplineUniLR(int pardir, int deg, 
		 Iterator kvec_start, const MeshLR* mesh)
    : pardir_(pardir), kvec_(kvec_start, kvec_start + deg + 2),
    mesh_(mesh), count_(0)
    {}

  /// Copy constructor
  BSplineUniLR(const BSplineUniLR& rhs);

  /// Swap the contents of two BSplineUniLRs
  void swap(BSplineUniLR& rhs) 
  {
    kvec_.swap(rhs.kvec_);
    std::swap(pardir_, rhs.pardir_);
  }

  ~BSplineUniLR() 
    {  
      //std::cout << "Delete LRBSpline " << this << std::endl;
    }; 

  /// Write the BSplineUniLR to a stream
  virtual void write(std::ostream& os) const;
  
  /// Read the BSplineUniLR from a stream
  virtual void read(std::istream& is);

  // ---------------------------
  // --- EVALUATION FUNCTION ---
  // ---------------------------

  /// Evaluate value of basis function in the parameter par
  double evalBasisFunc(double par) const;

  /// Evaluate one specified derivate in the parameter par
  double evalBasisFunction(double par, int deriv = 0,
			   bool at_end = false) const;

  /// Evaluate value and a number of derivatives in the parameter par.
  /// Note that the function is tested only up to and including deriv=3.
  /// For higher order derivatives use evalBasisFunction
  void evalBasisFunctions(double par, int deriv, double der[],
			  bool at_end = false) const;

 
  // -----------------------
  // --- QUERY FUNCTIONS ---
  // -----------------------

  // Parameter direction
  int pardir() const
  {
    return pardir_;
  }

  // Access the BSplineUni's knot vector in the given direction.  (The knot vectors
  // only contain incices to an external, shared vector of knot values).
  const std::vector<int>& kvec() const 
  {return kvec_;}
  const std::vector<int>& kvec_const() const 
  {return kvec_;}
  std::vector<int>& kvec()       
    {return kvec_;}

  // Get the polynomial degree of the spline.
  const int degree() const 
	{return (int)kvec().size() - 2;}  

  /// Get the index to the knot that defines the start (end) of the LRBSpline2D's support.
  // (The vector of the actual knot values is stored outside of the LRBSpline2D, as it 
  // is shared among many LRBSpline2Ds).
  const int suppMin() const 
  {return kvec().front();}
  const int suppMax() const 
  {return kvec().back();}

  /// Information about the parameter interval covered by this B-spline
  double min() const 
  { 
    return mesh_->kval(pardir_, kvec_[0]);
  };
  double max() const 
  { 
    return mesh_->kval(pardir_, kvec_[kvec_.size()-1]);
  };

  double knotval(int kn) const
  {
    return mesh_->kval(pardir_, kn);
  }

  // Count multiplicity in the ends of the B-spline
  int endmult(bool atstart) const;

  // Query whether the parameter speficied by the knots indexed by 'ix'
  // is covered by the support of this BSplineUniLR.  
  bool coversPar(int ix) const { 
    return (ix >= suppMin() &&  ix < suppMax());
  }

  // Check if the knot indexed by ix is used in the B-spline description
  bool useKnot(int ix) const {
    std::vector<int>::const_iterator it = 
      std::find(kvec_.begin(), kvec_.end(), ix);
    return (it != kvec_.end());
  }

  double getGrevilleParameter() const;

  bool overlaps(double pmin, double pmax) const; 

  void setMesh(const MeshLR* mesh)
  {
    mesh_ = mesh;
  }

  const MeshLR* getMesh()
  {
    return mesh_;
  }

  void setPardir(int pardir)
  {
    pardir_ = pardir;
  }

  void subtractKnotIdx(int del);

  void reverseParameterDirection();

  // -----------------
  // --- OPERATORS ---
  // -----------------

  // Operator defining a partial ordering of LRBSpline2Ds.
  int operator<(const BSplineUniLR& rhs) const;

  // Equality operator
  bool operator==(const BSplineUniLR& rhs) const;


  // Increase instance counter
  void incrCount() const
  {
    count_++;
  }

  // Decrease counter
  void decrCount() const
  {
    count_--;
  }

  int getCount() const
  {
    return count_;
  }

    private:
  int pardir_;
  std::vector<int> kvec_;
  const MeshLR *mesh_; // Information about global knot vectors and multiplicities

  mutable int count_;

    }; // end class LRBSpline2D

  inline std::ostream& operator<<(std::ostream& os, const BSplineUniLR& b) 
  {b.write(os); return os;}
  inline std::istream& operator>>(std::istream& is, BSplineUniLR& b) 
  {b.read (is); return is;}

  inline int BSplineUniLR::operator<(const BSplineUniLR& rhs) const
  { return compare_seq(kvec_.begin(), kvec_.end(), rhs.kvec_.begin(), rhs.kvec_.end()); }

}; // end namespace Go

#endif
