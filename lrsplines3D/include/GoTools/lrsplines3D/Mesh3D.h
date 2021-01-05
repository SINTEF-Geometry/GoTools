//===========================================================================
//                                                                           
// File: Mesh3D.h                                                            
//                                                                           
// Created: Mon Feb 25 11:07:09 2013                                         
//                                                                           
// Author: Sverre Briseid <Sverre.Briseid@sintef.no>
//                                                                           
// Revision: $Id$
//                                                                           
// Description:
//                                                                           
//===========================================================================

#ifndef _MESH3D_H
#define _MESH3D_H


#include <iostream>
#include <vector>
#include <array>
#include <assert.h>
#include "GoTools/geometry/Streamable.h"
#include "GoTools/utils/Array.h"
#include "GoTools/lrsplines2D/MeshLR.h"
#include "GoTools/lrsplines3D/Direction3D.h"
#include "GoTools/lrsplines3D/Mesh3DIterator.h"
#include "GoTools/lrsplines3D/IndexMesh3DIterator.h"

namespace Go
{

// // simple structure used for compact encoding of mesh topology.  
// struct GPos { 
//   // due to Visual Studio 2010 not supporting initializer lists, we have
//   // to make an explicit constructor here.
//   GPos(int i, int m) : ix(i), mult(m) {}
//   GPos() : ix(-1), mult(-1) {}

//   int ix; 
//   int mult;
// };  

  typedef Array<int, 2> Corner2D;

  // The domain is given by the lower left and upper right indices, as
  // well as the multiplicity of the local refinement.
  // @@sbr201302 For the 3D case the encoding is not as compact as the 2D case.
  // We use a vector of pairs denoting lower left and upper right of domains.
  // We may want to compress this structure using the same
  // approach as in the 2D case, defining the ll corner only and
  // letting the ur corner be deduced.
  struct GPos2D {
    GPos2D(int start1, int start2, int end1, int end2, int m)
      {
	ll[0] = start1;
	ll[1] = start2;
	ur[0] = end1;
	ur[1] = end2;
	mult = m;
      }
    GPos2D()
      : mult(-1)
      {}

    bool operator < (const GPos2D& other) const
      {
	if (ll[1] == other.ll[1])
	  return (ll[0] < other.ll[0]);
	else
	  return (ll[1] < other.ll[1]);
      }

    Corner2D ll;
    Corner2D ur;
    int mult;
  };

// =============================================================================
class Mesh3D : public MeshLR
// =============================================================================
{
public:
  // ---------------------------------------------------------
  // --- CONSTRUCTORS, READING, WRITING AND SWAP FUNCTIONS ---
  // ---------------------------------------------------------
  Mesh3D() {}; // empty mesh
  Mesh3D(std::istream& is); // read mesh from stream

  // Construct a full, 'tensor product' mesh, based on two knotvectors, expressed
  // in the provided ranges [kx_start, kx_end] and [ky_start, ky_end].  Multiplicities > 1
  // are allowed, and expressed by repeated values in the ranges.
  template<typename Iterator> 
  Mesh3D(Iterator kx_start, Iterator kx_end,
	 Iterator ky_start, Iterator ky_end,
	 Iterator kz_start, Iterator kz_end);
 
  // Construct a full, 'tensor product' mesh, based on two knotvectors, expressed
  // in the provided 1-D arrays 'xknots' and 'yknots'.  The arrays must have begin()
  // and end() member methods.  Multiplicities > 1 are allowed, and expressed by repeated
  // values.
  template<typename Array>
  Mesh3D(const Array& xknots,
	 const Array& yknots,
	 const Array& zknots);
  
  // Read the mesh from a stream
  virtual void read(std::istream& is);        

  // Write the mesh to a stream
  virtual void write(std::ostream& os) const; 

  // Swap two meshes
  void swap(Mesh3D& rhs);             

  // -----------------------
  // --- QUERY FUNCTIONS ---
  // -----------------------

  // The 'nu'-operator (c.f. definition in LR-spline paper).  Expresses, for a given consecutive
  // set of meshrectangles, the lowest multiplicity found therein (the lowest possible is 
  // 0, which can be interpreted as 'no meshrectangle').  The set of meshrectangles to be
  // examined is defined by the parameters 'd', 'ix', 'start' and 'end', and the resulting
  // lowest multiplicity is given as the return value.
  // The parameters mean:
  // d     - direction of the consecutive set of meshrectangles (are they lying on a row [YDIR] 
  //         or on a column [XDIR]).
  // ix    - the index of the column [XDIR] or row [YDIR] on which the meshrectangles are lying
  // start1 - the row (or column) index of the first meshrectangle of the consecutive set
  // end1   - the one-past-end index of the last meshrectangle of the consecutive set
  // start2 - the row (or column) index of the first meshrectangle of the consecutive set
  // end2   - the one-past-end index of the last meshrectangle of the consecutive set
  // If d==XDIR, start1 corresponds to YDIR. If d==YDIR, start1 corresponds to ZDIR. I.e. cyclic.
  int nu(Direction3D d, int ix,
	 int start1, int end1,
	 int start2, int end2) const;

  // Get the number of distinct knot valuess in a given direction (rows: YDIR, columns: XDIR).
  // Note that this is the number of _distinct_ knots, so multiplicities are not taken into
  // account.
  int numDistinctKnots(Direction3D d) const;
   int numDistinctKnots(int pardir) const;
 
  // Return the knot value for the knot with index 'ix' along direction 'd'
  // (rows: YDIR, columns: XDIR).
  double kval(Direction3D d, int ix) const;
  double kval(int pardir, int ix) const;

  // Get the lowest knot value (i.e. the first knot value in the knot vector), along a given
  // direction.
  double minParam(Direction3D d) const;

  // Get the highest knot value (i.e. the last knot value in the knot vector), along a given
  // direction.
  double maxParam(Direction3D d) const;

  // Get a pointer to the start of the knot vector in the given direction.
  const double* const knotsBegin(Direction3D d) const;
  virtual const double* const knotsBegin(int pardir) const;
  
  // Get a pointer to the one-past-end of the knot vector in the given direction.
  const double* const knotsEnd  (Direction3D d) const;
  virtual const double* const knotsEnd  (int pardir) const;

  std::vector<double> allKnots(Direction3D d) const;
  
  // Fetch the knot vector of the curve corresponding to a given row or column
  // Multiplicity is included
  // d  - determine whether to examine a row (YDIR) or a column (XDIR)
  // ix1 - index of row/column from which to fetch the knot vector
  // ix2 - index of row/column from which to fetch the knot vector
  std::vector<double> getKnots(Direction3D d, int ix1, int ix2, bool right=true) const;

#if 0
  // Determine the length of the longest k-meshrectangle with multiplicity 
  // (at least) 'mult' and with starting point at 'start'.
  // @@sbr201303 Does not seem to be in use, so no need to translate from 2D case.
  int extent(Direction3D d, int ix, int start, int mult) const;
#endif

  // Find the largest multiplicity of any of the meshrectangles excluding 
  // boundary
  // d  - direction (XDIR, YDIR or ZDIR)
  int largestInnerMult(Direction3D d) const; 

  // Find the largest multiplicity of any of the meshrectangles on a given row or column.
  // d  - determine whether to look at a row (YDIR) or column (XDIR)
  // ix - index of the row/column to examine.
  int largestMultInLine(Direction3D d, int ix) const; 

  // Find the minimum multiplicity of any of the meshrectangles on a given row or column.
  // d  - determine whether to look at a row (YDIR) or column (XDIR)
  // ix - index of the row/column to examine.
  int minMultInLine(Direction3D d, int ix) const; 

  // Fetch index of knot interval and modify parameter value if it is
  // very close to an existing knot (distance less than eps)
  int knotIntervalFuzzy(Direction3D d, double& par, double eps) const;

  // Fetch the index of a knot. The function returns -1 if no knot can be
  // found within an epsilon interval
  int getKnotIdx(Direction3D d, const double& par, double eps) const;

  // For a given row (or column) find all consecutive segments of meshrectangles with multiplicities
  // greater than or equal to a given threshold. Each found segment is represented as an integer pair,
  // representing the start index of the first meshrectangle in the segment and the one-past-end index
  // of the last meshrectangle in the segment.
  // d  - determine whether to examine a row (YDIR) or a column (XDIR)
  // ix - index of row/column to examine
  // threshold - the multiplicity threshold.  Default is one, which will give all segments of
  //             consecutive meshrectangles (considering those with multiplicity '0' to be nonexistent)
#if 0
  std::vector<std::pair<int, int> > segments(Direction3D dir, int ix, int threshold = 1) const;
#else
  std::vector<std::pair<Corner2D, Corner2D> > segments(Direction3D dir, int ix, int threshold = 1) const;
#endif

  std::vector<std::pair<std::vector<GPos2D>, int> >
    overlapRects(Direction3D dir, int start1,int end1, int start2, int end2) const;

  // Returns Mesh3DIterators referring to the first (begin) and
  // one-past-end (end) elements of the mesh. Mesh3DIterators can be
  // used to loop over elements in a mesh.
  Mesh3DIterator begin() const;
  Mesh3DIterator end() const;  

  // Returns IndexMesh3DIterators referring to the first
  // (indexMeshBegin) and to one-past-end (indexMeshEnd) element of
  // the index mesh (which is basically the mesh when counting the
  // knots with multiplicity). The IndexMesh3DIterator may be used for
  // looping over index elements. This is useful for check of linear
  // independence (where it must be examined which basis functions have
  // support on which index elements).
  IndexMesh3DIterator indexMeshBegin() const;
  IndexMesh3DIterator indexMeshEnd() const;

  // Index of first mesh rectangle in multiplicity vector
  int firstMeshVecIx(Direction3D d) const;
  // Index of last mesh rectangle in multiplicity vector
  int lastMeshVecIx(Direction3D d) const;
  
  // Return the specified vector of mrects, which provides detailed info on meshrectangles and
  // their multiplicities along a specified line
  // d  - direction of the specified line
  // ix - index of row/column to examine
  const std::vector<GPos2D>& mrects(Direction3D d, int ix) const;
  
  // Fetch sub mesh (using references) to the current mesh. Note that
  // the multiplicity of the start and end knots is as those for the
  // current mesh. It is no guarantee for knot multiplicity equal to
  // the order
  shared_ptr<Mesh3D> subMesh(int ix1, int ix2,
			     int iy1, int iy2,
			     int iz1, int iz2) const;

  // ----------------------
  // --- EDIT FUNCTIONS --- 
  // ----------------------

  // set the multiplicity of a consecutive set of meshrectangles.  
  // The consecutive set is specified by:
  // d  - the direction (a row: YDIR, a column: XDIR)
  // ix - the index of the row/column of the meshrectangles
  // start2 - the index to the start of the first consecutive meshrectangle along the line
  // end2   - the index to the one-past-end of the last consecutive meshrectangle along the line
  // start3 - the index to the start of the first consecutive meshrectangle along the line
  // end3   - the index to the one-past-end of the last consecutive meshrectangle along the line
    bool setMult(Direction3D d,
		 int ix,
		 int start1, int end1,
		 int start2, int end2,
		 int mult);

  // increment multiplicity of a consecutive set of meshrectangles by 'mult'
  // The consecutive set is specified by:
  // d  - the direction (a row: YDIR, a column: XDIR)
  // ix - the index of the row/column of the meshrectangles
  // start - the index to the start of the first consecutive meshrectangle along the line
  // end   - the index to the one-past-end of the last consecutive meshrectangle along the line
  void incrementMult(Direction3D d, int ix,
		     int start1, int end1,
		     int start2, int end2,
		     int mult);

  // Insert a line with X (or Y) fixed at 'kval', and with the
  // indicated multiplicity.  NB, 'kval' should be different from
  // any knot values already in the vector (otherwise it would not be
  // a new line).
  // Returns the index of the newly inserted line.
  // d    - specify whether the line should be parallel to y-axis (XDIR) or parallel to 
  //        x-axis (YDIR).  (A line parallel to the y-axis will yield a new knot in the 
  //        x-knotvector and vice versa).
  // kval - the knot value (i.e. parameter value in the 'fixed' direction).  Should be 
  //        different from any value already in the mesh in the given direction.
  // mult - the multiplicity of the meshrectangles on the new line.  (They will all have 
  //        the same multiplicity after insertion, but this can be changed with the 'setMult()' 
  //        and 'incrementMult()' member functions).
  int insertRectangle (Direction3D d, double kval, int mult = 0);

  // Change the parameter domain for the mesh.
  void setParameterDomain(double u1, double u2,
			  double v1, double v2,
			  double w1, double w2);

  // Using XDIR == 0, YDIR == 1, ZDIR == 2.
  void swapParameterDirection(int pardir1, int pardir2);

  // Using XDIR == 0, YDIR == 1, ZDIR == 2.
  void reverseParameterDirection(int pardir);

 private:

  // --------------------
  // --- PRIVATE DATA --- 
  // --------------------

  // store the values of the knots in x and y direction (only unique values - multiplicities
  // are accounted for in the grid itself.
  std::vector<double> knotvals_x_;
  std::vector<double> knotvals_y_;
  std::vector<double> knotvals_z_;

  std::vector<std::vector<GPos2D> > mrects_x_; // meshrectangles with x constant (|| with yz-plane)
                                               // (For 3D LRSplines this means rectangular segments.)
                                               // mrects_x_.size() == knotvals_x_.size().
                                               // The mesh rectangle consists of a number of rectangles
                                               // defined by their lower left and upper right indices, together with
                                               // the multiplicity of the rectangle.
                                               // If multiplicity alters, a new GPos2D follows.
                                               // @@sbr201304 We now only store GPos2D with a mult > 0, i.e.
                                               // we use an empty vector to indicate the mult in the parameter is 0.
  std::vector<std::vector<GPos2D> > mrects_y_; // meshrectangles with y constant (|| with xz-plane)
  std::vector<std::vector<GPos2D> > mrects_z_; // meshrectangles with z constant (|| with xy-plane)

  // -----------------------
  // --- PRIVATE METHODS ---
  // -----------------------

  // initializer function - used by all constructors
  template<typename Iterator>
  void init_(Iterator kx_start, Iterator kx_end,
	     Iterator ky_start, Iterator ky_end,
	     Iterator kz_start, Iterator kz_end);

  void consistency_check_() const ;
  const std::vector<GPos2D>& select_meshvec_(Direction3D d, int ix) const;
  std::vector<GPos2D>& select_meshvec_(Direction3D d, int ix);
  void select_meshvec(Direction3D d, int ix,
		      int start1, int end1,
		      int start2, int end2,
		      std::vector<GPos2D>& result) const;

  // Utility function that converts a knotvector (represented by the range [kvec_start, kvec_end) )
  // into a compact form consisting of a vector of all the unique values (return value of the 
  // function), and a vector containing all the corresponding multiplicities (argument 'm').
  template<typename Iterator>
  static std::vector<double> compactify_knotvec_(Iterator kvec_start, Iterator kvec_end, std::vector<int>& m);
  
}; // end class Mesh3D

// =============================================================================
template<typename Iterator>
Mesh3D::Mesh3D(Iterator kx_start, Iterator kx_end,
	       Iterator ky_start, Iterator ky_end,
	       Iterator kz_start, Iterator kz_end)
// =============================================================================
{
  init_(kx_start, kx_end, ky_start, ky_end, kz_start, kz_end);
}

// =============================================================================
template<typename Array>
Mesh3D::Mesh3D(const Array& xknots,
	       const Array& yknots,
	       const Array& zknots)
// =============================================================================
{
  init_(xknots.begin(), xknots.end(),
	yknots.begin(), yknots.end(),
	zknots.begin(), zknots.end());
}

// =============================================================================
template<typename Iterator>
void Mesh3D::init_(Iterator kx_start, Iterator kx_end,
		   Iterator ky_start, Iterator ky_end,
		   Iterator kz_start, Iterator kz_end)
// =============================================================================
{
  // saving the knotvals (NB: multiplicities should not be saved here, so sequence
  // should be strictly increasing).
  std::vector<int> mult_x, mult_y, mult_z; // will be used to store multiplicities
  knotvals_x_ = compactify_knotvec_(kx_start, kx_end, mult_x); // sequence of unique knot values returned,
  knotvals_y_ = compactify_knotvec_(ky_start, ky_end, mult_y); // while multiplicities go to 'mult_x'/'y'.
  knotvals_z_ = compactify_knotvec_(kz_start, kz_end, mult_z); // while multiplicities go to 'mult_x'/'y'.

  // We fill all the rectangles with the full domain since the inital knots are global.
  mrects_x_ = std::vector<std::vector<GPos2D> >(knotvals_x_.size(), std::vector<GPos2D>());
  for (int ki = 0; ki < knotvals_x_.size(); ++ki)
    if (mult_x[ki] > 0)
      mrects_x_[ki].push_back(GPos2D(0, 0, knotvals_y_.size() - 1, knotvals_z_.size() - 1, mult_x[ki]));
  mrects_y_ = std::vector<std::vector<GPos2D> >(knotvals_y_.size(), std::vector<GPos2D>());
  for (int ki = 0; ki < knotvals_y_.size(); ++ki)
    if (mult_y[ki] > 0)
      mrects_y_[ki].push_back(GPos2D(0, 0, knotvals_z_.size() - 1, knotvals_x_.size() - 1, mult_y[ki]));
  mrects_z_ = std::vector<std::vector<GPos2D> >(knotvals_z_.size(), std::vector<GPos2D>());//(1, GPos(0, 1)));
  for (int ki = 0; ki < knotvals_z_.size(); ++ki)
    if (mult_z[ki] > 0)
      mrects_z_[ki].push_back(GPos2D(0, 0, knotvals_x_.size() - 1, knotvals_y_.size() - 1, mult_z[ki]));

  // here we also checks that the grid is nonempty, which permits us to do the adjustments
  // further below.
  consistency_check_();

  // setting correct knot multiplicites
  for (int i = 0; i != numDistinctKnots(XDIR); ++i) 
    setMult(XDIR, i,
	    0, numDistinctKnots(YDIR) - 1,
	    0, numDistinctKnots(ZDIR) - 1,
	    mult_x[i]);
  for (int i = 0; i != numDistinctKnots(YDIR); ++i) 
    setMult(YDIR, i,
	    0, numDistinctKnots(ZDIR) - 1,
	    0, numDistinctKnots(XDIR) - 1,
	    mult_y[i]);
  for (int i = 0; i != numDistinctKnots(ZDIR); ++i) 
    setMult(ZDIR, i,
	    0, numDistinctKnots(XDIR) - 1,
	    0, numDistinctKnots(YDIR) - 1,
	    mult_z[i]);

}

// =============================================================================
template<typename Iterator>
std::vector<double> Mesh3D::compactify_knotvec_(Iterator kvec_start, Iterator kvec_end, std::vector<int>& mult) 
// =============================================================================
{
  // simple and inefficient implementation.  Optimization should not matter 
  // much here anyway, so the priority is on clear code.

  // Establishing the 'compact' knotvector without any explicit multiplicities
  std::vector<double> result;
  unique_copy(kvec_start, kvec_end, std::back_inserter(result));

  // counting multiplicities and keeping track of them in the 'mult' vector
  mult.clear();  
  for (auto i = result.begin(); i != result.end(); ++i) 
    mult.push_back(count(kvec_start, kvec_end, *i));

  return result;
}

// =============================================================================
inline double Mesh3D::minParam(Direction3D d) const
// =============================================================================
{
  return (d == XDIR) ? knotvals_x_.front() : (d == YDIR) ? knotvals_y_.front() : knotvals_z_.front();
}

// =============================================================================
inline double Mesh3D::maxParam(Direction3D d) const
// =============================================================================
{
  return (d == XDIR) ? knotvals_x_.back() : (d == YDIR) ? knotvals_y_.back() : knotvals_z_.back();
}

// =============================================================================
inline int Mesh3D::numDistinctKnots(Direction3D d) const 
// =============================================================================
{
  return (d == XDIR) ? (int)knotvals_x_.size() : (d == YDIR) ? (int)knotvals_y_.size() : (int)knotvals_z_.size();
};

// =============================================================================
inline int Mesh3D::numDistinctKnots(int pardir) const 
// =============================================================================
{
  return (pardir == 1) ? (int)knotvals_x_.size() : (pardir == 2) ? (int)knotvals_y_.size() : (int)knotvals_z_.size();
};

inline int Mesh3D::firstMeshVecIx(Direction3D d) const 
// =============================================================================
{
  return 0;
};

inline int Mesh3D::lastMeshVecIx(Direction3D d) const 
// =============================================================================
{
  return (d == XDIR) ? (int)mrects_x_.size()-1 : (d == YDIR) ? (int)mrects_y_.size()-1 : (int)mrects_z_.size()-1;
};

// =============================================================================
inline double Mesh3D::kval(Direction3D d, int ix) const
// =============================================================================
{
  return (d == XDIR) ? knotvals_x_[ix] : (d == YDIR) ? knotvals_y_[ix] : knotvals_z_[ix];
}

// =============================================================================
inline double Mesh3D::kval(int pardir, int ix) const
// =============================================================================
{
  return (pardir == 1) ? knotvals_x_[ix] : (pardir == 2) ? knotvals_y_[ix] : knotvals_z_[ix];
}

// =============================================================================
inline const std::vector<GPos2D>& Mesh3D::select_meshvec_(Direction3D d, int ix) const
// =============================================================================
{
  assert( (d == XDIR) || (d == YDIR) || (d == ZDIR));
  return (d == XDIR) ? mrects_x_[ix] : ((d == YDIR) ? mrects_y_[ix] : mrects_z_[ix]);
}

// =============================================================================
inline std::vector<GPos2D>& Mesh3D::select_meshvec_(Direction3D d, int ix)
// =============================================================================
{
  assert( (d == XDIR) || (d == YDIR) || (d == ZDIR));
  return (d == XDIR) ? mrects_x_[ix] : (d == YDIR) ? mrects_y_[ix] : mrects_z_[ix];
}

// =============================================================================
 inline void Mesh3D::select_meshvec(Direction3D d, int ix,
				    int start1, int end1,
				    int start2, int end2,
				    std::vector<GPos2D>& result) const
 // =============================================================================
{
  assert( (d == XDIR) || (d == YDIR) || (d == ZDIR));
  if (d == XDIR)
    {
      for (size_t ki=0; ki<mrects_x_[ix].size(); ++ki)
	{
	  if (mrects_x_[ix][ki].ll[0] > end1)
	    continue;
	  if (mrects_x_[ix][ki].ll[1] > end2)
	    continue;
	  if (mrects_x_[ix][ki].ur[0] < start1)
	    continue;
	  if (mrects_x_[ix][ki].ur[1] < start2)
	    continue;
	  result.push_back(mrects_x_[ix][ki]);
	}
    }
  else  if (d == YDIR)
    {
      for (size_t ki=0; ki<mrects_y_[ix].size(); ++ki)
	{
	  if (mrects_y_[ix][ki].ll[0] > end1)
	    continue;
	  if (mrects_y_[ix][ki].ll[1] > end2)
	    continue;
	  if (mrects_y_[ix][ki].ur[0] < start1)
	    continue;
	  if (mrects_y_[ix][ki].ur[1] < start2)
	    continue;
	  result.push_back(mrects_y_[ix][ki]);
	}
    }
  else 
    {
     for (size_t ki=0; ki<mrects_z_[ix].size(); ++ki)
	{
	  if (mrects_z_[ix][ki].ll[0] > end1)
	    continue;
	  if (mrects_z_[ix][ki].ll[1] > end2)
	    continue;
	  if (mrects_z_[ix][ki].ur[0] < start1)
	    continue;
	  if (mrects_z_[ix][ki].ur[1] < start2)
	    continue;
	  result.push_back(mrects_z_[ix][ki]);
	}
    }
}

//==============================================================================
inline const double* const Mesh3D::knotsBegin(Direction3D d) const
//==============================================================================
{
  return (d == XDIR) ? &knotvals_x_[0] : (d == YDIR) ? &knotvals_y_[0] : &knotvals_z_[0];
}

//==============================================================================
inline const double* const Mesh3D::knotsBegin(int pardir) const
//==============================================================================
{
  return (pardir == 1) ? &knotvals_x_[0] : (pardir == 2) ? &knotvals_y_[0] : &knotvals_z_[0];
}

//==============================================================================
inline const double* const Mesh3D::knotsEnd(Direction3D d) const
//==============================================================================
{
  return knotsBegin(d) + numDistinctKnots(d);
}

//==============================================================================
inline const double* const Mesh3D::knotsEnd(int pardir) const
//==============================================================================
{
  return knotsBegin(pardir) + numDistinctKnots(pardir);
}

// =============================================================================
inline std::vector<double> Mesh3D::allKnots(Direction3D d) const
// =============================================================================
{
  return (d == XDIR) ? knotvals_x_ : (d == YDIR) ? knotvals_y_ : knotvals_z_;
}

// =============================================================================
inline const std::vector<GPos2D>& Mesh3D::mrects(Direction3D d, int ix) const
// =============================================================================
{
  return select_meshvec_(d, ix);
}

// =============================================================================
inline Direction3D prev(Direction3D d) 
// =============================================================================
{ 
    return (d == XDIR) ? ZDIR : ((d == YDIR) ? XDIR : YDIR);
} 

// =============================================================================
inline Direction3D next(Direction3D d) 
// =============================================================================
{ 
    return (d == XDIR) ? YDIR : ((d == YDIR) ? ZDIR : XDIR);
} 

// defining streaming operators
inline std::ostream& operator<<(std::ostream& os, const GPos2D& g)  { return os << g.ll << " " << g.ur << " " << g.mult << " ";}
inline std::istream& operator>>(std::istream& is, GPos2D& g)        { return is >> g.ll >> g.ur >> g.mult;}
inline std::ostream& operator<<(std::ostream& os, const Mesh3D& m){ m.write(os); return os;}
inline std::istream& operator>>(std::istream& is, Mesh3D& m)      { m.read(is); return is;}

}; // end namespace Go

#endif // _MESH3D_H

