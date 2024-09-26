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

#include "GoTools/lrsplines3D/Mesh3D.h"
#include "GoTools/lrsplines3D/Mesh3DUtils.h"
#include "GoTools/utils/checks.h"
#include "GoTools/utils/StreamUtils.h"
#include "GoTools/lrsplines3D/Mesh3DIterator.h"
#include "GoTools/lrsplines3D/IndexMesh3DIterator.h"

#include <vector>
#include <set>
#include <assert.h>
#include <algorithm>
#include <stdexcept>

#include <array>

using std::vector;
using std::pair;
using std::find_if;
using std::max_element;
using std::min_element;

namespace Go
{

namespace { // anonymous namespace

// This function just checks the integrity of the vector of meshrectangles.
// - Verifies that the entries in the vector is properly sorted.
// - Verifies that there are no trivial entries (describing intervals of zero
//   length.
// - Verifies that multiplicities are >= 0.
bool mrvec_is_correct(const vector<GPos2D>& vec);

};


// =============================================================================
Mesh3D::Mesh3D(std::istream& is) {read(is); }
// =============================================================================

// =============================================================================
void Mesh3D::write(std::ostream& os) const
// =============================================================================
{
  object_to_stream(os, knotvals_x_); // #_global_knots_u, knot_u_0, ..., knot_u_n
  object_to_stream(os, knotvals_y_); // #_global_knots_v, knot_v_0, ..., knot_v_n
  object_to_stream(os, knotvals_z_); // #_global_knots_v, knot_v_0, ..., knot_v_n
  object_to_stream(os, mrects_x_); // #_mesh_rect_x 
  object_to_stream(os, mrects_y_);
  object_to_stream(os, mrects_z_);
}

// =============================================================================
void Mesh3D::read(std::istream& is)
// =============================================================================
{
  Mesh3D tmp;
  object_from_stream(is, tmp.knotvals_x_);
  object_from_stream(is, tmp.knotvals_y_);
  object_from_stream(is, tmp.knotvals_z_);
  object_from_stream(is, tmp.mrects_x_);
  object_from_stream(is, tmp.mrects_y_);
  object_from_stream(is, tmp.mrects_z_);
  //tmp.consistency_check_();
  swap(tmp);
}

// =============================================================================
void Mesh3D::swap(Mesh3D& rhs)
// =============================================================================
{
  std::swap(knotvals_x_, rhs.knotvals_x_);
  std::swap(knotvals_y_, rhs.knotvals_y_);
  std::swap(knotvals_z_, rhs.knotvals_z_);
  std::swap(mrects_x_,   rhs.mrects_x_);
  std::swap(mrects_y_,   rhs.mrects_y_);
  std::swap(mrects_z_,   rhs.mrects_z_);
}

// =============================================================================
void Mesh3D::consistency_check_() const
// =============================================================================
{
  assert(mrects_x_.size() > 1);
  assert(mrects_y_.size() > 1);
  assert(mrects_z_.size() > 1);

  assert(strictly_increasing(knotvals_x_));
  assert(strictly_increasing(knotvals_y_));
  assert(strictly_increasing(knotvals_z_));

  for (auto mv = mrects_x_.begin(); mv != mrects_x_.end(); ++mv)
    assert(mrvec_is_correct(*mv));
  for (auto mv = mrects_y_.begin(); mv != mrects_y_.end(); ++mv) 
    assert(mrvec_is_correct(*mv));
  for (auto mv = mrects_z_.begin(); mv != mrects_z_.end(); ++mv) 
    assert(mrvec_is_correct(*mv));

  // for (const auto& mv : mrects_x_) assert(mrvec_is_correct(mv));
  // for (const auto& mv : mrects_y_) assert(mrvec_is_correct(mv));
}

// =============================================================================
Mesh3DIterator Mesh3D::begin() const
// =============================================================================
{ 
  return Mesh3DIterator(*this, 0, 0, 0);
}

// =============================================================================
Mesh3DIterator Mesh3D::end() const
// =============================================================================
{ 
  return Mesh3DIterator(*this, numDistinctKnots(XDIR), numDistinctKnots(YDIR), numDistinctKnots(ZDIR));
}

// =============================================================================
IndexMesh3DIterator Mesh3D::indexMeshBegin() const
// =============================================================================
{ 
  return IndexMesh3DIterator(*this);
}

// =============================================================================
IndexMesh3DIterator Mesh3D::indexMeshEnd() const
// =============================================================================
{ 
  return IndexMesh3DIterator(*this, -1, -1, -1, -1);
}

// =============================================================================
int Mesh3D::nu(Direction3D d, int ix,
               int start1, int end1, int start2, int end2) const
// =============================================================================
{
  if (end1 >= numDistinctKnots(next(d)) || end2 >= numDistinctKnots(prev(d)))
    return 0; // proposed meshrectangle surpasses grid
  
  if (!(end1 > start1) || !(end2 > start2))
    return 0; // we can now safely assume that end1 > start1 and end2 > start2

  // We pick a layer of mesh rectangles.
  //const auto& mr = select_meshvec_(d, ix);
  std::vector<GPos2D> mr;
  select_meshvec(d, ix, start1, end1, start2, end2, mr);

  // return value
  int nu = 0; 

  // @obar. This algorithm is robust but slow - probably a lot of optimization can be done.
  // We make a look up table 'mults' and compare the 'nu' value when intersections are detected.
  std::set<int> vec_next, vec_prev;
  vec_next.insert(0); // first index
  vec_next.insert(this->numDistinctKnots(next(d))-1); // last index
  vec_prev.insert(0); // first index
  vec_prev.insert(this->numDistinctKnots(prev(d))-1); // last index
  for (auto i = mr.begin(); i != mr.end(); ++i)
    {
      vec_next.insert(i->ll[0]);
      vec_next.insert(i->ur[0]);
      vec_prev.insert(i->ll[1]);
      vec_prev.insert(i->ur[1]);
    }

  // Moves mesh rectangle indicis into vectors as it seems to be small benefit
  // in traversing vectors compared to sets
  vector<int> vec_next2(vec_next.begin(), vec_next.end());
  vector<int> vec_prev2(vec_prev.begin(), vec_prev.end());

  // Create multiplicity look up table
  vector<int> mults((vec_next.size()-1)*(vec_prev.size()-1),0);
  for (auto i = mr.begin(); i != mr.end(); ++i)
    {
      vector<int>::iterator it1 = std::find(vec_next2.begin(), vec_next2.end(), i->ll[0]);
      vector<int>::iterator it2 = std::find(vec_next2.begin(), vec_next2.end(), i->ur[0]);
      vector<int>::iterator it3 = std::find(vec_prev2.begin(), vec_prev2.end(), i->ll[1]);
      vector<int>::iterator it4 = std::find(vec_prev2.begin(), vec_prev2.end(), i->ur[1]);
      int ix_start = std::distance(vec_next2.begin(), it1);
      int jx_start = std::distance(vec_prev2.begin(), it3);
      int ix_end = std::distance(vec_next2.begin(), it2);
      int jx_end = std::distance(vec_prev2.begin(), it4);
      
      // int ix_start = std::distance(vec_next.begin(),vec_next.find(i->ll[0]));
      // int jx_start = std::distance(vec_prev.begin(),vec_prev.find(i->ll[1]));
      // int ix_end = std::distance(vec_next.begin(),vec_next.find(i->ur[0]));
      // int jx_end = std::distance(vec_prev.begin(),vec_prev.find(i->ur[1]));
      for (int ix=ix_start; ix!=ix_end; ++ix)
        {
          for (int jx=jx_start; jx!=jx_end; ++jx)
            {
              nu = std::max(nu,i->mult);
              mults[(vec_prev.size()-1)*ix+jx] = i->mult;
            }
        }
    }

  // Look up multiplicities if within the table partitions intersect the input rectangle
  int n_ix=0;
  int p_jx=0;
  bool no_intersections = true;
  for (n_ix=0; n_ix<(int)vec_next2.size()-1; ++n_ix)
    {
      int rstart1 = vec_next2[n_ix];
      int rend1 = vec_next2[n_ix+1];
      for (p_jx=0; p_jx<(int)vec_prev2.size()-1; ++p_jx)
        {
  	  int rstart2 = vec_prev2[p_jx];
  	  int rend2 = vec_prev2[p_jx+1];
          if ( end1 <= rstart1 || start1 >= rend1 || end2 <= rstart2 || start2 >= rend2) continue;
          else
            {
              no_intersections = false;
              nu = std::min(nu,mults[(vec_prev2.size()-1)*n_ix+p_jx]);
            }
        }
    }
  // for (auto it=vec_next.begin(); it!=std::prev(vec_next.end()); ++it, ++n_ix)
  //   {
  //     p_jx=0;
  //     for (auto jt=vec_prev.begin(); jt!=std::prev(vec_prev.end()); ++jt, ++p_jx)
  //       {
  //         int rstart1=*it, rstart2=*jt, rend1=*(std::next(it)), rend2=*(std::next(jt));
  //         if ( end1 <= rstart1 || start1 >= rend1 || end2 <= rstart2 || start2 >= rend2) continue;
  //         else
  //           {
  //             no_intersections = false;
  //             nu = std::min(nu,mults[(vec_prev.size()-1)*n_ix+p_jx]);
  //           }
  //       }
  //   }

  return no_intersections ? 0 : nu;
}

#if 0
// =============================================================================
int Mesh3D::extent(Direction3D d, int ix, int start, int mult) const
// =============================================================================
{
  const auto& mr = select_meshvec_(d, ix);
  auto p1 = find_if(mr.begin(), mr.end(), [start](const GPos2D& g) {return g.ix >= start;});
  if (p1 == mr.end() || p1->ix > start) --p1; // this will never break as long as start >= 0 
  auto p2 = find_if(p1, mr.end(), [mult](const GPos2D& g) {return g.mult < mult;});
  if (p2 == p1) return 0;
  else if (p2 == mr.end()) return numDistinctKnots(flip(d)) - 1 - start;
  else return p2->ix - start; 
}
#endif

// =============================================================================
void Mesh3D::incrementMult(Direction3D d, int ix,
			   int start1, int end1,
			   int start2, int end2,
			   int mult)
// =============================================================================
{
  // The current implementation is dead simple but likely inefficient.  A more
  // efficient implementation could be done at a later stage if this proves to 
  // be a relevant bottleneck (but that might be unlikely...).

  int curr_mult_low = nu(d, ix, start1, start1+1, start2, start2+1);
  int curr_mult_high = curr_mult_low;
  for (int j = start2; j < end2; ++j)
    for (int i = start1; i < end1; ++i)
      {
        const int curr_mult = nu(d, ix, i, i+1, j, j+1);
        if (curr_mult < curr_mult_low)
          curr_mult_low = curr_mult;
        else if (curr_mult > curr_mult_high)
          curr_mult_high = curr_mult;
        //	setMult(d, ix, i, i+1, j, j+1, cur_mult + mult); // set the new multiplicity to the old one plus 'mult'
      }

  // @@sbr201304 When this is triggered we must update the code to support this type of structure.
  assert(curr_mult_low == curr_mult_high);

  bool refined =
    setMult(d, ix, start1, end1, start2, end2, curr_mult_low + mult); // set the new multiplicity to the old one plus 'mult'
}
#if 1
// =============================================================================
bool Mesh3D::setMult(Direction3D d, int ix,
                     int start1, int end1,
                     int start2, int end2,
                     int mult)
// =============================================================================
{
  assert((end1 > start1) && (end2 > start2) && (mult >= 0)); // this is the contract for calling the function

  // We fetch our set of rectangles in the ix slice.
  std::vector<GPos2D>& mr = select_meshvec_(d, ix); 

  // We do not store 0-multiplicities in the rectangle slices.
  if (mr.size() == 0 && mult == 0)
    return false;

  if (mr.size() == 0)
    {
      mr.push_back(GPos2D(start1, start2, end1, end2, mult));
      return true;
    }

  // Do a pre run to remove mesh rectangles completely contained in the new mesh rectangle,
  // update the size of the new mesh rectangle by adding previous mesh rectangles if the
  // composition defines a rectangle, the domain of the new mesh rectangle can be updated,
  // and check if the new mesh rectangle is completely contained in a previous one.
  // Also look cases where it is sufficient to split the mesh rectangle into two pieces.
  // The idea is to reduce the number of mesh rectangles as much as possible to improve
  // performance when deducing information from the mesh rectangles (Mesh3D::nu)
  size_t ki;
  for (ki=0; ki<mr.size();)
    {
      int mrstart1 = mr[ki].ll[0];
      int mrstart2 = mr[ki].ll[1];
      int mrend1 = mr[ki].ur[0];
      int mrend2 = mr[ki].ur[1];
      int mrmult = mr[ki].mult;
      
      if ( mrend1 < start1 || mrstart1 > end1 || mrend2 < start2 || mrstart2 > end2)
	{
	  ++ki;
	  continue;
	}
      if (start1 >= mrstart1 && end1 <= mrend1 && start2 >= mrstart2 &&
	  end2 <= mrend2 && mult <= mrmult)
	{
	  // New mesh rectangle contained in current. Dismiss new.
	  break;
	}
    if (start1 <= mrstart1 && end1 >= mrend1 && start2 <= mrstart2 &&
	end2 >= mrend2 && mrmult <= mult)
      {
	// Current mesh rectange contained in new. Remove
	mr.erase(mr.begin()+ki);
      }
    else if (mrstart1 == start1 && mrend1 == end1 && mrmult == mult)
      {
	// New mesh rectangle extends current in direction 2. Merge current mesh
	// rectangle into new
    	start2 = std::min(start2, mrstart2);
    	end2 = std::max(end2, mrend2);
    	mr.erase(mr.begin()+ki);
      }
    else if (mrstart2 == start2 && mrend2 == end2 && mrmult == mult)
      {
	// New mesh rectangle extends current in direction 1. Merge current mesh
	// rectangle into new
    	start1 = std::min(start1, mrstart1);
    	end1 = std::max(end1, mrend1);
    	mr.erase(mr.begin()+ki);
      }
    else if (mrstart1 <= start1 && mrend1 >= end1 && mrmult == mult &&
	     !(start2 <= mrstart2 && end2 >= mrend2))
      {
	// New mesh rectangle is covered by the current one in direction 1 and extends
	// the current in only one end in direction 2. Remove the union from the new
	// mesh rectangle
	if (start2 > mrstart2)
	  start2 = mrend2;
	else if (end2 < mrend2)
	  end2 = mrstart2;
	++ki;
      }
    else if (mrstart2 <= start2 && mrend2 >= end2 && mrmult == mult &&
	     !(start1 <= mrstart1 && end1 >= mrend1)) 
      {
	// New mesh rectangle is covered by the current one in direction 2 and extends
	// the current in only one end in direction 1. Remove the union from the new
	// mesh rectangle
	if (start1 > mrstart1)
	  start1 = mrend1;
	else if (end1 < mrend1)
	  end1 = mrstart1;
	++ki;
      }
    else if (start1 <= mrstart1 && end1 >= mrend1 && mrmult == mult &&
    	     (!(mrstart2 < start2 && mrend2 > end2)))
      {
    	// Current mesh rectangle is covered by the new in direction 1 and extends
    	// the new in only one end in direction 2. Remove the union from the current
    	// mesh rectangle
    	if (mrstart2 > start2)
    	  mr[ki] = GPos2D(mrstart1, end2, mrend1, mrend2, mr[ki].mult);
    	else if (mrend2 < end2)
    	  mr[ki] = GPos2D(mrstart1, mrstart2, mrend1, start2, mr[ki].mult);
    	++ki;
      }
    else if (start2 <= mrstart2 && end2 >= mrend2 && mrmult == mult &&
    	     (!(mrstart1 < start1 && mrend1 > end1)))
      {
    	// Current mesh rectangle is covered by the new in direction 2 and extends
    	// the new in only one end in direction 1. Remove the union from the current
    	// mesh rectangle
    	if (mrstart1 > start1)
    	  mr[ki] = GPos2D(end1, mrstart2, mrend1, mrend2, mr[ki].mult);
    	else if (mrend1 < end1)
    	  mr[ki] = GPos2D(mrstart1, mrstart2, start1, mrend2, mr[ki].mult);
    	++ki;
      }
    else
      ++ki;
    }

  if (ki < mr.size())
    return false;   // Finished
   
  // Search through mesh rectangles for those affected
  //auto p = mr.end();
  bool detected = false;

  // A temporary variable that will store the new mesh rectangles
  vector<GPos2D> tmp_gpos_vec;
  tmp_gpos_vec.reserve(9);// 9 should be enough
  GPos2D tmp_gpos;
  for (ki=0; ki<mr.size();)
    {
      int mrstart1 = mr[ki].ll[0];
      int mrstart2 = mr[ki].ll[1];
      int mrend1 = mr[ki].ur[0];
      int mrend2 = mr[ki].ur[1];
      
      // A temporary variable that will be used to overwrite this mesh rectangle at the end
      tmp_gpos = mr[ki];

      // If the statement below is true, then there is no intersection: continue
      // if ( mrend1 <= start1 || mrstart1 >= end1 || mrend2 <= start2 || mrstart2 >= end2)
      if ( mrend1 <= start1 || mrstart1 >= end1 || mrend2 <= start2 || mrstart2 >= end2)
	{
	  ++ki;
	  continue;
	}

    if (mrstart1 >= start1 && mrend1 <= end1 && mrstart2 >= start2 && mrend2 <= end2)
      {
	// Previous mesh rectange contained in new. Remove previous
	mr.erase(mr.begin()+ki);
      }
    else
      {
	// Otherwise there is an intersection and we need to deal with it.
	// First step is to split to make a unified grid containing the intersecting rectangles
	vector<int> sort1 = {mrstart1,mrend1,start1,end1};
	vector<int> sort2 = {mrstart2,mrend2,start2,end2};
	std::sort(sort1.begin(),sort1.end());
	std::sort(sort2.begin(),sort2.end());
	// Run through the grid, checking for intersections
	for (int ix=0; ix!=3; ++ix) {
	  if (sort1[ix]==sort1[ix+1])
	    continue;
	  for (int jx=0; jx!=3; ++jx) {
	    if (sort1[ix]==sort1[ix+1] || sort2[jx]==sort2[jx+1])
	      continue;
	    // If the statement below is true, then there is an intersection: insert
	    if (!(sort1[ix+1] <= start1 || sort1[ix] >= end1 || sort2[jx+1] <= start2 || sort2[jx] >= end2)) {
	      // Just to avoid erasing the intersected mesh rectangle
	      if (!detected) {
		tmp_gpos = GPos2D(sort1[ix], sort2[jx], sort1[ix+1], sort2[jx+1], mult);
		// subsequent mesh rectangles will be inserted at the end of the vector
		detected = true;
	      }
	      else tmp_gpos_vec.push_back(GPos2D(sort1[ix], sort2[jx], sort1[ix+1], sort2[jx+1], mult));
	      continue;
	    }
	    if (!(sort1[ix+1] <= mrstart1 || sort1[ix] >= mrend1 || sort2[jx+1] <= mrstart2 || sort2[jx] >= mrend2)) {
	      if (!detected) {
		// subsequent mesh rectangles will be inserted at the end of the vector
		tmp_gpos = GPos2D(sort1[ix], sort2[jx], sort1[ix+1], sort2[jx+1], mr[ki].mult);
		detected = true;
	      }
	      else tmp_gpos_vec.push_back(GPos2D(sort1[ix], sort2[jx], sort1[ix+1], sort2[jx+1], mr[ki].mult));
	      continue;
	    }
	  }
	}
	// We can safely change the current mesh rectangle after we have finished with it
	// Changing this mesh rectangle avoids us erasing it from the vector
	mr[ki] = tmp_gpos;
	++ki;
      }
    }
  if (detected) {
    mr.insert(mr.end(),tmp_gpos_vec.begin(),tmp_gpos_vec.end());
  }
  else {
    // In this case, the new mesh rectangle does not intersect any existing
    // mesh rectangles, so we keep the old mesh rectangles and insert a
    // new one at the end.
    //result.insert(result.end(), mr.begin(), mr.end());
    mr.push_back(GPos2D(start1, start2, end1, end2, mult));
  }

  // Check mesh rectangles
  for (ki=0; ki<mr.size(); ++ki)
    {
      int mrstart1 = mr[ki].ll[0];
      int mrstart2 = mr[ki].ll[1];
      int mrend1 = mr[ki].ur[0];
      int mrend2 = mr[ki].ur[1];
      if (mrstart1 >= mrend1 || mrstart2 >= mrend2)
	{
	  std::cout << "Illegal mesh rectangle: (" << mrstart1;
	  std::cout << ", " << mrstart2 << "), (" << mrend1;
	  std::cout << ", " << mrend2 << ")" << std::endl;
	}
    }
 
  // verify that end contract of this function is fulfilled
  //assert(mrvec_is_correct(mr));
  return true;
}
#endif
#if 0
// =============================================================================
void Mesh3D::setMult(Direction3D d, int ix,
                     int start1, int end1,
                     int start2, int end2,
                     int mult)
// =============================================================================
{
  assert((end1 > start1) && (end2 > start2) && (mult >= 0)); // this is the contract for calling the function

  // We fetch our set of rectangles in the ix slice.
  auto& mr = select_meshvec_(d, ix); // auto == std::vector<GPos2D>

  // We do not store 0-multiplicities in the rectangle slices.
  if (mr.size() == 0 && mult == 0)
    return;

  // Search through mesh rectangles for those affected
  //auto p = mr.end();
  bool detected = false;

  // A temporary variable that will store the new mesh rectangles
  vector<GPos2D> tmp_gpos_vec;
  tmp_gpos_vec.reserve(9);// 9 should be enough
  GPos2D tmp_gpos;
  for (auto mr_it=mr.begin(); mr_it!=mr.end(); ++mr_it) {
    // A temporary variable that will be used to overwrite this mesh rectangle at the end
    tmp_gpos = *mr_it;

    int mrstart1 = mr_it->ll[0];
    int mrstart2 = mr_it->ll[1];
    int mrend1 = mr_it->ur[0];
    int mrend2 = mr_it->ur[1];

    // If the statement below is true, then there is no intersection: continue
    // if ( mrend1 <= start1 || mrstart1 >= end1 || mrend2 <= start2 || mrstart2 >= end2)
    if ( mrend1 < start1 || mrstart1 > end1 || mrend2 < start2 || mrstart2 > end2)
      continue;

    // Special cases
    if (start1 >= mrstart1 && end1 <= mrend1 && start2 >= mrstart2 && end2 <= mrend2)
      {
	// New mesh rectangle contained in previous. Dismiss new.
	detected = true;
	break;
      }
    if (mrstart1 >= start1 && mrend1 <= end1 && mrstart2 >= start2 && mrend2 <= end2)
      {
	// Previous mesh rectange contained in new. Update previous
	(*mr_it) = GPos2D(start1, start2, end1, end2, mult);
	detected = true;
	continue;
      }
    if (mrstart1 == start1 && mrend1 == end1)
      {
	// Equal boundaries in first direction. Update previous mesh rectangle
	(*mr_it) = GPos2D(start1, std::min(start2,mrstart2), end1, std::max(end2,mrend2), mult);
	detected = true;
	continue;
      }
    if (mrstart2 == start2 && mrend2 == end2)
      {
	// Equal boundaries in first direction. Update previous mesh rectangle
	(*mr_it) = GPos2D(std::min(start1,mrstart1), start2, std::max(end1,mrend1), end2, mult);
	detected = true;
	continue;
      }
 	
    // Otherwise there is an intersection and we need to deal with it.
    // First step is to split to make a unified grid containing the intersecting rectangles
    vector<int> sort1 = {mrstart1,mrend1,start1,end1};
    vector<int> sort2 = {mrstart2,mrend2,start2,end2};
    std::sort(sort1.begin(),sort1.end());
    std::sort(sort2.begin(),sort2.end());
    // Run through the grid, checking for intersections
    for (int ix=0; ix!=3; ++ix) {
      if (sort1[ix]==sort1[ix+1])
	continue;
      for (int jx=0; jx!=3; ++jx) {
        if (sort1[ix]==sort1[ix+1] || sort2[jx]==sort2[jx+1])
	  continue;
        // If the statement below is true, then there is an intersection: insert
        if (!(sort1[ix+1] <= start1 || sort1[ix] >= end1 || sort2[jx+1] <= start2 || sort2[jx] >= end2)) {
          // Just to avoid erasing the intersected mesh rectangle
          if (!detected) {
            tmp_gpos = GPos2D(sort1[ix], sort2[jx], sort1[ix+1], sort2[jx+1], mult);
            // subsequent mesh rectangles will be inserted at the end of the vector
            detected = true;
          }
          else tmp_gpos_vec.push_back(GPos2D(sort1[ix], sort2[jx], sort1[ix+1], sort2[jx+1], mult));
          continue;
        }
        if (!(sort1[ix+1] <= mrstart1 || sort1[ix] >= mrend1 || sort2[jx+1] <= mrstart2 || sort2[jx] >= mrend2)) {
          if (!detected) {
            // subsequent mesh rectangles will be inserted at the end of the vector
            tmp_gpos = GPos2D(sort1[ix], sort2[jx], sort1[ix+1], sort2[jx+1], mr_it->mult);
            detected = true;
          }
          else tmp_gpos_vec.push_back(GPos2D(sort1[ix], sort2[jx], sort1[ix+1], sort2[jx+1], mr_it->mult));
          continue;
        }
      }
    }
    // We can safely change the current mesh rectangle after we have finished with it
    // Changing this mesh rectangle avoids us erasing it from the vector
    *mr_it = tmp_gpos;
  }

  if (detected) {
    mr.insert(mr.end(),tmp_gpos_vec.begin(),tmp_gpos_vec.end());
  }
  else {
    // In this case, the new mesh rectangle does not intersect any existing
    // mesh rectangles, so we keep the old mesh rectangles and insert a
    // new one at the end.
    //result.insert(result.end(), mr.begin(), mr.end());
    mr.push_back(GPos2D(start1, start2, end1, end2, mult));
  }

  // verify that end contract of this function is fulfilled
  //assert(mrvec_is_correct(mr));
}
#endif
// =============================================================================
int Mesh3D::insertRectangle(Direction3D d, double kval, int mult)
// =============================================================================
{
  // Mesh3D does not impose any specific tolerance, but assumes that there is one
  // in place at a higher level.  It thus accepts new knot values that are 
  // arbitrarily close to existing ones, but protests if the new knot value is 
  // _exactly_ the same (bit-wise) as an existing one, since this breaks with the
  // principle that the knot-vectors should be strictly increasing.
  vector<double>& kvec = (d == XDIR) ? knotvals_x_ : ((d == YDIR)? knotvals_y_ : knotvals_z_);

  auto p = find_if(kvec.begin(), kvec.end(), [kval](double d) {return d >= kval;});
  if (*p == kval) 
    THROW("Knotvalue already in vector.");

  const int ix = (p - kvec.begin());  // this is the index of the line to insert
  
  auto& target = (d == XDIR) ? mrects_x_ : ((d == YDIR) ? mrects_y_ : mrects_z_);
  auto& other1  = (d == XDIR) ? mrects_y_ : ((d == YDIR) ? mrects_z_ : mrects_x_);
  auto& other2  = (d == XDIR) ? mrects_z_ : ((d == YDIR) ? mrects_x_ : mrects_y_);

  // @@sbr201304 Do we need to set a rectangle with multiplicity 0? Why not
  // let the vector be empty?

  // inserting new rectangle
  int last1 = other1.size() - 1;
  int last2 = other2.size() - 1;
  vector<GPos2D> ref_rect;
  if (mult > 0)
    ref_rect.push_back(GPos2D(0, 0, last1, last2, mult));
  target.insert(target.begin() + ix, ref_rect);//vector<GPos2D> (1, GPos2D(0, 0, last1, last2, mult)));
  kvec.insert(p, kval);
  
  // Adjust indexes in the other1 & other2.
  for (auto gvec_it = other1.begin(); gvec_it != other1.end(); ++gvec_it)
    for (auto g_it = gvec_it->begin(); g_it != gvec_it->end(); ++g_it)
      {
	if (g_it->ll[1] >= ix)
	  ++(g_it->ll[1]);
	if (g_it->ur[1] >= ix)
	  ++(g_it->ur[1]);
      }

  for (auto gvec_it = other2.begin(); gvec_it != other2.end(); ++gvec_it)
    for (auto g_it = gvec_it->begin(); g_it != gvec_it->end(); ++g_it)
      {
	if (g_it->ll[0] >= ix)
	  ++(g_it->ll[0]);
	if (g_it->ur[0] >= ix)
	  ++(g_it->ur[0]);
      }
  
  return ix;
}


// =============================================================================
void Mesh3D::setParameterDomain(double u1, double u2, double v1, double v2, double w1, double w2)
// =============================================================================
{
  double umin = minParam(XDIR);
  double umax = maxParam(XDIR);
  double vmin = minParam(YDIR);
  double vmax = maxParam(YDIR);
  double wmin = minParam(ZDIR);
  double wmax = maxParam(ZDIR);

  knotvals_x_[0] = u1;
  knotvals_x_[knotvals_x_.size()-1] = u2;
  knotvals_y_[0] = v1;
  knotvals_y_[knotvals_y_.size()-1] = v2;
  knotvals_z_[0] = w1;
  knotvals_z_[knotvals_z_.size()-1] = w2;

  // (u_old - umin)/(umax - umin) = (u_new - u1)/(u2 - u1).
  double u_quot = (u2 - u1)/(umax - umin);
  for (size_t ki = 1; ki < knotvals_x_.size() - 1; ++ki)
    knotvals_x_[ki] = (knotvals_x_[ki] - umin)*u_quot + u1;

  double v_quot = (v2 - v1)/(vmax - vmin);
  for (size_t ki = 1; ki < knotvals_y_.size() - 1; ++ki)
    knotvals_y_[ki] = (knotvals_y_[ki] - vmin)*v_quot + v1;

  double w_quot = (w2 - w1)/(wmax - wmin);
  for (size_t ki = 1; ki < knotvals_z_.size() - 1; ++ki)
    knotvals_z_[ki] = (knotvals_z_[ki] - wmin)*w_quot + w1;
}


// =============================================================================
void Mesh3D::swapParameterDirection(int pardir1, int pardir2)
// =============================================================================
{
  MESSAGE("LRSplineVolume::swapParameterDirection(): Not implemented yet.");
#if 0
  std::swap(knotvals_x_, knotvals_y_);
  std::swap(mrects_x_, mrects_y_);
#endif
}

// =============================================================================
void Mesh3D::reverseParameterDirection(int pardir)
// =============================================================================
{
  MESSAGE("LRSplineVolume::reverseParameterDirection(): Not implemented yet.");

#if 0  // We need the unique knots.
  // As well as the mrects_.
  if (dir_is_u)
    {
      MESSAGE("Reversing dir u.");
      vector<double> knotvals_x(knotvals_x_.size());
      double xmin = knotvals_x_.front();
      double xmax = knotvals_x_.back();
      for (size_t ki = 0; ki < knotvals_x_.size(); ++ki)
	knotvals_x[ki] = xmin + xmax - knotvals_x_[ki];
      knotvals_x_ = knotvals_x;

      // We must update the indices for the opposite direction.
      std::reverse(knotvals_x_.begin(), knotvals_x_.end());
      std::reverse(mrects_x_.begin(), mrects_x_.end());
      // We must also update the mrects_y_ indices.
      int last_ind = mrects_x_.size() - 1;
      for (auto iter = mrects_y_.begin(); iter != mrects_y_.end(); ++iter)
	{
	  for (auto iter2 = iter->begin(); iter2 != iter->end(); ++iter2)
	    {
//	      GPos2D = *iter2;
//	      int from_ind = iter2->ix;
	      int to_ind = (iter2 + 1 != iter->end()) ? iter2[1].ix : last_ind;
	      int new_from_ind = last_ind - to_ind;
	      iter2->ix = new_from_ind;
	    }

	  std::reverse(iter->begin(), iter->end());
	}
      
    }
  else
    {
      MESSAGE("Reversing dir v.");
      vector<double> knotvals_y(knotvals_y_.size());
      double ymin = knotvals_y_.front();
      double ymax = knotvals_y_.back();
      for (size_t ki = 0; ki < knotvals_y_.size(); ++ki)
	knotvals_y[ki] = ymin + ymax - knotvals_y_[ki];
      knotvals_y_ = knotvals_y;

      // We must update the indices for the opposite direction.
      std::reverse(knotvals_y_.begin(), knotvals_y_.end());
      std::reverse(mrects_y_.begin(), mrects_y_.end());

      // We must also update the mrects_y_ indices.
      int last_ind = mrects_y_.size() - 1;
      for (auto iter = mrects_x_.begin(); iter != mrects_x_.end(); ++iter)
	{
	  for (auto iter2 = iter->begin(); iter2 != iter->end(); ++iter2)
	    {
//	      GPos2D = *iter2;
//	      int from_ind = iter2->ix;
	      int to_ind = (iter2 + 1 != iter->end()) ? iter2[1].ix : last_ind;
	      int new_from_ind = last_ind - to_ind;
	      iter2->ix = new_from_ind;
	    }

	  std::reverse(iter->begin(), iter->end());
	}
    }
#endif
}


// =============================================================================
vector<double> Mesh3D::getKnots(Direction3D d, int ix1, int ix2, bool right) const
// =============================================================================
{
  vector<int> knot_idx;
  int beg = firstMeshVecIx(d);
  int end = lastMeshVecIx(d);
  
  int orto1_min, orto1_max, orto2_min, orto2_max;
  if ((right || ix1 == firstMeshVecIx(next(d))) &&
      (!(ix1 == lastMeshVecIx(next(d)))))
    {
      orto1_min = ix1;
      orto1_max = ix1 + 1;
    }
  else
    {
      orto1_min = ix1 - 1;
      orto1_max = ix1;
    }
  if ((right || ix2 == firstMeshVecIx(prev(d))) &&
      (!(ix2 == lastMeshVecIx(prev(d)))))
    {
      orto2_min = ix2;
      orto2_max = ix2 + 1;
    }
  else
    {
      orto2_min = ix2 - 1;
      orto2_max = ix2;
    }
  
  for (int pos = beg; pos <= end; ++pos) 
    {
      knot_idx.insert(knot_idx.end(), nu(d, pos, orto1_min, orto1_max, orto2_min, orto2_max), pos); 
      // 0 multiplicities will not be inserted
    }
  
  vector<double> knots(knot_idx.size());
  for (size_t k1=0; k1<knot_idx.size(); ++k1)
    knots[k1] = kval(d, knot_idx[k1]);

  return knots;
}

// =============================================================================
vector<pair<Corner2D, Corner2D> > Mesh3D::segments(Direction3D dir, int ix, int threshold) const
//vector<pair<int, int> > Mesh3D::segments(Direction3D d, int ix, int threshold) const
// =============================================================================
{
  const auto& mvec = select_meshvec_(dir, ix);
  vector<pair<Corner2D, Corner2D> > result;
  const int BLANK = -1; // use this flag to indicate uninitialized value
  int start1 = BLANK;
  int start2 = BLANK;
  for (auto i = mvec.begin(); i != mvec.end(); ++i)
    if (i->mult >= threshold) {
      result.emplace_back(pair<Corner2D, Corner2D>(i->ll, i->ur));
    }

  return result;
}

// =============================================================================
vector<pair<vector<GPos2D>, int> >
Mesh3D::overlapRects(Direction3D dir, int start1, int end1,
		     int start2, int end2) const

// =============================================================================
{
  int num = numDistinctKnots(dir);
  vector<pair<vector<GPos2D>,int> > result(num);  
  for (int ix=0; ix<num; ++ix)
    {
      vector<GPos2D>  curr_res;  
      const auto& mvec = select_meshvec_(dir, ix);
      for (auto i = mvec.begin(); i != mvec.end(); ++i)
	{
	  if (i->ll[0] >= end1 || i->ll[1] >= end2 ||
	      i->ur[0] <= start1 || i->ur[1] <= start2)
	    continue;
	  curr_res.emplace_back(GPos2D(i->ll[0], i->ll[1],
				       i->ur[0], i->ur[1], i->mult));
 	}
      result[ix] = std::make_pair(curr_res, ix);
    }
  return result;
}

// =============================================================================
int Mesh3D::largestInnerMult(Direction3D d) const
// =============================================================================
{
  int maxmult = 0;
  int nknots = numDistinctKnots(d);
  for (int ix=1; ix<nknots-1; ++ix)
    {
      const auto& mvec = select_meshvec_(d, ix);
      if (mvec.size() > 0)
	{
	  int currmult = max_element(mvec.begin(), 
				     mvec.end(), 
				     [](const GPos2D& a, const GPos2D& b) 
				     {return a.mult < b.mult;} )->mult;
	  maxmult = std::max(maxmult, currmult);
	}
    }
  return maxmult;
}

// =============================================================================
int Mesh3D::largestMultInLine(Direction3D d, int ix) const
// =============================================================================
{
  const auto& mvec = select_meshvec_(d, ix);
  assert( ! mvec.empty());
  return max_element(mvec.begin(), 
		     mvec.end(), 
		     [](const GPos2D& a, const GPos2D& b) {return a.mult < b.mult;} )->mult;
}

// =============================================================================
int Mesh3D::minMultInLine(Direction3D d, int ix) const
// =============================================================================
{
  const auto& mvec = select_meshvec_(d, ix);
  assert( ! mvec.empty());
  return min_element(mvec.begin(), 
		     mvec.end(), 
		     [](const GPos2D& a, const GPos2D& b) {return a.mult < b.mult;} )->mult;
}

// =============================================================================
  int Mesh3D::knotIntervalFuzzy(Direction3D d, double& par, double eps) const
// =============================================================================
{
  int ix = Mesh3DUtils::last_nonlarger_knotvalue_ix(*this, d, par);
  const double *st = knotsBegin(d);
  int nmb = numDistinctKnots(d); 
  if (par - st[ix] < eps)
    par = st[ix];
  else if (ix < nmb-1 && st[ix+1] - par < eps)
    {
      par = st[ix+1];
      while (ix < nmb-1 && st[ix] == st[ix+1])
	++ix;
    }

  return ix;
}

// =============================================================================
  int Mesh3D::getKnotIdx(Direction3D d, const double& par, double eps) const
// =============================================================================
{
  int ix = Mesh3DUtils::last_nonlarger_knotvalue_ix(*this, d, par);
  const double *st = (d == XDIR) ? &knotvals_x_[0] :
    ((d == YDIR) ? &knotvals_y_[0] : &knotvals_z_[0]);
  int nmb = (d == XDIR) ? (int)knotvals_x_.size() :
    ((d == YDIR) ? (int)knotvals_y_.size() : (int)knotvals_z_.size());
  if (par - st[ix] < eps)
    return ix;
  else if (ix < nmb-1 && st[ix+1] - par < eps)
    return ix+1;
  else
    return -1;
}

// =============================================================================
shared_ptr<Mesh3D> Mesh3D::subMesh(int ix1, int ix2,
				   int iy1, int iy2,
				   int iz1, int iz2) const
// =============================================================================
  {
    shared_ptr<Mesh3D> submesh(new Mesh3D());

    // Copy knot vectors
    submesh->knotvals_x_.insert(submesh->knotvals_x_.end(), 
				knotvals_x_.begin()+ix1, knotvals_x_.begin()+ix2+1);
    submesh->knotvals_y_.insert(submesh->knotvals_y_.end(), 
			       knotvals_y_.begin()+iy1, knotvals_y_.begin()+iy2+1);
    submesh->knotvals_z_.insert(submesh->knotvals_z_.end(), 
			       knotvals_z_.begin()+iz1, knotvals_z_.begin()+iz2+1);

    // Transfer knot multiplicity information. Note that the domain in all
    // parameter directions is involved for all parameter directions
    // 1. parameter direction
    int ki;
    for (ki=ix1; ki<=ix2; ++ki)
      {
	vector<GPos2D> x_orth;
	vector<GPos2D> mr = select_meshvec_(XDIR, ki);
	size_t kj;

	// Find first multiplicity instance. The instance may be connected to
	// a previous knot in the 2. or 3. parameter direction
	for (kj=0; kj<mr.size(); ++kj)
	  {
	    if (mr[kj].ur[0] <= iy1 || mr[kj].ll[0] >= iy2 ||
		mr[kj].ur[1] <= iz1 || mr[kj].ll[1] >= iz2)
	      continue;
	    x_orth.push_back(GPos2D(std::max(mr[kj].ll[0],iy1)-iy1,
				    std::max(mr[kj].ll[1],iz1)-iz1,
				    std::min(mr[kj].ur[0],iy2)-iy1,
				    std::min(mr[kj].ur[1],iz2)-iz1,
				    mr[kj].mult));
	    // if (mr[kj].ll[0] >= iy1 && mr[kj].ur[0] <= iy2 &&
	    // 	mr[kj].ll[1] >= iz1 && mr[kj].ur[1] <= iz2)
	    //   {	// We start indexing at (y1, z1)
	    // 	x_orth.push_back(GPos2D(mr[kj].ll[0]-iy1, mr[kj].ll[1]-iz1,
	    // 				mr[kj].ur[0]-iy1, mr[kj].ur[1]-iz1,
	    // 				mr[kj].mult));
	    //   }
	  }
	submesh->mrects_x_.push_back(x_orth);
      }

    // 2. parameter direction
    for (ki=iy1; ki<=iy2; ++ki)
      {
	vector<GPos2D> y_orth;
	vector<GPos2D> mr = select_meshvec_(YDIR, ki);

	size_t kj;

	// Find first multiplicity instance. The instance may be connected to
	// a previous knot in the 1. or 3. parameter direction
	for (kj=0; kj<mr.size(); ++kj)
	  {
	    if (mr[kj].ur[0] <= iz1 || mr[kj].ll[0] >= iz2 ||
		mr[kj].ur[1] <= ix1 || mr[kj].ll[1] >= ix2)
	      continue;
	    y_orth.push_back(GPos2D(std::max(mr[kj].ll[0],iz1)-iz1,
				    std::max(mr[kj].ll[1],ix1)-ix1,
				    std::min(mr[kj].ur[0],iz2)-iz1,
				    std::min(mr[kj].ur[1],ix2)-ix1,
				    mr[kj].mult));
	  }
	submesh->mrects_y_.push_back(y_orth);
      }

    // 3. parameter direction
    for (ki=iz1; ki<=iz2; ++ki)
      {
	vector<GPos2D> z_orth;
	vector<GPos2D> mr = select_meshvec_(ZDIR, ki);

	size_t kj;
	// Find first multiplicity instance. The instance may be connected to
	// a previous knot in the 1. or 2. parameter direction
	for (kj=0; kj<mr.size(); ++kj)
	  {
	    if (mr[kj].ur[0] <= ix1 || mr[kj].ll[0] >= ix2 ||
		mr[kj].ur[1] <= iy1 || mr[kj].ll[1] >= iy2)
	      continue;
	    z_orth.push_back(GPos2D(std::max(mr[kj].ll[0],ix1)-ix1,
				    std::max(mr[kj].ll[1],iy1)-iy1,
				    std::min(mr[kj].ur[0],ix2)-ix1,
				    std::min(mr[kj].ur[1],iy2)-iy1,
				    mr[kj].mult));
	  }
	submesh->mrects_z_.push_back(z_orth);
      }

    return submesh;
  }


// =============================================================================
namespace { // anonymous namespace
// =============================================================================

// =============================================================================
bool mrvec_is_correct(const vector<GPos2D>& v)
// =============================================================================
{
  GPos2D prev_gpos = v.front(); // this would always be smaller than first index, which is zero
  for (auto m = v.begin() + 1; m != v.end(); ++m) {
#if 0
    if      (m->ix <= prev_ix) return false; // should be _strictly_ increasing.
    else if (m->mult < 0)      return false; // multiplicities should never be negative
#endif
    // The second index should never decline.
    if (prev_gpos.ll[1] > m->ll[1])
      return false;
    // If equal second index, the first index should increase.
    else if ((prev_gpos.ll[1] == m->ll[1]) && (prev_gpos.ll[0]>= m->ll[0]))
      return false;
    else if (m->mult < 0)
      return false; // multiplicities should never be negative
    prev_gpos = *m;
  }
  return true;
}
  
}; // end anonymous namespace

}; // end namespace Go

