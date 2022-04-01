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

#include "GoTools/lrsplines2D/Mesh2D.h"
#include "GoTools/lrsplines2D/Mesh2DUtils.h"
#include "GoTools/utils/checks.h"
#include "GoTools/utils/StreamUtils.h"
#include "GoTools/lrsplines2D/Mesh2DIterator.h"
#include "GoTools/lrsplines2D/IndexMesh2DIterator.h"

#include <vector>
#include <assert.h>
#include <algorithm>
#include <stdexcept>
//#include <iostream>
//#include <sstream>
//#include <string>
#include <array>
//#include <cctype>

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
bool mrvec_is_correct(const vector<GPos>& vec);

};


// =============================================================================
Mesh2D::Mesh2D(std::istream& is) {read(is); }
// =============================================================================

// =============================================================================
  Mesh2D::Mesh2D(const std::vector<double>& xknots,
		 const std::vector<double>& yknots,
		 const std::vector<std::vector<int> >& mrvecx,
		 const std::vector<std::vector<int> >& mrvecy)
    : knotvals_x_(xknots), knotvals_y_(yknots)
// =============================================================================
  {
    mrects_x_.resize(xknots.size());
    mrects_y_.resize(yknots.size());

    // Create matrix for knot multiplicity corresponding to xknots
    vector<vector<int> > xmults(xknots.size());
    for (size_t ki=0; ki<xmults.size(); ++ki)
      xmults[ki].resize(mrvecy.size(), 0);
    for (size_t kj=0; kj<mrvecy.size(); ++kj)
      {
	vector<int> mult;
	vector<int> ixy = compactify_ixvec_(mrvecy[kj].begin(), mrvecy[kj].end(),
					    mult);
	for (size_t kr=0; kr<ixy.size(); ++kr)
	  xmults[ixy[kr]][kj] = mult[kr];
      }
    
    // Create matrix for knot multiplicity corresponding to yknots
    vector<vector<int> > ymults(yknots.size());
    for (size_t ki=0; ki<ymults.size(); ++ki)
      ymults[ki].resize(mrvecx.size(), 0);
    for (size_t kj=0; kj<mrvecx.size(); ++kj)
      {
	vector<int> mult;
	vector<int> ixx = compactify_ixvec_(mrvecx[kj].begin(), mrvecx[kj].end(),
					    mult);
	for (size_t kr=0; kr<ixx.size(); ++kr)
	  ymults[ixx[kr]][kj] = mult[kr];
      }

    // Collect mesh rectangles
    for (size_t ki=0; ki<xmults.size(); ++ki)
      {
	GPos pos0(0, xmults[ki][0]);
	mrects_x_[ki].push_back(pos0);
	for (size_t kj=1; kj<xmults[ki].size(); ++kj)
	  {
	    GPos pos1((int)kj, xmults[ki][kj]);
	    if (pos1.mult != pos0.mult)
	      {
		mrects_x_[ki].push_back(pos1);
		pos0 = pos1;
	      }
	  }
      }
	      
    for (size_t ki=0; ki<ymults.size(); ++ki)
      {
	GPos pos0(0, ymults[ki][0]);
	mrects_y_[ki].push_back(pos0);
	for (size_t kj=1; kj<ymults[ki].size(); ++kj)
	  {
	    GPos pos1((int)kj, ymults[ki][kj]);
	    if (pos1.mult != pos0.mult)
	      {
		mrects_y_[ki].push_back(pos1);
		pos0 = pos1;
	      }
	  }
      }
	      
	    

    int stop_break = 1;
  }
  
// =============================================================================
void Mesh2D::write(std::ostream& os) const
// =============================================================================
{
  object_to_stream(os, knotvals_x_); // #_global_knots_u, knot_u_0, ..., knot_u_n
  object_to_stream(os, knotvals_y_); // #_global_knots_v, knot_v_0, ..., knot_v_n
  object_to_stream(os, mrects_x_); // #_mesh_rect_x 
  object_to_stream(os, mrects_y_);
}

// =============================================================================
void Mesh2D::read(std::istream& is)
// =============================================================================
{
  Mesh2D tmp;
  object_from_stream(is, tmp.knotvals_x_);
  object_from_stream(is, tmp.knotvals_y_);
  object_from_stream(is, tmp.mrects_x_);
  object_from_stream(is, tmp.mrects_y_);
  tmp.consistency_check_();
  swap(tmp);
}

// =============================================================================
void Mesh2D::swap(Mesh2D& rhs)
// =============================================================================
{
  std::swap(knotvals_x_, rhs.knotvals_x_);
  std::swap(knotvals_y_, rhs.knotvals_y_);
  std::swap(mrects_x_,   rhs.mrects_x_);
  std::swap(mrects_y_,   rhs.mrects_y_);
}

// =============================================================================
void Mesh2D::consistency_check_() const
// =============================================================================
{
    // This function is called from the (public) read() function. We can not trust the user to fulfill the
    // contract for calling the function. Hence it is not a good idea to use assert which, when triggered,
    // results in a core dump.
    if (mrects_x_.size() < 2)
        THROW("We need at least 2 mesh values in the x direction!");
    if (mrects_y_.size() < 2)
        THROW("We need at least 2 mesh values in the y direction!");

    if (!strictly_increasing(knotvals_x_))
        THROW("The knotvals_x_ should be strictly increasing!");
    if (!strictly_increasing(knotvals_y_))
        THROW("The knotvals_y_ should be strictly increasing!");

    for (auto mv = mrects_x_.begin(); mv != mrects_x_.end(); ++mv)
        if (!mrvec_is_correct(*mv))
            THROW("Found illegal mrvec!");
    for (auto mv = mrects_y_.begin(); mv != mrects_y_.end(); ++mv) 
        if (!mrvec_is_correct(*mv))
            THROW("Found illegal mrvec!");

  // for (const auto& mv : mrects_x_) assert(mrvec_is_correct(mv));
  // for (const auto& mv : mrects_y_) assert(mrvec_is_correct(mv));
}

// =============================================================================
Mesh2DIterator Mesh2D::begin() const
// =============================================================================
{ 
  return Mesh2DIterator(*this, 0, 0);
}

// =============================================================================
Mesh2DIterator Mesh2D::end() const
// =============================================================================
{ 
  return Mesh2DIterator(*this, numDistinctKnots(XFIXED), numDistinctKnots(YFIXED));
}

// =============================================================================
IndexMesh2DIterator Mesh2D::indexMeshBegin() const
// =============================================================================
{ 
  return IndexMesh2DIterator(*this);
}

// =============================================================================
IndexMesh2DIterator Mesh2D::indexMeshEnd() const
// =============================================================================
{ 
  return IndexMesh2DIterator(*this, -1, -1, -1, -1);
}

// =============================================================================
int Mesh2D::nu(Direction2D d, int ix, int start, int end) const
// =============================================================================
{
  if (end >= numDistinctKnots(flip(d))) return 0; // proposed meshrectangle surpasses grid
  const auto& mr = select_meshvec_(d, ix);
  if (!(end > start)) return 0; // we can now safely assume that end > start

  int result = mr[0].mult;
  for (auto i = mr.begin(); i != mr.end(); ++i) 
    if      (i->ix <= start) result = i->mult;
    else if (i->ix >= end)   break; // finished
    else if (i->mult == 0)   return 0; // gap encountered - nu is zero
    else                     result = std::min(result, i->mult);
  
  return result;
}

// =============================================================================
int Mesh2D::extent(Direction2D d, int ix, int start, int mult) const
// =============================================================================
{
  const auto& mr = select_meshvec_(d, ix);
  auto p1 = find_if(mr.begin(), mr.end(), [start](const GPos& g) {return g.ix >= start;});
  if (p1 == mr.end() || p1->ix > start) --p1; // this will never break as long as start >= 0 
  auto p2 = find_if(p1, mr.end(), [mult](const GPos& g) {return g.mult < mult;});
  if (p2 == p1) return 0;
  else if (p2 == mr.end()) return numDistinctKnots(flip(d)) - 1 - start;
  else return p2->ix - start; 
}

// =============================================================================
  void Mesh2D::incrementMult(Direction2D d, int ix, int start, int end, int mult)
// =============================================================================
{
  // The current implementation is dead simple but likely inefficient.  A more
  // efficient implementation could be done at a later stage if this proves to 
  // be a relevant bottleneck (but that might be unlikely...).

  for (int i = start; i < end; ++i) {
    const int cur_mult = nu(d, ix, i, i+1);
    // set the new multiplicity to the old one plus 'mult'
    bool refined = 
      setMult(d, ix, i, i+1, cur_mult + mult);
  }
}

// =============================================================================
  bool Mesh2D::setMult(Direction2D d, int ix, int start, int end, int mult)
// =============================================================================
{
  if ((end <= start) || (mult < 0))
    THROW("Negative multiplicity or end less than start");
  //assert((end > start) && (mult >= 0)); // this is the contract for calling the function
  auto& mr = select_meshvec_(d, ix); // auto == std::vector<GPos>
  const int last_pos = numDistinctKnots(flip(d)) - 1;
  vector<GPos> result;
  result.reserve(mr.size() + 1);

  // keeping preceeding GPos (not affected by the change)
  auto p = find_if(mr.begin(), mr.end(), [start](GPos& g) {return g.ix >= start;});
  result.insert(result.end(), mr.begin(), p);
  
  // inserting the new GPos
  if (result.empty() || result.back().mult != mult) 
    result.push_back(GPos(start, mult));
  if (end < last_pos) {
    int tmp = nu(d, ix, end, end+1);
    if (tmp != mult)
      result.push_back(GPos(end, tmp));
  }

  // inserting succeeding GPoses (not affected by the change)
  p = find_if(p, mr.end(), [end](GPos& g) {return g.ix > end;});
  result.insert(result.end(), p, mr.end());

  // if (result.size() <= mr.size())
  //   return false;
  mr.swap(result);

  // verify that end contract of this function is fulfilled
  assert(mrvec_is_correct(mr));
  return
    true;
}

// =============================================================================
  int Mesh2D::insertLine(Direction2D d, double kval, int mult)
// =============================================================================
{
  // Mesh2D does not impose any specific tolerance, but assumes that there is one
  // in place at a higher level.  It thus accepts new knot values that are 
  // arbitrarily close to existing ones, but protests if the new knot value is 
  // _exactly_ the same (bit-wise) as an existing one, since this breaks with the
  // principle that the knot-vectors should be strictly increasing.
  vector<double>& kvec = (d == XFIXED) ? knotvals_x_ : knotvals_y_;

  auto p = find_if(kvec.begin(), kvec.end(), [kval](double d) {return d >= kval;});
  if (*p == kval) 
    THROW("Knotvalue already in vector.");

  const int ix = (p - kvec.begin());  // this is the index of the line to insert
  
  auto& target = (d == XFIXED) ? mrects_x_ : mrects_y_;
  auto& other  = (d == XFIXED) ? mrects_y_ : mrects_x_;

  // inserting new line
  target.insert(target.begin() + ix, vector<GPos> (1, GPos(0, mult)));
  kvec.insert(p, kval);
  
  // adjust indexes in the other 
  for (auto gvec_it = other.begin(); gvec_it != other.end(); ++gvec_it)
    for (auto g_it = gvec_it->begin(); g_it != gvec_it->end(); ++g_it)
      if (g_it->ix >= ix) ++(g_it->ix);
  
  return ix;
}


// =============================================================================
void Mesh2D::setParameterDomain(double u1, double u2, double v1, double v2)
// =============================================================================
{
  double umin = minParam(XFIXED);
  double umax = maxParam(XFIXED);
  double vmin = minParam(YFIXED);
  double vmax = maxParam(YFIXED);

  knotvals_x_[0] = u1;
  knotvals_x_[knotvals_x_.size()-1] = u2;
  knotvals_y_[0] = v1;
  knotvals_y_[knotvals_y_.size()-1] = v2;

  // (u_old - umin)/(umax - umin) = (u_new - u1)/(u2 - u1).
  double u_quot = (u2 - u1)/(umax - umin);
  for (size_t ki = 1; ki < knotvals_x_.size() - 1; ++ki)
    knotvals_x_[ki] = (knotvals_x_[ki] - umin)*u_quot + u1;

  double v_quot = (v2 - v1)/(vmax - vmin);
  for (size_t ki = 1; ki < knotvals_y_.size() - 1; ++ki)
    knotvals_y_[ki] = (knotvals_y_[ki] - vmin)*v_quot + v1;
}


// =============================================================================
void Mesh2D::swapParameterDirection()
// =============================================================================
{
  std::swap(knotvals_x_, knotvals_y_);
  std::swap(mrects_x_, mrects_y_);
}

// =============================================================================
void Mesh2D::reverseParameterDirection(bool dir_is_u)
// =============================================================================
{
  // We need the unique knots.
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
//	      GPos = *iter2;
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
//	      GPos = *iter2;
//	      int from_ind = iter2->ix;
	      int to_ind = (iter2 + 1 != iter->end()) ? iter2[1].ix : last_ind;
	      int new_from_ind = last_ind - to_ind;
	      iter2->ix = new_from_ind;
	    }

	  std::reverse(iter->begin(), iter->end());
	}
    }
}

// =============================================================================
int Mesh2D::removeUnusedLines(Direction2D d)
// =============================================================================
{
  // Identify unused parameter lines
  vector<int> empty_ix;
  for (int i = 0; i != numDistinctKnots(d); ++i)
    if (largestMultInLine(d, i) == 0)
      empty_ix.push_back(i);

  vector<double>& kvals       =   (d==XFIXED) ? knotvals_x_ : knotvals_y_;
  vector<vector<GPos>>&mrects =   (d==XFIXED) ? mrects_x_   : mrects_y_;
  vector<vector<GPos>>&mr_other = (d==XFIXED) ? mrects_y_   : mrects_x_;
  
  // Remove unused parameter lines
  int num_removed = int(empty_ix.size());
  for (int j = num_removed - 1; j >= 0; --j) {
    kvals.erase(kvals.begin() + empty_ix[j]);
    mrects.erase(mrects.begin() + empty_ix[j]);
    for (auto it = mr_other.begin(); it != mr_other.end(); ++it)
      for (auto  it2 = it->begin(); it2 != it->end(); ++it2)
       if (it2->ix >= empty_ix[j])
         --(it2->ix);
  }
  return num_removed;
}
  

// =============================================================================
  vector<double> Mesh2D::getKnots(Direction2D d, int ix, bool right) const
// =============================================================================
{
  vector<int> knot_idx;
  int beg = firstMeshVecIx(d);
  int end = lastMeshVecIx(d);
  
  int orto_min, orto_max;
  if ((right || ix == firstMeshVecIx(flip(d))) &&
      (!(ix == lastMeshVecIx(flip(d)))))
    {
      orto_min = ix;
      orto_max = ix + 1;
    }
  else
    {
      orto_min = ix - 1;
      orto_max = ix;
    }
  
  for (int pos = beg; pos <= end; ++pos) 
    {
      knot_idx.insert(knot_idx.end(), nu(d, pos, orto_min, orto_max), pos); 
      // 0 multiplicities will not be inserted
    }
  
  vector<double> knots(knot_idx.size());
  for (size_t k1=0; k1<knot_idx.size(); ++k1)
    knots[k1] = kval(d, knot_idx[k1]);

  return knots;
}

// =============================================================================
vector<pair<int, int> > Mesh2D::segments(Direction2D d, int ix, int threshold) const
// =============================================================================
{
  const auto& mvec = select_meshvec_(d, ix);
  vector<pair<int, int> > result;
  const int BLANK = -1; // use this flag to indicate uninitialized value
  int start = BLANK;
  for (auto i = mvec.begin(); i != mvec.end(); ++i)
    if (i->mult >= threshold && start == BLANK) start = i->ix;
    else if (i->mult < threshold && start != BLANK) {
      result.emplace_back(pair<int, int>(start, i->ix));
      start = BLANK;
    }

  if (start != BLANK) 
    result.emplace_back(pair<int, int>(start, numDistinctKnots(flip(d)) - 1));

  return result;
}

// =============================================================================
vector<pair<int, int> > Mesh2D::zeroSegments(Direction2D d, int ix) const
// =============================================================================
{
  const auto& mvec = select_meshvec_(d, ix);
  vector<pair<int, int> > result;
  const int BLANK = -1; // use this flag to indicate uninitialized value
  int start = BLANK;
  for (auto i = mvec.begin(); i != mvec.end(); ++i)
    if (i->mult == 0 && start == BLANK) start = i->ix;
    else if (i->mult > 0 && start != BLANK) {
      result.emplace_back(pair<int, int>(start, i->ix));
      start = BLANK;
    }

  if (start != BLANK) 
    result.emplace_back(pair<int, int>(start, numDistinctKnots(flip(d)) - 1));

  return result;
}

// =============================================================================
int Mesh2D::largestInnerMult(Direction2D d) const
// =============================================================================
{
  int maxmult = 0;
  int nknots = numDistinctKnots(d);
  for (int ix=1; ix<nknots-1; ++ix)
    {
      const auto& mvec = select_meshvec_(d, ix);
      assert( ! mvec.empty());
      int currmult = max_element(mvec.begin(), 
				 mvec.end(), 
				 [](const GPos& a, const GPos& b) 
				 {return a.mult < b.mult;} )->mult;
      maxmult = std::max(maxmult, currmult);
    }
  return maxmult;
}

// =============================================================================
int Mesh2D::largestMultInLine(Direction2D d, int ix) const
// =============================================================================
{
  const auto& mvec = select_meshvec_(d, ix);
  assert( ! mvec.empty());
  return max_element(mvec.begin(), 
		     mvec.end(), 
		     [](const GPos& a, const GPos& b) {return a.mult < b.mult;} )->mult;
}

// =============================================================================
int Mesh2D::minMultInLine(Direction2D d, int ix) const
// =============================================================================
{
  const auto& mvec = select_meshvec_(d, ix);
  assert( ! mvec.empty());
  return min_element(mvec.begin(), 
		     mvec.end(), 
		     [](const GPos& a, const GPos& b) {return a.mult < b.mult;} )->mult;
}

// =============================================================================
  int Mesh2D::knotIntervalFuzzy(Direction2D d, double& par, double eps) const
// =============================================================================
{
  int ix = Mesh2DUtils::last_nonlarger_knotvalue_ix(*this, d, par);
  const double *st = (d == XFIXED) ? &knotvals_x_[0] : &knotvals_y_[0];
  int nmb = (d == XFIXED) ? (int)knotvals_x_.size() : (int)knotvals_y_.size();
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
  int Mesh2D::getKnotIdx(Direction2D d, double& par, double eps) const
// =============================================================================
{
  int ix = Mesh2DUtils::last_nonlarger_knotvalue_ix(*this, d, par);
  const double *st = (d == XFIXED) ? &knotvals_x_[0] : &knotvals_y_[0];
  int nmb = (d == XFIXED) ? (int)knotvals_x_.size() : (int)knotvals_y_.size();
  if (par - st[ix] < eps)
    return ix;
  else if (ix < nmb-1 && st[ix+1] - par < eps)
    return ix+1;
  else
    return -1;
}

// =============================================================================
shared_ptr<Mesh2D>  Mesh2D::subMesh(int ix1, int ix2, int iy1, int iy2) const
// =============================================================================
  {
    shared_ptr<Mesh2D> submesh(new Mesh2D());

    // Copy knot vectors
    submesh->knotvals_x_.insert(submesh->knotvals_x_.end(), 
				knotvals_x_.begin()+ix1, knotvals_x_.begin()+ix2+1);
    submesh->knotvals_y_.insert(submesh->knotvals_y_.end(), 
			       knotvals_y_.begin()+iy1, knotvals_y_.begin()+iy2+1);

    // Transfer knot multiplicity information. Note that the domain in both
    // parameter direction is involved for both parameter directions
    // 1. parameter direction
    int ki;
    for (ki=ix1; ki<=ix2; ++ki)
      {
	vector<GPos> ydir;
	vector<GPos> mr = select_meshvec_(XFIXED, ki);
	int mult;
	size_t kj;

	// Find first multiplicity instance. The instance may be connected to
	// a previous knot in the 2. parameter direction
	for (kj=0; kj<mr.size(); ++kj)
	  {
	    if (mr[kj].ix <= iy1)
	      mult = mr[kj].mult;
	    else if (mr[kj].ix > iy1)
	      break;
	  }
	// Note that the knot values in the new mesh starts from index zero
	ydir.push_back(GPos(0, mult));

	// Copy sucsessive multiplicity instances and modify knot index
	for(; kj<mr.size(); ++kj)
	  {
	    if (mr[kj].ix < iy2)
	      ydir.push_back(GPos(mr[kj].ix-iy1, mr[kj].mult));
	    else
	      break;
	  }
	submesh->mrects_x_.push_back(ydir);
      }

    // 2. parameter direction
    for (ki=iy1; ki<=iy2; ++ki)
      {
	vector<GPos> xdir;
	vector<GPos> mr = select_meshvec_(YFIXED, ki);
	int mult;
	size_t kj;
	  for (kj=0; kj<mr.size(); ++kj)
	  {
	    if (mr[kj].ix <= ix1)
	      mult = mr[kj].mult;
	    else if (mr[kj].ix > ix1)
	      break;
	  }
	xdir.push_back(GPos(0, mult));
	for(; kj<mr.size(); ++kj)
	  {
	    if (mr[kj].ix < ix2)
	      xdir.push_back(GPos(mr[kj].ix-ix1, mr[kj].mult));
	    else
	      break;
	  }
	submesh->mrects_y_.push_back(xdir);
      }

    return submesh;
  }

// =============================================================================
namespace { // anonymous namespace
// =============================================================================

// =============================================================================
bool mrvec_is_correct(const vector<GPos>& v)
// =============================================================================
{
  int prev_ix = -1; // this would always be smaller than first index, which is zero
  for (auto m = v.begin(); m != v.end(); ++m) {
    if      (m->ix <= prev_ix) return false; // should be _strictly_ increasing.
    else if (m->mult < 0)      return false; // multiplicities should never be negative
    prev_ix = m->ix;
  }
  return true;
}



}; // end anonymous namespace

}; // end namespace Go

