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

#include "GoTools/lrsplines3D/LRSplineVolume.h"
#include "GoTools/lrsplines3D/LinDepUtils3D.h"
#include "GoTools/lrsplines3D/Direction3D.h"
#include <algorithm>
#include <iostream>
#include <fstream>

#define DEBUG

using namespace std;
using namespace Go;
namespace
{
  // Private functions
  bool overloadedMeshRectangles(const LRSplineVolume& vol);
  
  int overloadedKnotTuples(const LRSplineVolume& vol);
  
  bool initializeOverload(const LRSplineVolume& vol);

}; // end anonumous namespace


  // Defines a structure that represents a mesh rectangle in three
  // parametric dimensions. NOTE: the mesh rectangle refers to the
  // tensor product expansion of the LR mesh. The tensor product
  // expanded mesh have numDistinctKnots(XDIR) x
  // numDistinctKnots(YDIR) x numDistinctKnots(ZDIR) "distinct" mesh
  // rectangles, and each one
  // of these in a given Direction2D is repeated according to its
  // multiplicity, i.e., between 0 and degree(Direction2D)+1 times. The
  // operator< is needed for sorting when used in an STL map. Note
  // that the ordering of mesh rectangles resulting from the sort
  // implementation effectively corresponds to four nested loops:
  // direction (outer-most), v-knots, u-knots, and multiplicity
  // (inner-most).
struct MeshRectangle {
    Direction3D dir; // direction of mesh rectangle
    int umin;               // index of lower parameter knot in u
    int vmin;               // index of lower parameter knot in v
    int wmin;               // index of lower parameter knot in w
    int mult;               // multiplicity of mesh rectangle
    inline bool operator<( const MeshRectangle rhs) const {
      return 
        (dir  < rhs.dir)  ? true  : // compare direction
        (dir  > rhs.dir)  ? false :
        (umin < rhs.umin) ? true  : // compare u-knot value
        (umin > rhs.umin) ? false :
        (vmin < rhs.vmin) ? true  : // compare v-knot value
        (vmin > rhs.vmin) ? false :
        (wmin < rhs.wmin) ? true  : // compare w-knot value
        (wmin > rhs.wmin) ? false :
        (mult < rhs.mult) ? true  : // compare multiplicity
        (mult > rhs.mult) ? false :
                            false;  // all the same!
    }
  };

struct KnotTuple
{
  int u_ix_;
  int v_ix_;
  int w_ix_;
  vector<LRBSpline3D*> bsplines_;

  KnotTuple(int u_ix, int v_ix, int w_ix, LRBSpline3D *bspl)
  {
    u_ix_ = u_ix;
    v_ix_ = v_ix;
    w_ix_ = w_ix;
    bsplines_.push_back(bspl);
  }

  bool hasKnotTuple(int u_ix, int v_ix, int w_ix)
  {
    if (u_ix == u_ix_ && v_ix == v_ix_ && w_ix == w_ix_)
      return true;
    else
      return false;
  }

  void addBspline(LRBSpline3D *bspl)
  {
    bsplines_.push_back(bspl);
  }
};

typedef unsigned int Index;
typedef map<MeshRectangle, Index> MeshRectangleIndexMap;

//============================================================================
// Fetch all unpeelable B-splines
// That is: B-splines where all corresponding elements are overloaded.
// An element is overloaded if it lies in the support of at least n B-splines
// where n = (degree + 1 in 1. parameter direction) x (degree + 1 in second
// parameter direction) x (degree + 1 in third arameter direction).
// In addition must all of these B-splines be overloaded.
// That is: All elements in the support of these B-splines lies in the support
// of more than n B-splines
// Finally the number of candidate B-splines is at least minnmb and these
//  B-splines have overloaded mesh line segments
vector<LRBSpline3D*> LinDepUtils3D::fetchUnpeelable( const LRSplineVolume& vol,
						   int minnmb)
//============================================================================
{
  vector<LRBSpline3D*> fun;
  
  // Initialize elements
  bool overload = initializeOverload(vol);
  if (!overload)
    return fun;

#ifdef DEBUG
  std::ofstream of1("overload_cand2_0.g2");
  for (auto it2=vol.basisFunctionsBegin(); it2!=vol.basisFunctionsEnd(); ++it2)
    {
      bool curr = it2->second->getOverload();
      if (curr)
  	{
  	  //std::cout << it2->second.get() << std::endl;
  	  LRBSpline3D *cand = it2->second.get();
  	  of1 << "410 1 0 4 255 0 0 255" << std::endl;
  	  of1 << "12" << std::endl;
  	  of1 << cand->umin() << " " << cand->vmin() << cand->wmin();
  	  of1 << cand->umax() << " " << cand->vmin() << cand->wmin() << std::endl;
  	  of1 << cand->umin() << " " << cand->vmin() << cand->wmin();
  	  of1 << cand->umin() << " " << cand->vmax() << cand->wmin() << std::endl;
  	  of1 << cand->umax() << " " << cand->vmin() << cand->wmin();
  	  of1 << cand->umax() << " " << cand->vmax() << cand->wmin() << std::endl;
  	  of1 << cand->umin() << " " << cand->vmax() << cand->wmin();
  	  of1 << cand->umax() << " " << cand->vmax() << cand->wmin() << std::endl;
  	  of1 << cand->umin() << " " << cand->vmin() << cand->wmax();
  	  of1 << cand->umax() << " " << cand->vmin() << cand->wmax() << std::endl;
  	  of1 << cand->umin() << " " << cand->vmin() << cand->wmax();
  	  of1 << cand->umin() << " " << cand->vmax() << cand->wmax() << std::endl;
  	  of1 << cand->umax() << " " << cand->vmin() << cand->wmax();
  	  of1 << cand->umax() << " " << cand->vmax() << cand->wmax() << std::endl;
  	  of1 << cand->umin() << " " << cand->vmax() << cand->wmax();
  	  of1 << cand->umax() << " " << cand->vmax() << cand->wmax() << std::endl;
  	  of1 << cand->umin() << " " << cand->vmin() << cand->wmin();
  	  of1 << cand->umin() << " " << cand->vmin() << cand->wmax() << std::endl;
   	  of1 << cand->umax() << " " << cand->vmin() << cand->wmin();
  	  of1 << cand->umax() << " " << cand->vmin() << cand->wmax() << std::endl;
  	  of1 << cand->umin() << " " << cand->vmax() << cand->wmin();
  	  of1 << cand->umin() << " " << cand->vmax() << cand->wmax() << std::endl;
  	  of1 << cand->umax() << " " << cand->vmax() << cand->wmin();
  	  of1 << cand->umax() << " " << cand->vmax() << cand->wmax() << std::endl;
    	}
     }
  std::cout << std::endl;
#endif

  bool changed = true;
  while (changed)
    {
      changed = false;

      overload = overloadedKnotTuples(vol);
      
      // Reset element flag
      for (auto it1=vol.elementsBegin(); it1!=vol.elementsEnd(); ++it1)
	{
	  bool curr = it1->second->getOverload();
	  if (curr)
	    {
	      // Element2D::resetOverload sets the overload flag in the current
	      // element and removes the overload flag in the supporting B-splines
	      // if the element is not overloaded
	      bool found = it1->second->resetOverload();
	      if (found)
		overload = true;
	      if (curr != found)
		changed = true;
	    }
	}
      if (!overload)
	break;
    }

#ifdef DEBUG
  std::ofstream of2("overload_cand2_1.g2");
  for (auto it2=vol.basisFunctionsBegin(); it2!=vol.basisFunctionsEnd(); ++it2)
    {
      bool curr = it2->second->getOverload();
      if (curr)
  	{
  	  std::cout << it2->second.get() << std::endl;
  	  LRBSpline3D *cand = it2->second.get();
  	  of2 << "410 1 0 4 255 0 0 255" << std::endl;
  	  of2 << "12" << std::endl;
  	  of2 << cand->umin() << " " << cand->vmin() << cand->wmin();
  	  of2 << cand->umax() << " " << cand->vmin() << cand->wmin() << std::endl;
  	  of2 << cand->umin() << " " << cand->vmin() << cand->wmin();
  	  of2 << cand->umin() << " " << cand->vmax() << cand->wmin() << std::endl;
  	  of2 << cand->umax() << " " << cand->vmin() << cand->wmin();
  	  of2 << cand->umax() << " " << cand->vmax() << cand->wmin() << std::endl;
  	  of2 << cand->umin() << " " << cand->vmax() << cand->wmin();
  	  of2 << cand->umax() << " " << cand->vmax() << cand->wmin() << std::endl;
  	  of2 << cand->umin() << " " << cand->vmin() << cand->wmax();
  	  of2 << cand->umax() << " " << cand->vmin() << cand->wmax() << std::endl;
  	  of2 << cand->umin() << " " << cand->vmin() << cand->wmax();
  	  of2 << cand->umin() << " " << cand->vmax() << cand->wmax() << std::endl;
  	  of2 << cand->umax() << " " << cand->vmin() << cand->wmax();
  	  of2 << cand->umax() << " " << cand->vmax() << cand->wmax() << std::endl;
  	  of2 << cand->umin() << " " << cand->vmax() << cand->wmax();
  	  of2 << cand->umax() << " " << cand->vmax() << cand->wmax() << std::endl;
  	  of2 << cand->umin() << " " << cand->vmin() << cand->wmin();
  	  of2 << cand->umin() << " " << cand->vmin() << cand->wmax() << std::endl;
   	  of2 << cand->umax() << " " << cand->vmin() << cand->wmin();
  	  of2 << cand->umax() << " " << cand->vmin() << cand->wmax() << std::endl;
  	  of2 << cand->umin() << " " << cand->vmax() << cand->wmin();
  	  of2 << cand->umin() << " " << cand->vmax() << cand->wmax() << std::endl;
  	  of2 << cand->umax() << " " << cand->vmax() << cand->wmin();
  	  of2 << cand->umax() << " " << cand->vmax() << cand->wmax() << std::endl;
   	}
     }
  std::cout << std::endl;
#endif
  
  if (!overload)
    return fun;

  int nmb = 0;
  for (auto it2=vol.basisFunctionsBegin(); it2!=vol.basisFunctionsEnd(); ++it2)
    {
      bool curr = it2->second->getOverload();
      if (curr)
	++nmb;
    }
#ifdef DEBUG
  std::cout << "Nmb overloaded pre meshrec: " << nmb << std::endl;
#endif

  changed = true;
  while (changed)
    {
      changed = false;

      overload = (nmb >= minnmb) ? overloadedMeshRectangles(vol) : false;
      if (!overload)
	break;
      
      overload = overloadedKnotTuples(vol);
      
      // Reset element flag
      for (auto it1=vol.elementsBegin(); it1!=vol.elementsEnd(); ++it1)
	{
	  bool curr = it1->second->getOverload();
	  if (curr)
	    {
	      // Element2D::resetOverload sets the overload flag in the current
	      // element and removes the overload flag in the supporting B-splines
	      // if the element is not overloaded
	      bool found = it1->second->resetOverload();
	      if (found)
		overload = true;
	      if (curr != found)
		changed = true;
	    }
	}

      if (!overload)
	break;
    }
  
  if (overload)
    {
      // Collect overloaded Bsplines
      for (auto it2=vol.basisFunctionsBegin(); it2!=vol.basisFunctionsEnd(); ++it2)
	{
	  bool curr = it2->second->getOverload();
	  if (curr)
	    fun.push_back(it2->second.get());
	}
    }

  // Reset flags
  for (auto it1=vol.elementsBegin(); it1!=vol.elementsEnd(); ++it1)
    it1->second->eraseOverload();
  for (auto it2=vol.basisFunctionsBegin(); it2!=vol.basisFunctionsEnd(); ++it2)
    it2->second->eraseOverload();

  return fun;
}
  
//==============================================================================
// Given a set of non-peelable B-splines, check if they can be combined in
// linear dependence relations
//
void LinDepUtils3D::checkOverloaded(int minNmb, vector<LRBSpline3D*>& funs,
				    vector<vector<LRBSpline3D*> >& lindep)
//==============================================================================
{
  // Check input
  size_t nmb_funs = funs.size();
  
  if (funs.size() < minNmb)
    return;

#ifdef DEBUG
  std::ofstream of("overloaded.g2");
  for (size_t ki=0; ki<funs.size(); ++ki)
    {
      of << "410 1 0 4 255 0 0 255" << std::endl;
      of << "4" << std::endl;
      of << funs[ki]->umin() << " " << funs[ki]->vmin() << funs[ki]->wmin();
      of << funs[ki]->umax() << " " << funs[ki]->vmin() << funs[ki]->wmin() << std::endl;
      of << funs[ki]->umin() << " " << funs[ki]->vmin() << funs[ki]->wmin();
      of << funs[ki]->umin() << " " << funs[ki]->vmax() << funs[ki]->wmin() << std::endl;
      of << funs[ki]->umax() << " " << funs[ki]->vmin() << funs[ki]->wmin();
      of << funs[ki]->umax() << " " << funs[ki]->vmax() << funs[ki]->wmin() << std::endl;
      of << funs[ki]->umin() << " " << funs[ki]->vmax() << funs[ki]->wmin();
      of << funs[ki]->umax() << " " << funs[ki]->vmax() << funs[ki]->wmin() << std::endl;
      of << funs[ki]->umin() << " " << funs[ki]->vmin() << funs[ki]->wmax();
      of << funs[ki]->umax() << " " << funs[ki]->vmin() << funs[ki]->wmax() << std::endl;
      of << funs[ki]->umin() << " " << funs[ki]->vmin() << funs[ki]->wmax();
      of << funs[ki]->umin() << " " << funs[ki]->vmax() << funs[ki]->wmax() << std::endl;
      of << funs[ki]->umax() << " " << funs[ki]->vmin() << funs[ki]->wmax();
      of << funs[ki]->umax() << " " << funs[ki]->vmax() << funs[ki]->wmax() << std::endl;
      of << funs[ki]->umin() << " " << funs[ki]->vmax() << funs[ki]->wmax();
      of << funs[ki]->umax() << " " << funs[ki]->vmax() << funs[ki]->wmax() << std::endl;
      of << funs[ki]->umin() << " " << funs[ki]->vmin() << funs[ki]->wmin();
      of << funs[ki]->umin() << " " << funs[ki]->vmin() << funs[ki]->wmax() << std::endl;
      of << funs[ki]->umax() << " " << funs[ki]->vmin() << funs[ki]->wmin();
      of << funs[ki]->umax() << " " << funs[ki]->vmin() << funs[ki]->wmax() << std::endl;
      of << funs[ki]->umin() << " " << funs[ki]->vmax() << funs[ki]->wmin();
      of << funs[ki]->umin() << " " << funs[ki]->vmax() << funs[ki]->wmax() << std::endl;
      of << funs[ki]->umax() << " " << funs[ki]->vmax() << funs[ki]->wmin();
      of << funs[ki]->umax() << " " << funs[ki]->vmax() << funs[ki]->wmax() << std::endl;
     }
#endif
  // To ensure correct nesting level
  // and collect zero depth B-splines
  for (size_t ki=0; ki<funs.size(); ++ki)
    {
      funs[ki]->setNestLevel(-1);
      funs[ki]->computeNestLevel();
      int nestdepth = funs[ki]->getNestLevel();
      if (nestdepth == 0)
	{
	  vector<LRBSpline3D*> zerodepth;
	  zerodepth.push_back(funs[ki]);
	  lindep.push_back(zerodepth);
	}
    }

  // Distribute higher depth B-splines
   for (size_t ki=0; ki<funs.size(); ++ki)
    {
      int nestdepth = funs[ki]->getNestLevel();
      if (nestdepth == 0)
	continue;

      size_t kj;
      for (kj=0; kj<lindep.size(); ++kj)
	{
	  if (lindep[kj][0]->covers(funs[ki]))
	    {
	      lindep[kj].push_back(funs[ki]);
	      break;
	    }
	}
      if (kj == lindep.size())
	std::cout << "Stop here" << std::endl;
    }

   // Check for overlap between linear dependence groups
   for (size_t ki=0; ki<lindep.size(); ++ki)
     {
       for (size_t kj=0; kj<lindep[ki].size(); ++kj)
	 {
	   if (lindep[ki][kj]->getNestLevel() > 0)
	     break;
	   for (size_t kr=ki+1; kr<lindep.size(); )
	     {
	       size_t kh;
	       int depth = 0;
	       for (kh=0; kh<lindep[kr].size(); ++kh)
		 {
		   depth = lindep[kr][kh]->getNestLevel();
		   if (depth > 0)
		     break;
		   if (lindep[ki][kj]->overlaps(lindep[kr][kh]))
		     {
		       lindep[ki].insert(lindep[ki].begin()+kj, lindep[kr].begin(),
					 lindep[kr].begin()+kh);
		       lindep[ki].insert(lindep[ki].end(), lindep[kr].begin()+kh+1,
					 lindep[kr].end());
		       lindep.erase(lindep.begin()+kr);
		       break;
		     }
		 }
	       if (kh==lindep[kr].size() || depth > 0)
		 ++kr;
	     }
	 }
     }

   for (int ka=(int)lindep.size()-1; ka>=0; --ka)
     if (lindep[ka].size() < minNmb)
       lindep.erase(lindep.begin() + ka);
   
#ifdef DEBUG
  std::cout << "Number of linear dependency sources: " << lindep.size() << std::endl;
  for (size_t kj=0; kj<lindep.size(); ++kj)
    {
      for (size_t kr=0; kr<lindep[kj].size(); ++kr)
	{
	  if (lindep[kj][kr]->getNestLevel() > 0)
	    break;
	  std::cout << lindep[kj][kr]->umin() << " " << lindep[kj][kr]->umax() << " ";
	  std::cout << lindep[kj][kr]->vmin() << " " << lindep[kj][kr]->vmax() << " ";
	  std::cout << lindep[kj][kr]->wmin() << " " << lindep[kj][kr]->wmax() << std::endl;
	}
      std::cout << lindep[kj].size() << std::endl;
    }
#endif
}



namespace {


//==============================================================================
// Check if identified overloaded B-splines have overloaded knot tuples
//

int overloadedKnotTuples(const LRSplineVolume& vol)
//==============================================================================
{
  vector<LRBSpline3D*> overload;
  vector<KnotTuple> knot_tuples;
  int numbspl = 0;
  for (auto it2=vol.basisFunctionsBegin(); it2!=vol.basisFunctionsEnd(); ++it2)
    {
      bool curr = it2->second->getOverload();
      if (!curr)
	continue;

      numbspl++;
      overload.push_back(it2->second.get());

      vector<int> kvec1 = it2->second->kvec(XDIR);
      vector<int> kvec2 = it2->second->kvec(YDIR);
      vector<int> kvec3 = it2->second->kvec(XDIR);
     for (size_t ki=0; ki<kvec1.size(); ++ki)
	for (size_t kj=0; kj<kvec2.size(); ++kj)
	  for (size_t kh=0; kh<kvec3.size(); ++kh)
	    {
	      size_t kr;
	      for (kr=0; kr<knot_tuples.size(); ++kr)
		if (knot_tuples[kr].hasKnotTuple(kvec1[ki], kvec2[kj], kvec3[kh]))
		  break;
	      if (kr == knot_tuples.size())
		knot_tuples.push_back(KnotTuple(kvec1[ki], kvec2[kj], kvec3[kh],
						it2->second.get()));
	      else
		knot_tuples[kr].addBspline(it2->second.get());
	    }
    }

  for (size_t ki=0; ki<knot_tuples.size(); ++ki)
    {
      if (knot_tuples[ki].bsplines_.size() == 1 &&
	  knot_tuples[ki].bsplines_[0]->getOverload())
	{
	  knot_tuples[ki].bsplines_[0]->eraseOverload();
	  numbspl--;
	}
    }

  return (numbspl > 0);
}

//==============================================================================
// First round: identify B-splines with only overloaded elements
// 
bool initializeOverload(const LRSplineVolume& vol)
//==============================================================================
{
  // Initialize elements
  bool overload = false;
  int nmb_el_init = 0;
  int expected_nmb = (vol.degree(XDIR)+1)*(vol.degree(YDIR)+1)*(vol.degree(ZDIR)+1);
  for (auto it1=vol.elementsBegin(); it1!=vol.elementsEnd(); ++it1)
    {
      bool found = it1->second->initOverload(expected_nmb);
      if (found)
	{
	  overload = true;
	  nmb_el_init++;
	}
    }

  if (!overload)
    return false;

  // Initialize Bsplines
  overload = false;
  int nmb_bspl_init = 0;
  for (auto it2=vol.basisFunctionsBegin(); it2!=vol.basisFunctionsEnd(); ++it2)
    {
      bool found = it2->second->checkOverload();
      if (found)
	{
	  overload = true;
	  nmb_bspl_init++;
	}
    }

  return overload;
}

//==============================================================================
// Check if identified overloaded B-splines have overloaded mesh rectangles.
// 
bool overloadedMeshRectangles(const LRSplineVolume& vol)
//==============================================================================
{
#ifdef DEBUG
  std::cout << "Overloaded mesh rectangles, start." << std::endl;
#endif
  Direction3D ds[3] = {XDIR, YDIR, ZDIR};
  const vector<Direction3D> DirSeq(ds,ds+3);
  MeshRectangleIndexMap meshrectangles;
  vector<vector<size_t> > bspl;
  int numbspl = 0;
  vector<LRBSpline3D*> overload;
  for (auto it2=vol.basisFunctionsBegin(); it2!=vol.basisFunctionsEnd(); ++it2)
    {
      bool curr = it2->second->getOverload();
      if (!curr)
	continue;

      numbspl++;
      overload.push_back(it2->second.get());
      
      // Collect associated mesh rectangles
      map<Direction3D, vector<int> > kvec;
      map<Direction3D, vector<int> > kmul;
      map<Direction3D, vector<int> > kall;
      for (vector<Direction3D>::const_iterator ixy=DirSeq.begin();
	   ixy!=DirSeq.end(); ++ixy) {
	vector<int> kvec_all = it2->second->kvec(*ixy);
	kvec[*ixy] = kvec_all;
	vector<int>::const_iterator it = unique( kvec[*ixy].begin(), kvec[*ixy].end() );
	kvec[*ixy].resize( it - kvec[*ixy].begin() );
	kmul[*ixy].resize( kvec[*ixy].size() );
	for ( vector<int>::iterator it_mu=kmul[*ixy].begin(), it_uk=kvec[*ixy].begin();
              (it_mu!=kmul[*ixy].end() && it_uk!=kvec[*ixy].end());
	      ++it_mu, ++it_uk) {
	  *it_mu = (int) count(kvec_all.begin(), kvec_all.end(), *it_uk) - 1; // Note: we ensure zero-based multiplicity.
	}
	int kmin = kvec[*ixy].front();
	int kmax = kvec[*ixy].back();
	kall[*ixy].resize( kmax - kmin ); // Note: we exclude the last knot index.
	int kv = kmin;
	for (vector<int>::iterator it_al=kall[*ixy].begin();
	     it_al!=kall[*ixy].end(); ++it_al)
	  *it_al = kv++;
      }
      
      // Using the above information for this B-spline, we proceed
      // to update the MeshPart/column-wise storage of the 
      // incidence matrix.
      for (vector<Direction3D>::const_iterator 
	     ixy=DirSeq.begin(); ixy!=DirSeq.end(); ++ixy)
	{ // Loop: over directions
	  // When doing mesh rectangles in a given Direction2D, we must
	  // loop first over the (p+2) knot indices of the B-spline in
	  // that Direction2D, and then for each of these loop over ALL
	  // but the last knot indices spanned by the B-spline.
	  vector<int> KV = (*ixy==YDIR) ? kvec[YDIR] : kall[YDIR];
	  vector<int> KU = (*ixy==XDIR) ? kvec[XDIR] : kall[XDIR];
	  vector<int> KW = (*ixy==ZDIR) ? kvec[ZDIR] : kall[ZDIR];
	  for (unsigned int ikw=0; ikw!=KW.size(); ++ikw)
	    { // Loop: over w-knots
	      for (unsigned int ikv=0; ikv!=KV.size(); ++ikv)
		{ // Loop: over v-knots
		  for (unsigned int iku=0; iku!=KU.size(); ++iku)
		    {// Loop over u-knots
		      int nmb = kmul[*ixy][(*ixy==XDIR) ? iku :
					   ((*ixy==YDIR) ? ikv : ikw)];
		      for (int nu=0; nu<=nmb; ++nu)
			{ // Loop: over multiplicity
			  MeshRectangle tmprec;
			  tmprec.dir = *ixy;
			  tmprec.wmin = KW[ikw];
			  tmprec.vmin = KV[ikv];
			  tmprec.umin = KU[iku];
			  tmprec.mult = nu;

			  // Check if the mesh rectangle exists already
			  auto it = meshrectangles.find(tmprec);
			  if (it == meshrectangles.end())
			    {
			      // Insert new meshrectangle
			      meshrectangles[tmprec] = bspl.size();
			      vector<size_t> tmp;
			      tmp.push_back(overload.size()-1);
			      bspl.push_back(tmp);
			    }
			  else
			    {
			      bspl[it->second].push_back(overload.size()-1);
			    }
			}
		    }
		}
	    }
	}
    }

#ifdef DEBUG
  std::cout << "Overloaded mesh rectangles, middle. Found: " << numbspl << std::endl;
#endif
  // Remove mesh rectangles that contain less than two overloaded B-splines
  // and reset the overload flag for associated B-splines
  bool changed = true;
  vector<bool> on(bspl.size(), true);
  vector<size_t> num(bspl.size());
  for (size_t ki=0; ki<bspl.size(); ++ki)
    num[ki] = bspl[ki].size();

  size_t num2 = bspl.size();
  while (changed)
    {
#ifdef DEBUG
      std::cout << num2 << ", ";
#endif
      changed = false;
      for (size_t ki=0; ki<bspl.size(); ++ki)
	{
	  if (!on[ki])
	    continue;
	  if (num[ki] < 2)
	    {
	      for (size_t kj=0; kj<num[ki]; ++kj)
		{
		  for (size_t kr=0; kr<bspl.size(); ++kr)
		    {
		      if (!on[ki])
			continue;
		      if (kr == ki)
			continue;
		      //bspl[kr].remove(bspl[ki][kj]);
		      auto it = std::find(bspl[kr].begin(), bspl[kr].begin()+num[kr], bspl[ki][kj]);
		      if (it != bspl[kr].begin()+num[kr])
			{
			  std::swap(bspl[kr][num[kr]-1], *it);
			  --num[kr];
			}
			//bspl[kr].erase(it);
		    }
		  overload[bspl[ki][kj]]->eraseOverload();
		}
	      changed = true;
	      on[ki] = false;
	      -num2;
	    }
	}
    }
#ifdef DEBUG
  std::cout << std::endl << "Overloaded mesh rectangles, finish" << std::endl;
#endif
  return (bspl.size() > 0);
}



}; // end anonymous namespace
