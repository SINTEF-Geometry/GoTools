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

#include "GoTools/lrsplines2D/LRSplineSurface.h"
#include "GoTools/lrsplines2D/LinDepUtils.h"
#include "GoTools/lrsplines2D/Direction2D.h"
#include <algorithm>
#include <iostream>

using namespace std;
using namespace Go;


//#define DEBUG

  // Defines a structure that represents a mesh rectangle in two
  // parametric dimensions. NOTE: the mesh rectangle refers to the
  // tensor product expansion of the LR mesh. The tensor product
  // expanded mesh have numDistinctKnots(XFIXED) x
  // numDistinctKnots(YFIXED) "distinct" mesh rectangles, and each one
  // of these in a given Direction2D is repeated according to its
  // multiplicity, i.e., between 0 and degree(Direction2D)+1 times. The
  // operator< is needed for sorting when used in an STL map. Note
  // that the ordering of mesh rectangles resulting from the sort
  // implementation effectively corresponds to four nested loops:
  // direction (outer-most), v-knots, u-knots, and multiplicity
  // (inner-most).
struct MeshRectangle {
    Direction2D dir; // direction of mesh rectangle
    int vmin;               // index of lower parameter knot in v
    int umin;               // index of lower parameter knot in u
    int mult;               // multiplicity of mesh rectangle
    inline bool operator<( const MeshRectangle rhs) const {
      return 
        (dir  < rhs.dir)  ? true  : // compare direction
        (dir  > rhs.dir)  ? false :
        (vmin < rhs.vmin) ? true  : // compare v-knot value
        (vmin > rhs.vmin) ? false :
        (umin < rhs.umin) ? true  : // compare u-knot value
        (umin > rhs.umin) ? false :
        (mult < rhs.mult) ? true  : // compare multiplicity
        (mult > rhs.mult) ? false :
                            false;  // all the same!
    }
  };


  //============================================================================
  // Type definitions to shorten the subsequent notation. 

  //   Index: Used to loop over Elements (represented by pointers
  //     to elements) and B-splines (represented by pointers to
  //     LRBSpline2Ds) through the IncidenceMatrix, and through
  //     SwitchVectors.

  //   ElementIndexMap: Maps each Element (represented by a
  //     pointer to an Element2D) to a linear Index.

  //   BsplineIndexMap: Maps each B-spline (represented by a pointer
  //     to a LRBSpline2D) to a linear Index.

  //   MeshRectangleIndexMap: Maps each mesh rectangle (represented by
  //     a the MeshRectangle structure) to a linear Index.

  //   SwitchVector: Used to construct two "pseudo-bool" vectors, one
  //     for each (Index of) MeshPart (element or mesh rectangle) and
  //     one for each (Index of) B-spline. A value of ON means that
  //     whatever the Index refers to is "ON" (potentially relevant
  //     for linear dependence), while a value of OFF means that
  //     whatever the Index refers to is "OFF" (irrelevant for linear
  //     dependence).

  //   IncidenceMatrix: holds the information contained in the
  //     incidence matrices, cf. [Dokken, Lyche & Pettersen,
  //     2012]. The first incidence matrix essentially maps B-splines
  //     to elements and vice versa. The data can be stored
  //     Element/Row-wise, i.e., with the i'th entry holding a vector
  //     containing the Indices of the B-splines with support on the
  //     i'the element, or B-spline/Row-wise, i.e., with the i'th
  //     entry holding a vector containing the Indices of the Elements
  //     on which the i'th B-spline has support. The second incidence
  //     matrix essentially maps B-splines to mesh rectangles and vice
  //     versa.

  //  DirSeq: Constant vector with ordering of directions to loop over
  //     when processing MeshRectangles.

  //============================================================================
  typedef unsigned int Index;
  typedef std::map<Element2D*, Index> ElementIndexMap;
  typedef std::map<LRBSpline2D*, Index> BsplineIndexMap;
  typedef map<MeshRectangle, Index> MeshRectangleIndexMap;
  typedef std::vector<short > SwitchVector;
  const SwitchVector::value_type ON = 1, OFF = 0;
  typedef vector<vector<Index> > IncidenceMatrix;
  Direction2D ds[2] = {XFIXED, YFIXED};
  const vector<Direction2D> DirSeq(ds,ds+2);

  //============================================================================
  // Constructs a map from an Element to a simple Index.
  //============================================================================
  ElementIndexMap construct_elementindex_map( const LRSplineSurface & lrs)
  {
    ElementIndexMap result;
    Index i = 0;
    for (LRSplineSurface::ElementMap::const_iterator it_el=lrs.elementsBegin(); 
        it_el!=lrs.elementsEnd(); ++it_el)
      result[it_el->second.get()] = i++;
    return result;
  }

  //============================================================================
  // Constructs a map from a B-spline to a simple Index.
  // ============================================================================
  BsplineIndexMap construct_bsplineindex_map( const LRSplineSurface & lrs )
  {
    BsplineIndexMap result;
    Index i = 0;
    for (LRSplineSurface::BSplineMap::const_iterator it_bs=lrs.basisFunctionsBegin(); 
        it_bs!=lrs.basisFunctionsEnd(); ++it_bs)
      result[it_bs->second.get()] = i++;
    return result;
  }

  //============================================================================
  // Computes the Element/Column-wise storage (m1) and the
  // Bspline/Row-wise storage (b1) of the first incidence matrix
  // (FIM).
  // ============================================================================
  void construct_first_incidence_matrix( const LRSplineSurface& lrs,
      const ElementIndexMap& MImap, 
      const BsplineIndexMap& BSmap,
      IncidenceMatrix& m1, 
      IncidenceMatrix& b1 )
  {
    // Initialize the Element/Column-wise storage (m1) of the FIM, the
    // B-spline/Row-wise storage (b1) of FIM, and fill in data by
    // looping over first elements and then B-splines.
    m1.resize( lrs.numElements() );
    b1.resize( lrs.numBasisFunctions() );
    for (LRSplineSurface::ElementMap::const_iterator it_el=lrs.elementsBegin();
        it_el!=lrs.elementsEnd(); ++it_el) { // Loop 0: over elements.
      Index in_el = MImap.at(it_el->second.get());
      for (vector<LRBSpline2D*>::const_iterator it_bs = 
	     it_el->second->supportBegin();
	   it_bs != it_el->second->supportEnd(); ++it_bs)
	{
	  Index in_bs = BSmap.at((*it_bs));
	  // We store this (Element, B-spline)-pair in the m1
	  // representation and its transpose b1 representation.
	  m1[in_el].push_back(in_bs);
	  b1[in_bs].push_back(in_el);
	}
    }
    return;
  }

  //============================================================================
  // Counts the number of unpeeled/ON entries of a given SwitchVector.
  // ============================================================================
  Index numEntriesOn( const SwitchVector& onoff )
  {
    Index result = 0;
    for (Index i=0; i!=onoff.size(); ++i)
      if (onoff.at(i)==ON) // If this entry is ON ...
	result++;          // ... add one to the counter.
    return result;
  }

  //============================================================================
  // Identifies overloaded elements and B-splines, and peels
  // non-overloaded B-splines. This is done by updating the element
  // (ele_is_on) and B-spline (fun_is_on) SwitchVectors. After return,
  // all overloaded elements iE will have ele_is_on[iE]=ON, and all
  // overloaded B-splines iB wil have fun_is_on[iB]=ON, while all
  // other entries of ele_is_on and fun_is_on will be OFF.
  // ============================================================================
  void peel_nonoverloaded_functions( const Index maxNumFun,
      const IncidenceMatrix& m1, const IncidenceMatrix& b1,
      SwitchVector& ele_is_on, SwitchVector& fun_is_on) 
  {
    // We first initialize the SwitchVector of elements as OFF
    // (ele_is_on), and initialize the SwitchVector of B-splines as
    // OFF (fun_is_on), and then start a loop over Elements.
    ele_is_on.resize( m1.size(), OFF );
    fun_is_on.resize( b1.size(), OFF );
    for (Index in_el=0; in_el!=m1.size(); ++in_el) {
      // If this element is overloaded, we turn the element ON. Also,
      // we turn all the B-splines with support on it ON, noting that
      // in this way we (probably) overestimate the number of
      // B-splines that are ACTUALLY ON. In this way, we avoid
      // checking all functions later on, altough we (probably) have
      // to switch some functions back OFF.
      if ( m1.at(in_el).size()>maxNumFun ) {
        ele_is_on.at(in_el) = ON;
        for (Index it_bs=0; it_bs!=m1.at(in_el).size(); ++it_bs)
          fun_is_on.at(m1.at(in_el).at(it_bs)) = ON;
      }
    }
    // Count the number of overloaded elements, and make return if
    // there are no overloaded elements.
    if (numEntriesOn(ele_is_on)==0)
      return;
    // We update the B-splines SwitchVector, since we have (probably)
    // switched too many ON above. A B-spline is overloaded (here ON)
    // if all elements in its support are overloaded (here ON).
    for (Index in_bs=0; in_bs!=fun_is_on.size(); ++in_bs)
      if (fun_is_on.at(in_bs)==ON)
        for (Index in_el=0; in_el!=b1.at(in_bs).size(); ++in_el)
          if (ele_is_on.at(b1.at(in_bs).at(in_el))==OFF) {
            // Element in_el which is in the support of function in_bs
            // is OFF (here NOT overloaded), so this function in_bs
            // should also be OFF (here NOT overloaded).
            fun_is_on.at(in_bs) = OFF;
            break;
          }
    // Finally, check if all B-splines are OFF (here NOT overloaded),
    // and return this.
    return;
  }

  //============================================================================
  // Peels overloaded B-splines of a given LR spline using the
  // supplied incidence matrix. The incidence matrix may contain
  // information from elements (M1), from mesh rectangles (M2), or
  // from both (M=[M1,M2]), and below we refer to this "thing"
  // commonly as "MeshPart". The peeling is done by updating the
  // element (ele_is_on) and B-spline (fun_is_on)
  // SwitchVectors. MeshParts iP with only one overloaded B-spline is
  // switched OFF (ele_is_on[iP]=OFF), and the corresponding
  // overloaded B-splines iB are also switched OFF
  // (fun_is_on[iB]=OFF).
  // ============================================================================
  void peel_overloaded_functions( const IncidenceMatrix& m1,
      SwitchVector& ele_is_on,
      SwitchVector& fun_is_on ) 
  {
    bool try_new_peel;
    do { // Loop 0: over peeling iterations.
      try_new_peel = false; // Start by assuming this is the last peel.
      for (Index in_el=0; in_el!=m1.size(); ++in_el ) // Loop 1: over all MeshParts.
        if (ele_is_on.at(in_el)==ON) {
          // This MeshPart is ON, so lets count how many B-splines with
          // support on it that are also ON.
          int num_fun_on = 0; // Number of functions that are ON
          Index in_bs_on;     // Index of last function that was ON
          for (Index in_bs=0; in_bs!=m1[in_el].size(); ++in_bs ) // Loop 2A: over B-splines on MeshPart.
            if (fun_is_on[m1[in_el][in_bs]]==ON) { 
              // This B-spline is ON, so store its index and
              // increment the counter.
              in_bs_on = in_bs;
              if ( (num_fun_on++) > 1 )
                // This MeshPart is covered by more than one ON
                // B-splines, so we cannot peel it, and therefore we
                // skip it.
                break;
              }
          // end loop 2A.
          if (num_fun_on==0) {
            // We know that no ON B-spline covers this MeshPart, so we
            // can forget about it.
            ele_is_on[in_el] = OFF;
          } else if (num_fun_on==1) {
            // We know that exactly one ON B-spline covers this
            // MeshPart, so we CAN peel it. We do this my switching the
            // MeshPart and the B-spline OFF. Also, make sure we try
            // another peel.
            ele_is_on[in_el] = OFF;
            fun_is_on[m1[in_el][in_bs_on]] = OFF;
            try_new_peel = true;
          }
        }
      // end loop 1.
    } while (try_new_peel); // end loop 0
    return;
  }

  //============================================================================
  // Constructs a map from a MeshRectangle to a simple Index.
  //============================================================================
  MeshRectangleIndexMap construct_meshrectangleindex_map( const LRSplineSurface& lrs )
  {
    MeshRectangleIndexMap result;
    Index i = 0;
    for (vector<Direction2D>::const_iterator dirit=DirSeq.begin(); dirit!=DirSeq.end(); ++dirit)
      for (int vknot=0; vknot!=lrs.mesh().numDistinctKnots(YFIXED); ++vknot)
        for (int uknot=0; uknot!=lrs.mesh().numDistinctKnots(XFIXED); ++uknot) {
          int fixed = (*dirit==XFIXED) ? uknot : vknot;
		  int start = (*dirit==XFIXED) ? vknot : uknot;
		  for (int mult = 0; mult<=lrs.mesh().nu(*dirit,fixed,start,start+1); ++mult) {
			  MeshRectangle t_meshrec;
			  t_meshrec.dir = *dirit; t_meshrec.vmin = vknot; t_meshrec.umin = uknot; t_meshrec.mult = mult;
			  result[t_meshrec] = i++;
		  }
		}
		return result;
  }

  //============================================================================
  // Constructs relevant parts of the full incidence matrix M=[M1,M2],
  // by gluing the second mesh rectangle incidence matrix (M2) onto
  // the existing first element incidence matrix (M1) as given in
  // input (M), and prolongs the Mesh SwitchVector by gluing the mesh
  // rectangle switches onto the existing element switches, as given
  // in the input (ele_is_on). Only those rows of the full incidence
  // matrix corresponding to B-splines that are ON are constructed,
  // while all other are skipped, i.e., kept OFF.
  // ============================================================================
  void construct_unpeeled_rows_of_full_incidence_matrix ( const LRSplineSurface lrs, 
    const ElementIndexMap MImap, const MeshRectangleIndexMap MRImap, 
    const SwitchVector fun_is_on, IncidenceMatrix& m, 
    SwitchVector& ele_is_on )
  {
    // Resize the incidence matrix to contain both elements (M1) and
    // mesh rectangles (M2) information, and likewise resize MeshPart
    // SwitchVector to contain both elements (first part) and mesh
    // rectangles (last part).
    m.resize( MImap.size() + MRImap.size() );
    ele_is_on.resize( MImap.size() + MRImap.size(), OFF);
    Index in_bs = 0;
    for (LRSplineSurface::BSplineMap::const_iterator it_bs =
	   lrs.basisFunctionsBegin();  
	 it_bs!=lrs.basisFunctionsEnd(); ++it_bs, ++in_bs) { // Loop: over B-splines
      if (fun_is_on[in_bs]==ON) {
        // This B-spline is ON, so for each parametric direction, we
        // extract its local knot vector indices, find the unique
        // knots (kvec), their multiplicity (kmul), and ALL knots
        // indices in the other parametric direction (kall)
        map<Direction2D, vector<int> > kvec;
        map<Direction2D, vector<int> > kmul;
        map<Direction2D, vector<int> > kall;
        for (vector<Direction2D>::const_iterator ixy=DirSeq.begin();
            ixy!=DirSeq.end(); ++ixy) {
          vector<int> kvec_all = it_bs->second->kvec(*ixy);
          kvec[*ixy] = kvec_all;
          vector<int>::const_iterator it = unique( kvec[*ixy].begin(), kvec[*ixy].end() );
          kvec[*ixy].resize( it - kvec[*ixy].begin() );
          kmul[*ixy].resize( kvec[*ixy].size() );
          for ( vector<int>::iterator it_mu=kmul[*ixy].begin(), it_uk=kvec[*ixy].begin();
              (it_mu!=kmul[*ixy].end() && it_uk!=kvec[*ixy].end());
		++it_mu, it_uk++ ) {
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
        // to update the MeshPart/column-wise storage of the full
        // incidence matrix.
        for (vector<Direction2D>::const_iterator 
            ixy=DirSeq.begin(); ixy!=DirSeq.end(); ++ixy) { // Loop: over directions
          // When doing mesh rectangles in a given Direction2D, we must
          // loop first over the (p+2) knot indices of the B-spline in
          // that Direction2D, and then for each of these loop over ALL
          // but the last knot indices spanned by the B-spline.
          vector<int> KV = (*ixy==YFIXED) ? kvec[YFIXED] : kall[YFIXED];
          vector<int> KU = (*ixy==XFIXED) ? kvec[XFIXED] : kall[XFIXED];
		  for (unsigned int ikv=0; ikv!=KV.size(); ++ikv) { // Loop: over v-knots
			for (unsigned int iku=0; iku!=KU.size(); ++iku) {// Loop over u-knots
			  for (int nu=0; nu<=kmul[*ixy][(*ixy==XFIXED) ? iku : ikv]; ++nu) { // Loop: over multiplicity
				MeshRectangle tmprec;
				tmprec.dir = *ixy; tmprec.vmin = KV[ikv]; tmprec.umin = KU[iku]; tmprec.mult = nu;
				Index in_mr = (int)MImap.size() + MRImap.at(tmprec);
				m.at(in_mr).push_back(in_bs); // Add this function to this mesh rectangle.
		  	    ele_is_on.at(in_mr) = ON;     // Switch this mesh rectangle ON.
			  }
		    }
		  }
		}
      }
    }
    return;
  }

  //============================================================================
  // Tests for potential linear dependence among the LR B-splines of a
  // given LR spline. To be more precise: Tests whether a given LR
  // spline is peelable. We say an LR spline is peelable if the
  // incidence matrix contains NO non-zero entries after peeling it,
  // i.e., if ALL the overloaded LR B-splines are peelable,
  // cf. [Dokken, Lyche & Pettersen, 2012]. The LR spline PEELABILITY
  // is a SUFFICIENT condition for linear INDEPENDENCE of the LR
  // B-splines, or equivalently, the LR spline UNPEELABILITY is a
  // NECESSARY condition for linear DEPENDENCE of the LR B-splines. If
  // an LR spline IS peelable, the LR B-splines are linearly
  // INdependent. If the LR spline is NOT peelable, the LR B-splines
  // MAY be linearly dependent, and further investigations are
  // required to determine whether some of the LR B-splines are
  // ACTUALLY part of a linear dependence relation.
  // ============================================================================
  bool LinDepUtils::isPeelable( const LRSplineSurface& lrs ) {
    vector<LRBSpline2D*> funs = unpeelableBasisFunctions(lrs);
    return (funs.size()==0);
  }

  //============================================================================
  // Finds unpeelable overloaded LR B-splines of a given LR spline. An
  // LR B-spline is unpeelable, if it cannot be peeled as described in
  // [Dokken, Lyche & Pettersen, 2012]. These unpeeplable LR B-splines
  // corresponds to rows in the incidence matrix with non-zero
  // entries.
  // ============================================================================
  vector<LRBSpline2D*> LinDepUtils::unpeelableBasisFunctions ( 
    const LRSplineSurface& lrs )
  {
    // Initiate empty vector of unpeelable functions (result),
    // construct maps from elements (MImap) and B-splines (BSmap) into
    // linear indices, and set up the first incidence matrix, stored
    // both element/column-wise (m1) and B-spline/row-wise (b1).
    vector<LRBSpline2D*> result;
    const ElementIndexMap MImap = construct_elementindex_map( lrs );
    const BsplineIndexMap BSmap = construct_bsplineindex_map( lrs );
    IncidenceMatrix m1, b1;
    construct_first_incidence_matrix( lrs, MImap, BSmap, m1, b1);
    // Define the maximum number of B-splines allowed on an element
    // defining overloading, initiate empty element (ele_is_on) and
    // B-spline (fun_is_on) SwitchVectors, peel away all
    // non-overloaded B-splines, and make return if we have peeled
    // all, i.e., if there are NO overloaded B-splines. NOTE: this is
    // a more aggresive initialization of the peeling than what is
    // stated in the algorithm in [Dokken, Lyche & Pettersen]!
    SwitchVector ele_is_on;
    SwitchVector fun_is_on;
    const Index maxNumFun = (lrs.degree(XFIXED)+1) * (lrs.degree(YFIXED)+1);
    peel_nonoverloaded_functions( maxNumFun, m1, b1, ele_is_on, fun_is_on );
    if ( (numEntriesOn(fun_is_on))==0 )
      return result;
    // If we have reached this point, there ARE overloaded B-splines.
    // We therefore try to peel as many of the these as possible based
    // on elements ONLY, and make return, if we have peeled all, i.e.,
    // if there are NO element-wise unpeelable overloaded B-splines.
    peel_overloaded_functions( m1, ele_is_on, fun_is_on );
    if ( (numEntriesOn(fun_is_on))==0 )
      return result;
    // If we have reached this point, there are element-wise
    // unpeelable overloaded B-splines. We therefore try peeling based
    // on BOTH elements AND mesh rectangles. First, construct the map
    // from MeshRectangles to a linear index, then pad (ON rows of)
    // the second incidence matrix onto (ON rows of) the first
    // incidence matrix to get (ON rows of) the full incidence matrix,
    // and likewise pad the element SwitchVector onto the
    // MeshRectangle SwitchVector to get the full MeshPart
    // SwitchVector, and finally try peeling the full system, and make
    // return, if we have peeled all, i.e., if there are NO unpeelable
    // overloaded B-splines.
    const MeshRectangleIndexMap MRImap = construct_meshrectangleindex_map( lrs );    
    construct_unpeeled_rows_of_full_incidence_matrix( lrs, MImap, MRImap, fun_is_on, m1, ele_is_on );
    peel_overloaded_functions( m1, ele_is_on, fun_is_on );
    if ( (numEntriesOn(fun_is_on))==0 )
      return result;
    // If we have reached this point, there ARE unpeelable overloaded
    // B-splines. We therefore return these.
    Index int_bs = 0;
    for (BsplineIndexMap::const_iterator it_bs=BSmap.begin();
	 it_bs!=BSmap.end(); ++it_bs)   // Loop over: B-splines
      if (fun_is_on[int_bs++]==ON)      // If this function is ON ...
        result.push_back(it_bs->first); // ... add it to the list of functions
    return result;
  }

//============================================================================
// Fetch all unpeelable B-splines
// That is: B-splines where all corresponding elements are overloaded.
// An element is overloaded if it lies in the support of at least n B-splines
// where n = (degree + 1 in 1. parameter direction) x (degree + 1 in second
// parameter direction). In addition must all of these B-splines be overloaded.
// That is: All elements in the support of these B-splines lies in the support
// of more than n B-splines
// Finally the number of candidate B-splines is at least minnmb and these
//  B-splines have overloaded mesh line segments
vector<LRBSpline2D*> LinDepUtils::fetchUnpeelable( const LRSplineSurface& surf,
						   int minnmb)
//============================================================================
{
  vector<LRBSpline2D*> fun;
  
  // Initialize elements
  bool overload = false;
  int nmb_el_init = 0;
  int expected_nmb = (surf.degree(XFIXED)+1)*(surf.degree(YFIXED)+1);
  for (auto it1=surf.elementsBegin(); it1!=surf.elementsEnd(); ++it1)
    {
      bool found = it1->second->initOverload(expected_nmb);
      if (found)
	{
	  overload = true;
	  nmb_el_init++;
	}
    }

  if (!overload)
    return fun;

  // Initialize Bsplines
  overload = false;
  int nmb_bspl_init = 0;
  for (auto it2=surf.basisFunctionsBegin(); it2!=surf.basisFunctionsEnd(); ++it2)
    {
      bool found = it2->second->checkOverload();
      if (found)
	{
	  overload = true;
	  nmb_bspl_init++;
	}
    }

  if (!overload)
    return fun;

#ifdef DEBUG
  std::ofstream of1("overload_cand2_0.g2");
  for (auto it2=surf.basisFunctionsBegin(); it2!=surf.basisFunctionsEnd(); ++it2)
    {
      bool curr = it2->second->getOverload();
      if (curr)
  	{
  	  std::cout << it2->second.get() << std::endl;
  	  LRBSpline2D *cand = it2->second.get();
  	  of1 << "410 1 0 4 255 0 0 255" << std::endl;
  	  of1 << "4" << std::endl;
  	  of1 << cand>umin() << " " << cand>vmin() << " 0 ";
  	  of1 << cand>umax() << " " << cand>vmin() << " 0" << std::endl;
  	  of1 << cand>umin() << " " << cand>vmin() << " 0 ";
  	  of1 << cand>umin() << " " << cand>vmax() << " 0" << std::endl;
  	  of1 << cand>umax() << " " << cand>vmin() << " 0 ";
  	  of1 << cand>umax() << " " << cand>vmax() << " 0" << std::endl;
  	  of1 << cand>umin() << " " << cand>vmax() << " 0 ";
  	  of1 << cand>umax() << " " << cand>vmax() << " 0" << std::endl;
  	}
     }
  std::cout << std::endl;
  writeg2Mesh(surf, of1);
#endif

  bool changed = true;
  while (changed)
    {
      changed = false;

      // Reset element flag
      overload = false;
      for (auto it1=surf.elementsBegin(); it1!=surf.elementsEnd(); ++it1)
	{
	  bool curr = it1->second->getOverload();
	  if (curr)
	    {
	      bool found = it1->second->resetOverload();
	      if (found)
		overload = true;
	      if (curr != found)
		changed = true;
	    }
	}

      // if (!overload)
      // 	break;
      // if (!changed)
      // 	break;

      // // Reset Bspline flag
      // overload = false;
      // for (auto it2=surf.basisFunctionsBegin(); it2!=surf.basisFunctionsEnd(); ++it2)
      // 	{
      // 	  bool curr = it2->second->getOverload();
      // 	  if (curr)
      // 	    {
      // 	      bool found = it2->second->checkOverload();
      // 	      if (found)
      // 		overload = true;
      // 	      if (curr != found)
      // 		changed = true;
      // 	    }
      // 	}
    }

#ifdef DEBUG
  std::ofstream of2("overload_cand2_1.g2");
  for (auto it2=surf.basisFunctionsBegin(); it2!=surf.basisFunctionsEnd(); ++it2)
    {
      bool curr = it2->second->getOverload();
      if (curr)
  	{
  	  std::cout << it2->second.get() << std::endl;
  	  LRBSpline2D *cand = it2->second.get();
  	  of2 << "410 1 0 4 255 0 0 255" << std::endl;
  	  of2 << "4" << std::endl;
  	  of2 << cand>umin() << " " << cand>vmin() << " 0 ";
  	  of2 << cand>umax() << " " << cand>vmin() << " 0" << std::endl;
  	  of2 << cand>umin() << " " << cand>vmin() << " 0 ";
  	  of2 << cand>umin() << " " << cand>vmax() << " 0" << std::endl;
  	  of2 << cand>umax() << " " << cand>vmin() << " 0 ";
  	  of2 << cand>umax() << " " << cand>vmax() << " 0" << std::endl;
  	  of2 << cand>umin() << " " << cand>vmax() << " 0 ";
  	  of2 << cand>umax() << " " << cand>vmax() << " 0" << std::endl;
  	}
     }
  std::cout << std::endl;
  writeg2Mesh(surf, of2);
#endif
  
  if (!overload)
    return fun;

  int nmb = 0;
  for (auto it2=surf.basisFunctionsBegin(); it2!=surf.basisFunctionsEnd(); ++it2)
    {
      bool curr = it2->second->getOverload();
      if (curr)
	++nmb;
    }
#ifdef DEBUG
  std::cout << "Nmb overloaded pre meshrec: " << nmb << std::endl;
#endif
  
  overload = (nmb >= minnmb) ? overloadedMeshRectangles(surf) : false;
  
  if (overload)
    {
      // Collect overloaded Bsplines
      for (auto it2=surf.basisFunctionsBegin(); it2!=surf.basisFunctionsEnd(); ++it2)
	{
	  bool curr = it2->second->getOverload();
	  if (curr)
	    fun.push_back(it2->second.get());
	}
    }

  // Reset flags
  for (auto it1=surf.elementsBegin(); it1!=surf.elementsEnd(); ++it1)
    it1->second->eraseOverload();
  for (auto it2=surf.basisFunctionsBegin(); it2!=surf.basisFunctionsEnd(); ++it2)
    it2->second->eraseOverload();

  return fun;
}


//==============================================================================
// Check if identified overloaded B-splines have overloaded mesh rectangles.
// 
bool LinDepUtils::overloadedMeshRectangles(const LRSplineSurface& surf)
//==============================================================================
{
#ifdef DEBUG
  std::cout << "Overloaded mesh rectangles, start." << std::endl;
#endif
  Direction2D ds[2] = {XFIXED, YFIXED};
  const vector<Direction2D> DirSeq(ds,ds+2);
  MeshRectangleIndexMap meshrectangles;
  vector<vector<size_t> > bspl;
  int numbspl = 0;
  vector<LRBSpline2D*> overload;
  for (auto it2=surf.basisFunctionsBegin(); it2!=surf.basisFunctionsEnd(); ++it2)
    {
      bool curr = it2->second->getOverload();
      if (!curr)
	continue;

      numbspl++;
      overload.push_back(it2->second.get());
      
      // Collect associated mesh rectangles
      map<Direction2D, vector<int> > kvec;
      map<Direction2D, vector<int> > kmul;
      map<Direction2D, vector<int> > kall;
      for (vector<Direction2D>::const_iterator ixy=DirSeq.begin();
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
      for (vector<Direction2D>::const_iterator 
	     ixy=DirSeq.begin(); ixy!=DirSeq.end(); ++ixy)
	{ // Loop: over directions
	  // When doing mesh rectangles in a given Direction2D, we must
	  // loop first over the (p+2) knot indices of the B-spline in
	  // that Direction2D, and then for each of these loop over ALL
	  // but the last knot indices spanned by the B-spline.
	  vector<int> KV = (*ixy==YFIXED) ? kvec[YFIXED] : kall[YFIXED];
	  vector<int> KU = (*ixy==XFIXED) ? kvec[XFIXED] : kall[XFIXED];
	  for (unsigned int ikv=0; ikv!=KV.size(); ++ikv)
	    { // Loop: over v-knots
	      for (unsigned int iku=0; iku!=KU.size(); ++iku)
		{// Loop over u-knots
		  int nmb = kmul[*ixy][(*ixy==XFIXED) ? iku : ikv];
		  for (int nu=0; nu<=nmb; ++nu)
		    { // Loop: over multiplicity
		      MeshRectangle tmprec;
		      tmprec.dir = *ixy;
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

//==============================================================================
// Given a set of non-peelable B-splines, check if they can be combined in
// linear dependence relations
//
void LinDepUtils::checkOverloaded(int minNmb, vector<LRBSpline2D*>& funs,
				  vector<vector<LRBSpline2D*> >& lindep)
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
      of << funs[ki]->umin() << " " << funs[ki]->vmin() << " 0 ";
      of << funs[ki]->umax() << " " << funs[ki]->vmin() << " 0" << std::endl;
      of << funs[ki]->umin() << " " << funs[ki]->vmin() << " 0 ";
      of << funs[ki]->umin() << " " << funs[ki]->vmax() << " 0" << std::endl;
      of << funs[ki]->umax() << " " << funs[ki]->vmin() << " 0 ";
      of << funs[ki]->umax() << " " << funs[ki]->vmax() << " 0" << std::endl;
      of << funs[ki]->umin() << " " << funs[ki]->vmax() << " 0 ";
      of << funs[ki]->umax() << " " << funs[ki]->vmax() << " 0" << std::endl;
    }
#endif
  for (size_t ki=0; ki<funs.size(); ++ki)
    {
      // Check for previous identification
      size_t kj, kh;
      for (kj=0; kj<lindep.size(); ++kj)
	{
	  for (kh=0; kh<lindep[kj].size(); ++kh)
	    if (lindep[kj][kh] == funs[ki])
	      break;
	  if (kh < lindep[kj].size())
	    break;
	}
      if (kj < lindep.size())
	continue;
      
      // Find all B-splines with support completely inside the support of
      // this B-spline and count how many are overloaded
      vector<LRBSpline2D*> inside;
      inside.push_back(funs[ki]);
      for (auto it1=funs[ki]->supportedElementBegin();
	   it1!=funs[ki]->supportedElementEnd(); ++it1)
	{
	  for (auto it2=(*it1)->supportBegin(); it2!=(*it1)->supportEnd(); ++it2)
	    {
	      if (*it2 == funs[ki])
		continue;
	      if (funs[ki]->covers(*it2))
		{
		  // Check for overload
		  for (kh=0; kh<funs.size(); ++kh)
		    if ((*it2) == funs[kh])
		      break;
		  if (kh < funs.size())
		    {
		      // Check for multiplicity
		      size_t kr;
		      for (kr=0; kr<inside.size(); ++kr)
			if ((*it2) == inside[kr])
			  break;
		      if (kr == inside.size())
			inside.push_back(*it2);
		    }
		}
	    }
	}
      if ((int)inside.size() >= minNmb)
	lindep.push_back(inside);
    }
  
  // Remove internal cases
  for (size_t kr=0; kr<lindep.size(); )
    {
      size_t kh;
      for (kh=0; kh<lindep.size(); ++kh)
	{
	  if (kh == kr)
	    continue;
	  size_t kv;
	  for (kv=1; kv<lindep[kh].size(); ++kv)
	    if (lindep[kr][0] == lindep[kh][kv])
	      break;
	  if (kv < lindep[kh].size())
	    break;
	}
      if (kh < lindep.size())
	lindep.erase(lindep.begin()+kr);
      else
	{
	  ++kr;
	}
    }
  
#ifdef DEBUG
  std::cout << "Number of linear dependency sources: " << lindep.size() << std::endl;
  for (size_t kj=0; kj<lindep.size(); ++kj)
    {
      std::cout << lindep[kj][0]->umin() << " " << lindep[kj][0]->umax() << " ";
      std::cout << lindep[kj][0]->vmin() << " " << lindep[kj][0]->vmax() << std::endl;
    }
#endif
}

