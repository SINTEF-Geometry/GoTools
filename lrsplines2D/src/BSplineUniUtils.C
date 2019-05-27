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

#include "GoTools/utils/checks.h"
#include "GoTools/lrsplines2D/BSplineUniUtils.h"

//------------------------------------------------------------------------------

using std::vector;
using std::unique_ptr;

namespace Go
{
//==============================================================================
  bool BSplineUniUtils::identify_bsplineuni(const BSplineUniLR* bspline, 
					    vector<unique_ptr<BSplineUniLR> >& bspline_vec,
					    int& ix)
//==============================================================================
  {
    if (bspline_vec.size() == 0)
      return false;

    if (ix < 0 || ix >= (int)bspline_vec.size())
      ix = (int)bspline_vec.size()/2;

    int comp = (*bspline < *bspline_vec[ix]);    
    if (comp == 0)
      return true;

    // Search nearby
    if (ix > 0 && *bspline == *bspline_vec[ix-1])
      {
	--ix;
	return true;
      }
    if (ix < (int)bspline_vec.size()-1 && *bspline == *bspline_vec[ix+1])
      {
	++ix;
	return true;
      }

    // Search first 
    if (*bspline == *bspline_vec[0])
      {
	ix = 0;
	return true;
      }


    // Midpoint search
    int ix1 = (comp < 0) ? 0 : ix+1;
    int ix2 = (comp < 0) ? ix-1 : (int)bspline_vec.size()-1;
     while (ix2 > ix1 && (ix != ix1 || ix2-ix1>1))
      {
	ix = (ix1 + ix2)/2;
	comp = (*bspline < *bspline_vec[ix]);
	if (comp == 0)
	  return true;
	if (comp > 0)
	  ix1 = ix;
	if (comp < 0)
	  ix2 = ix-1;
      }
    ix2 = std::min(std::max(ix2, 0), (int)bspline_vec.size()-1);
    if (*bspline == *bspline_vec[ix2])
      {
	ix = ix2;
	return true;
      }
    return false;
  }

//==============================================================================
  bool BSplineUniUtils::identify_bsplineuni(vector<int>::const_iterator start,
					    vector<int>::const_iterator end,
					    vector<unique_ptr<BSplineUniLR> >& bspline_vec,
					    int& ix)
//==============================================================================
  {
    if (bspline_vec.size() == 0)
      return false;

    if (ix < 0 || ix >= (int)bspline_vec.size())
      ix = (int)bspline_vec.size()/2;

    int comp = compare_seq(start, end, bspline_vec[ix]->kvec_const().begin(),
			   bspline_vec[ix]->kvec_const().end());
    if (comp == 0)
      return true;

    // Search nearby
    if (ix > 0 && 
	compare_seq(start, end, bspline_vec[ix-1]->kvec_const().begin(),
		    bspline_vec[ix-1]->kvec_const().end()) == 0)
      {
	--ix;
	return true;
      }
    if (ix < (int)bspline_vec.size()-1 && 
	compare_seq(start, end, bspline_vec[ix+1]->kvec_const().begin(),
		    bspline_vec[ix+1]->kvec_const().end()) == 0)
      {
	++ix;
	return true;
      }

    // Search first 
    if (compare_seq(start, end, bspline_vec[0]->kvec_const().begin(),
		    bspline_vec[0]->kvec_const().end()) == 0)
      {
	ix = 0;
	return true;
      }


    // Midpoint search
    int ix1 = (comp < 0) ? 0 : ix+1;
    int ix2 = (comp < 0) ? ix-1 : (int)bspline_vec.size()-1;
    while (ix2 > ix1 && (ix != ix1 || ix2-ix1>1))
      {
	ix = (ix1 + ix2)/2;
	comp = compare_seq(start, end, bspline_vec[ix]->kvec_const().begin(),
			       bspline_vec[ix]->kvec_const().end());
	if (comp == 0)
	  return true;
	if (comp > 0)
	  ix1 = ix;
	if (comp < 0)
	  ix2 = ix-1;
      }
    ix2 = std::min(std::max(ix2, 0), (int)bspline_vec.size()-1);
    if (compare_seq(start, end, bspline_vec[ix2]->kvec_const().begin(),
		    bspline_vec[ix2]->kvec_const().end()) == 0)
      {
	ix = ix2;
	return true;
      }
    return false;
  }

//==============================================================================
  void BSplineUniUtils::insert_univariate(vector<unique_ptr<BSplineUniLR> >& bspline_vec,
					  BSplineUniLR* bspline, 
					  int& ix)
//==============================================================================
  {
    if (bspline_vec.size() == 0)
      {
	bspline_vec.push_back(unique_ptr<BSplineUniLR>(bspline));
	ix = 0;
      }
    else
      {
	int ix2;
	int comp = 0;
	for (ix2=std::min((int)bspline_vec.size()-1,ix+1); ix2>=0; --ix2)
	  {
	    comp = ((*bspline) < (*bspline_vec[ix2]));
	    if (comp >= 0)
	      {
		ix = ix2+1;
		break;
	      }
	  }
	if (comp > 0)
	  bspline_vec.insert(bspline_vec.begin()+ix, unique_ptr<BSplineUniLR>(bspline));
	else if (ix2 < 0 && comp < 0)
	  {
	    ix = 0;
	    bspline_vec.insert(bspline_vec.begin(), unique_ptr<BSplineUniLR>(bspline));
	  }
	else if (comp != 0)
	  {
	    for (ix2=ix; ix2<(int)bspline_vec.size(); ++ix2)
	      {
		comp = ((*bspline) < (*bspline_vec[ix2]));
		if (comp >= 0)
		  {
		    ix = ix2+1;
		    break;
		  }
		if (comp > 0)
		  bspline_vec.insert(bspline_vec.begin()+ix, unique_ptr<BSplineUniLR>(bspline));
	      }
	  }
      }
  }

}; // end namespace Go
