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

#ifndef _TRACECONTOURSTYPEDEFS_H
#define _TRACECONTOURSTYPEDEFS_H

#include <vector>
#include <functional>
#include <algorithm>
#include "GoTools/geometry/SplineCurve.h"

namespace Go
{
  /// Common typedefs used both by SSurfTraceIsocontours and LRTraceIsocontours
  
  using CurvePtr = std::shared_ptr<const SplineCurve>;

  // The first curve pointer of the pair represents a 2D curve in the parameter
  // domain of the investigated LR spline function.  The second curve pointer
  // represents the corresponding 3D curve, if requested.
  using CurveVec = std::vector<std::pair<CurvePtr, CurvePtr>>;


  // ----------------------------------------------------------------------------
  // The following utility template function ought not to have any overhead in
  // optimized mode, due to the 'move' semantics of STL containers.
  template<typename Container, typename R, typename T>
  std::vector<R> apply_transform(const Container& input,
				 const std::function<R(T)>& trans_fun) 
  // ----------------------------------------------------------------------------
  {
    std::vector<R> tmp(input.size());
    std::transform(input.begin(), input.end(), tmp.begin(), trans_fun);
    return tmp;
  }



  // ----------------------------------------------------------------------------
  template<typename Iterator, typename Arg, typename Val>
  std::vector<Val> apply_transform(const Iterator& i1,
				   const Iterator& i2, const std::function<Val(Arg)>& fun)
    // ----------------------------------------------------------------------------
  {
    const size_t numel = size_t(i2-i1);
    std::vector<Val> result(numel);
    std::transform(i1, i2, result.begin(), fun);
    return result;
  }

  // ----------------------------------------------------------------------------
  template<typename T>
  std::vector<T> insert_back(const std::vector<T>& v, T elem)
  // ----------------------------------------------------------------------------
  {
    std::vector<T> result(v);
    result.push_back(elem);
    return result;
  }

  // ----------------------------------------------------------------------------
  template<typename T>
  std::vector<T> reverse_vec(std::vector<T> v)
  // ----------------------------------------------------------------------------
  {
    reverse(v.begin(), v.end()); return v;
  }

  // ----------------------------------------------------------------------------
  template<typename T>
  std::vector<T> merge_vec(std::vector<T> v1, const std::vector<T>& v2)
  // ----------------------------------------------------------------------------
  {
    v1.insert(v1.end(), v2.begin(), v2.end()); return v1;
  }
};

#endif
