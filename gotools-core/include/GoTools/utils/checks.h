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

#ifndef CHECKS_H_
#define CHECKS_H_

#include <algorithm>
#include <iostream>
#include <functional>


//==============================================================================
// THIS FILE CONTAINS VARIOUS GENERIC RANGE CHECK FUNCTIONS (WHETHER A RANGE IS
// INCREASING/DECREASING, WHETHER INTERVALS OVERLAP, AND WHETHER ONE RANGE IS
// 'SMALLER' THAN ANOTHER
//==============================================================================

namespace Go
{

// =============================================================================
template<typename Iterator> bool strictly_increasing(Iterator begin, Iterator end)
// =============================================================================
{
  return is_sorted(begin, end, std::less_equal<decltype(*begin)>());
  //return is_sorted(begin, end, [](int next, int cur) {return next <= cur;});
}

// =============================================================================
template<typename Array> bool strictly_increasing(const Array& A)
// =============================================================================
{
  return strictly_increasing(A.begin(), A.end());
}


// =============================================================================
template<typename Iterator> bool weakly_increasing(Iterator begin, Iterator end)
// =============================================================================
{
  return is_sorted(begin, end, std::less<decltype(*begin)>());
}

// =============================================================================
template<typename Array> bool weakly_increasing(const Array& A)
// =============================================================================
{
  return weakly_increasing(A.begin(), A.end());
}

// =============================================================================
template<typename Iterator> bool strictly_decreasing(Iterator begin, Iterator end)
// =============================================================================
{
  return is_sorted(begin, end, std::greater_equal<decltype(*begin)>());
};

// =============================================================================
template<typename Array> bool strictly_decreasing(const Array& A)
// =============================================================================
{
  return strictly_decreasing(A.begin(), A.end());
};

// =============================================================================
template<typename Iterator> bool weakly_decreasing(Iterator begin, Iterator end)
// =============================================================================
{
  return is_sorted(begin, end, std::greater<decltype(*begin)>());
};

// =============================================================================
template<typename Array> bool weakly_decreasing(const Array& A)
// =============================================================================
{
  return weakly_decreasing(A.begin(), A.end());
};

// =============================================================================
template<typename ValueType>
bool nondecreasing(ValueType a, ValueType b, ValueType c) 
// =============================================================================
{
  return a <= b && b <= c;
}

// =============================================================================
template<typename ValueType>
bool nonincreasing(ValueType a, ValueType b, ValueType c) 
// =============================================================================
{
  return a >= b && b >= c;
}

// =============================================================================
template<typename ValueType>
bool strictly_decreasing(ValueType a, ValueType b, ValueType c)
// =============================================================================
{
  return a > b && b > c;
}

// =============================================================================
template<typename ValueType>
bool strictly_increasing(ValueType a, ValueType b, ValueType c)
// =============================================================================
{
  return a < b && b < c;
}

// =============================================================================
// checks if the two intervals [front1, back1] and [front2, back2] has nonzero
// overlap
template<typename ValueType>
bool interval_overlap(ValueType front1, ValueType back1, 
		      ValueType front2, ValueType back2)
// =============================================================================
{
  return ( (back1 - front2) * (back2 - front1) > 0 );
}

// =============================================================================
// Return -1 if first sequence "smaller" than second.  Return 1 if second sequence
// "smaller" than first.  Return 0 if they are equal.
template<typename Iterator> int compare_seq(Iterator begin_1, Iterator end_1, 
					    Iterator begin_2, Iterator end_2)
// =============================================================================
{
  const int len1 = (int)(end_1 - begin_1);
  const int len2 = (int)(end_2 - begin_2);
  if (len1 != len2) return (len1 < len2) ? -1 : 1;

  // sequences are of equal length.  Compare elements
  for ( ; begin_1 != end_1; ++begin_1, ++ begin_2) {
    if (*begin_1 != *begin_2) return (*begin_1 < *begin_2) ? -1 : 1;
  }
  // sequences are equal
  return 0;
}



}; // end namespace Go

#endif


