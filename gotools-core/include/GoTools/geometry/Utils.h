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

#ifndef _UTILS_H
#define _UTILS_H


#include "GoTools/geometry/SplineCurve.h"
#include "GoTools/utils/errormacros.h"
#include <math.h>
#include <algorithm>
#include <ctype.h>


namespace Go
{


    /// Iterator traits classes are provided for the Go namespace, so that 
    /// we won't need to include \code <iterator> \endcode or similar headers.

template <class Iterator>
struct go_iterator_traits {
  typedef typename Iterator::value_type        value_type;
  typedef typename Iterator::difference_type   difference_type;
  typedef typename Iterator::pointer           pointer;
  typedef typename Iterator::reference         reference;
};

template <class T>
struct go_iterator_traits<T*> {
  typedef T                          value_type;
  typedef int                        difference_type;
  typedef T*                         pointer;
  typedef T&                         reference;
};


template <class T>
struct go_iterator_traits<const T*> {
  typedef T                          value_type;
  typedef int                        difference_type;
  typedef const T*                   pointer;
  typedef const T&                   reference;
};

/// Namespace for some utility functions
namespace Utils
{
/// sum finds the sum of the elements
template <typename ForwardIterator>
inline typename go_iterator_traits<ForwardIterator>::value_type 
sum(ForwardIterator first,
    ForwardIterator last)
{
    typename go_iterator_traits<ForwardIterator>::value_type sum = 0;
    for (; first != last; ++first)
        sum += *first;
    return sum;
}

/// sum_squared finds the squared sum of the elements
template <typename ForwardIterator>
inline typename go_iterator_traits<ForwardIterator>::value_type
sum_squared(ForwardIterator first,
            ForwardIterator last)
{
    typename go_iterator_traits<ForwardIterator>::value_type sum = 0;
    for (; first != last; ++first)
        sum += (*first)*(*first);
    return sum;
}

/// distance_squared
template <typename ForwardIterator>
inline typename go_iterator_traits<ForwardIterator>::value_type
distance_squared(ForwardIterator first1,
                 ForwardIterator last1,
                 ForwardIterator first2)
{
    typename go_iterator_traits<ForwardIterator>::value_type sum = 0;
    for (; first1 != last1; ++first1, ++first2)
        sum += (*first1 - *first2)*(*first1 - *first2);
    return sum;
}

/// normalize makes the length of a vector 1.0
template <typename ForwardIterator>
inline void
normalize(ForwardIterator first,
          ForwardIterator last)
{
    typename go_iterator_traits<ForwardIterator>::value_type d
        = sqrt(sum_squared(first, last));
    d = 1.0/d;
    for (; first != last; ++first)
        (*first) *= d;
}

/// inner product
template <typename ForwardIterator>
inline typename go_iterator_traits<ForwardIterator>::value_type
inner(ForwardIterator first,
      ForwardIterator last,
      ForwardIterator second)
{
    typename go_iterator_traits<ForwardIterator>::value_type sum = 0;
    for (; first != last; ++first, ++second)
        sum += (*first)*(*second);
    return sum;
}

/// eat white space
template <typename InputStream>
inline InputStream& eatwhite(InputStream& is)
{
    char c;
    while (is.get(c)) {
        if (isspace(c)==0) {
            is.putback(c);
            break;
        }
    }
    return is;
}

} // End of namespace Utils

} // End of namespace Go


#endif
