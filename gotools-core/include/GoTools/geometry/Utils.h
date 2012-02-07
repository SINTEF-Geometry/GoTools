// $Id: Utils.h,v 1.19 2007-12-04 16:12:01 jbt Exp $


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

#if defined(_MSC_VER) && (_MSC_VER < 1600)

    /// sum finds the sum of the elements
    inline double
	sum(double* first,
	    double* last)
	{
	    double sum = 0.0;
	    for (; first != last; ++first)
		sum += *first;
	    return sum;
	}

    /// sum_squared finds the squared sum of the elements
    inline double
	sum_squared(double* first,
		    double* last)
	{
	    double sum = 0.0;
	    for (; first != last; ++first)
		sum += (*first)*(*first);
	    return sum;
	}
    /// distance_squared
    inline double
	distance_squared(const double* first1,
			 const double* last1,
			 const double* first2)
	{
	    double sum = 0;
	    for (; first1 != last1; ++first1, ++first2)
		sum += (*first1 - *first2)*(*first1 - *first2);
	    return sum;
	}    

    /// normalize makes the length of a vector 1.0
    inline void
	normalize(double* first,
		  double* last)
	{
	    double d
		= sqrt(sum_squared(first, last));
	    d = 1.0/d;
	    for (; first != last; ++first)
		(*first) *= d;
	}

    /// inner product
    inline double
	inner(double* first,
	      double* last,
	      double* second)
	{
	    double sum = 0;
	    for (; first != last; ++first, ++second)
		sum += (*first)*(*second);
	    return sum;
	}
#else

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
#endif // _MSC_VER < 1600

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

} // End of namespace Go


#endif
