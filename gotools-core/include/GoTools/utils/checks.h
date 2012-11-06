#ifndef CHECKS_H_
#define CHECKS_H_

#include <algorithm>
#include <iostream>


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
  const int len1 = end_1 - begin_1;
  const int len2 = end_2 - begin_2;
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


