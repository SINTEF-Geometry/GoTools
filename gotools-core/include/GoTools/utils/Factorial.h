//===========================================================================
//                                                                           
// File: Factorial.h                                                         
//                                                                           
// Created: Fri May 27 14:37:29 2005                                         
//                                                                           
// Author: Odd A. Andersen <Odd.Andersen@sintef.no>
//                                                                           
// Revision: $Id: Factorial.h,v 1.4 2006-04-19 09:15:50 jbt Exp $
//                                                                           
// Description:
//                                                                           
//===========================================================================

#ifndef _FACTORIAL_H
#define _FACTORIAL_H


namespace Go {


/// Compile-time factorial calculations.


/// Calculate the factorial of a given integer N.
template<int N>
struct Factorial {
    enum { value = N * Factorial<N-1>::value };
};


/// \internal
template <>
struct Factorial<1> {
    enum { value = 1 };
};


/// Compute the inverse of the factorial of a given integer N, where
/// T is typically 'float' or 'double'.
template<typename T, int N>
struct InverseFactorial {
    static inline T val()
    {
	return T(1) / T(Factorial<N>::value);
    }
};
 

};
#endif // _FACTORIAL_H

