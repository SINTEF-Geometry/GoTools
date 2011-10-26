//===========================================================================
//                                                                           
// File: binom.h                                                             
//                                                                           
// Created: Fri May  4 12:55:18 2001                                         
//                                                                           
// Author: Atgeirr F Rasmussen <atgeirr@sintef.no>
//                                                                           
// Revision: $Id: binom.h,v 1.13 2005-06-09 07:29:45 oan Exp $
//                                                                           
// Description:
//                                                                           
//===========================================================================

#ifndef _BINOM_H
#define _BINOM_H

#include <stdexcept>
#include <vector>

namespace Go {


/// Computes the binomial coefficient: n! / (i! (n-i)!)
inline double binom(int n, int i)
{
    if (i < 0 || i > n)
	return 0.0;

    static std::vector<std::vector<double> > pascals_triangle(10);
    static bool first_time = true;
    if (first_time) {
	first_time = false;
	pascals_triangle.resize(1);
	pascals_triangle[0].resize(1);
	pascals_triangle[0][0] = 1;
    }
    int old_size = (int)pascals_triangle.size();
    if (old_size < n+1) {
	// We must expand the triangle
	pascals_triangle.resize(n+1);
	// Compute the terms of the new rows
	for (int nr = old_size; nr < n+1; ++nr) {
	    pascals_triangle[nr].resize(nr+1);
	    pascals_triangle[nr][0] = 1.0;
	    for (int j = 1; j < nr; ++j) {
		pascals_triangle[nr][j]
		    = pascals_triangle[nr-1][j-1] + pascals_triangle[nr-1][j];
	    }
	    pascals_triangle[nr][nr] = 1.0;
	}
    }
    return pascals_triangle[n][i];
}

/// computes n! (n factorial)
inline double factorial(int n)
{
    double res = 1;
    for (int i = 2; i <= n; ++i) {
 	res *= i;
    }
    return res;
}

/// computes the trinomial coefficient: n! / (i! j! (n-i-j)!)
inline double trinomial(int n, int i, int j)
{
    if (i < 0 || i > n || j < 0 || j > n)
	return 0;

    return binom(n, i) * binom(n-i, j);
}


/// computes the quadrinomial coefficient: n! / (i! j! k! (n-i-j-k)!)
inline double quadrinomial(int n, int i, int j, int k)
{
    if (i < 0 || i > n || j < 0 || j > n || k < 0 || k > n)
	return 0;

    return binom(n, i) * binom(n-i, j) * binom(n-i-j, k);
}

}; // end Go

#endif // _BINOM_H

