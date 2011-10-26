//==========================================================================
//                                                                          
// File: Binomial.h                                                          
//                                                                          
// Created: Mon Nov  4 13:18:19 2002                                         
//                                                                          
// Author: Jan B. Thomassen <jbt@sintef.no>
//                                                                          
// Revision: $Id: Binomial.h,v 1.11 2006-12-13 12:49:28 jbt Exp $
//                                                                          
// Description:
//                                                                          
//==========================================================================

#ifndef _BINOMIAL_H
#define _BINOMIAL_H


#include <vector>


namespace Go {


/** 
 * Class that computes the binomial coefficients.
 */

class Binomial {
public:
    /// Useful iterator
    typedef std::vector<double>::iterator iter;

    /// Default constructor. 
    Binomial() : pascals_triangle(10)
    {
	for (int i = 0; i < 10; ++i)
	    pascals_triangle[i].reserve(10);
	pascals_triangle.resize(1);
	pascals_triangle[0].resize(1);
	pascals_triangle[0][0] = 1.0;
    }
    /// Constructor that creates an initial Pascal's triangle to some
    /// specified level.
    /// \param n the level up to which you want to make the triangle
    Binomial(int n)
    {
        init(n);
    }

    /// Evaluates a binomial coefficient. If the
    /// current table is too small, the table will be expanded
    /// accordingly. 
    /// \param n an integer
    /// \param i an integer
    /// \return \f[ {n \choose i} \f]
    double operator() (int n, int i)
    {
	if (i < 0 || i > n)
	    return 0.0;
	if (int(pascals_triangle.size()) < n+1)
	    expand(n);

	return pascals_triangle[n][i];
    }

    /// Gives an iterator to the first element of the vector containing the
    /// n'th line of Pascal's triangle. Expands the table if necessary.
    /// \param n an integer
    /// \return an iterator iter such that iter[i] = operator()(n,i).
    iter operator[] (int n)
    {
	if (int(pascals_triangle.size()) < n+1)
	    expand(n);

	return pascals_triangle[n].begin();
    }

    /// Evaluates a trinomial coefficient. This function should be
    /// considered a temporary solution.
    /// \param n an integer
    /// \param i an integer
    /// \param j an integer
    /// \return \f[\frac{n!}{i!j!(n-i-j)!}.\f]
    double trinomial(int n, int i, int j)
    {
	if (i < 0 || i > n || j < 0 || j > n)
	    return 0.0;
	if (int(pascals_triangle.size()) < n+1)
	    expand(n);

	return pascals_triangle[n][i] * pascals_triangle[n-i][j];
    }

    /// Evaluates a quadrinomial coefficient. This function should be
    /// considered a temporary solution.
    /// \param n an integer
    /// \param i an integer
    /// \param j an integer
    /// \param k an integer
    /// \return \f[\frac{n!}{i!j!k!(n-i-j-k)!.}\f]
    double quadrinomial(int n, int i, int j, int k)
    {
	if (i < 0 || i > n || j < 0 || j > n || k < 0 || k > n)
	    return 0;
	if (int(pascals_triangle.size()) < n+1)
	    expand(n);
	
	return pascals_triangle[n][i]
	    * pascals_triangle[n-i][j]
	    * pascals_triangle[n-i-j][k];
    }

private:
    void init(int n);
    void expand(int n);

    std::vector<std::vector<double> > pascals_triangle;

};


} // namespace Go


#endif // _BINOMIAL_H

