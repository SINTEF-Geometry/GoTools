/*****************************************************************************/
/*                                                                           */
/*                                                                           */
/* (c) Copyright 1998 by                                                     */
/*     SINTEF, Oslo, Norway                                                  */
/*     All rights reserved. See the copyright.h for more details.            */
/*                                                                           */
/*****************************************************************************/

#ifndef PRVEC_H
#define PRVEC_H

#include <vector>
#include <iostream>

/*<PrVec-syntax: */

/** PrVec -  This class represents a vector.
 *
 */
class PrVec
{
protected:

    std::vector<double> a_;

public:
    /// Default constructor
    PrVec()
    {}
    /// Constructor. Constructs a vector with n elements initialized
    /// by fill_with.
    PrVec(int n, double fill_with = 0.0)
	: a_(n, fill_with)
    {}

    /// Constructor generating av vector from a range of elements 
    /// specified with two iterators.
    template <typename InputIterator>
    PrVec(InputIterator begin, InputIterator end)
	: a_(begin, end)
    {}

    /// Change size of vector
    void redim(int n, double fill_with = 0.0);

    /// Query size of vector
    int size() const {return (int)a_.size();}

    /// Element access
          double& operator () (int i) {return a_[i];}
    /// Element access
    const double& operator () (int i) const {return a_[i];}
    /// Element access
          double& operator [] (int i) {return a_[i];}
    /// Element access
    const double& operator [] (int i) const {return a_[i];}

    /// Compute inner product with another vector of the same length.
    double inner(const PrVec& x);

    /// Read vector elements from stream 'is'.
    void read(std::istream& is);

    /// Write vector elements, separated with spaces, to stream 'os'
    void print(std::ostream& os);
};


/*>PrVec-syntax: */

/*Class:PrVec

Name:              PrVec
Syntax:	           @PrVec-syntax
Keywords:
Description:       
Member functions:
                   
Constructors:
Files:
Example:

See also:
Developed by:      SINTEF Applied Mathematics, Oslo, Norway
Author:	           Michael Floater, SINTEF
Date:              Dec. 98
*/

#endif // PRVEC_H





