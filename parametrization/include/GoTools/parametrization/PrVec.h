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





