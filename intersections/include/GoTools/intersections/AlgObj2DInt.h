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

#ifndef _ALGOBJ2DINT_H
#define _ALGOBJ2DINT_H


#include "GoTools/intersections/AlgObjectInt.h"
#include "GoTools/utils/Array.h"
#include <vector>


namespace Go {


/// Struct that represents a monomial in two variables. By this we
/// mean a term in a polynomial of the form \f$ax^iy^j\f$. This struct
/// is used by AlgObj2DInt.

struct Alg2DElem {

    /// The factor \a a in \f$ax^iy^j\f$
    double factor_;

    /// The degrees \a i and \a j in \f$ax^iy^j\f$
    Array<int, 2> degrees_; // (deg_x, deg_y).

    /// Constructor
    Alg2DElem(double factor, int degree_x, int degree_y)
	: factor_(factor), degrees_(Array<int, 2>(degree_x, degree_y))
    {}
};


/// Class for 2-dimensional algebraic intersection objects.

class AlgObj2DInt : public AlgObjectInt {
public:
    /// Constructor.
    /// \param degree the total degree of the algebraic expression,
    /// i.e. the maximum sum of exponents of a factor.
    AlgObj2DInt(int degree);

    /// Constructor.
    /// \param terms the terms in the algebraic expression,
    /// i.e. elements on the form \f$a_{ij} x^i y^j\f$.
    AlgObj2DInt(const std::vector<Alg2DElem>& terms);

    /// Destructor.
    virtual ~AlgObj2DInt();

    /// Get the number of terms in the algebraic object.
    /// \return The number of terms in the object.
    int numTerms()
    { return (int)terms_.size(); }

    /// Get a term from the algebraic object.
    /// \return the term with index \a index
    Alg2DElem term(int index);

protected:

    int degree_;
    std::vector<Alg2DElem> terms_;
    bool power_basis_;

};


} // namespace Go


#endif // _ALGOBJ2DINT_H

