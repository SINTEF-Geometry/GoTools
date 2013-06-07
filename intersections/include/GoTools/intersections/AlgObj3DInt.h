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

#ifndef _ALGOBJ3DINT_H
#define _ALGOBJ3DINT_H


#include "GoTools/intersections/AlgObjectInt.h"
#include "GoTools/utils/Array.h"
#include "GoTools/implicitization/BernsteinTetrahedralPoly.h"
#include "GoTools/utils/BaryCoordSystem.h"
#include <vector>


namespace Go {


/// Struct that represents a monomial in three variables. By this we
/// mean a term in a polynomial of the form \f$ax^iy^jx^k\f$. This
/// struct is used by AlgObj3DInt.

struct Alg3DElem {

    /// The factor \a a in \f$ax^iy^jx^k\f$
    double factor_;

    /// The degrees \a i, \a j and \a k in \f$ax^iy^jx^k\f$
    Array<int, 3> degrees_; // (deg_x, deg_y, deg_z).

    /// Constructor
    Alg3DElem(double factor, int degree_x, int degree_y, int degree_z)
	: factor_(factor), 
	  degrees_(Array<int, 3>(degree_x, degree_y, degree_z))
    {}
};


/// Class for 3-dimensional algebraic intersection objects. Supports
/// two different representations: Polynomials on power basis, or
/// Bernstein polynomials on a tetrahedron.

class AlgObj3DInt : public AlgObjectInt {
public:
    /// Constructor.
    /// \param degree the total degree of the algebraic expression,
    /// i.e. the maximum sum of exponents of a term.
    AlgObj3DInt(int degree);

    /// Constructor.
    /// \param terms the terms in the algebraic expression,
    /// i.e. elements on the form \f$a_{ijk} x^i y^j z^k\f$.
    AlgObj3DInt(const std::vector<Alg3DElem>& terms);

    /// Constructor.
    /// We define the algebraic expression using another
    /// representation than the standard power basis formulation. This
    /// will (typically) be the result from an approximative
    /// implicitization of a spline surface.
    /// \param implicit the implicit object.
    /// \param bc the barycentric coordinate system for the
    /// representation.
    AlgObj3DInt(const BernsteinTetrahedralPoly& implicit,
		const BaryCoordSystem3D& bc);

    /// Destructor.
    virtual ~AlgObj3DInt();

    /// Get the number of terms in the algebraic object.
    /// \return The number of terms in the object.
    int numTerms()
    { return (int)terms_.size(); }

    /// Get the degree of the algebraic object
    /// \return the degree
    int degree()
    { return degree_; }

    /// Get the corresponding term from the algebraic object.
    /// \param index the index of the term in question. Indexing
    /// starts at 0.
    Alg3DElem term(int index);

    /// Get the implicit representation of the object.
    /// \param impl the implicit representation.
    /// \param bc the corresponding coordinate system.
    void getImplicit(BernsteinTetrahedralPoly& impl, BaryCoordSystem3D& bc)
    {
	impl = implicit_;
	bc = bc_;
    }

    /// Verify whether we are using a standard power basis
    /// representation for the object.
    /// \return True if we are using a standard power basis
    /// representation.  For typical algebraic objects like spheres
    /// this will be the case, but not for approximative
    /// implicitizations of spline surfaces.
    bool usingPowerBasis()
    { return power_basis_; }


protected:
    int degree_; // Total degree, i.e. the highest sum of the degrees
		 // of a term in the polynomial.

    std::vector<Alg3DElem> terms_;
    bool power_basis_; // The other option is Bernstein basis &
		       // barycentric coordinates.

    // @@sbr Not sure which to use!
    // Should the approximative impclicitization take place inside
    // this function?  Would assume that it was to be calculated on
    // the outside and only stored here.  Could then use class in
    // Implicitization, which depends on that structure not to change.
    // But tu implement easy support for converting from power to
    // Bernstein basis this approach makes sense.
    BernsteinTetrahedralPoly implicit_; // deg &
					// (deg+1)*(deg+2)*(deg+3)/6
					// basis elements.
    BaryCoordSystem3D bc_; // 4 3D-corners. Requires boundingbox for
			   // alg obj.

private:

//     // We need a utility function to change from power basis to
//     // Bernstein, using the format in implicitization.
//     BernsteinTetrahedralPoly powerToBernstein();

//     // We convert x to the Bernstein basis.
//     SplineSurface unitSf();

};


} // namespace Go


#endif // _ALGOBJ3DINT_H

