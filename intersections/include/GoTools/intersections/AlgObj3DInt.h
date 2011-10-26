//===========================================================================
//                                                                           
// File: AlgObj3DInt.h                                                       
//                                                                           
// Created: Mon Jan 24 11:03:55 2005                                         
//                                                                           
// Author: Sverre Briseid <Sverre.Briseid@sintef.no>
//                                                                           
// Revision: $Id: AlgObj3DInt.h,v 1.11 2006-05-04 12:19:00 jbt Exp $
//                                                                           
// Description:
//                                                                           
//===========================================================================


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

