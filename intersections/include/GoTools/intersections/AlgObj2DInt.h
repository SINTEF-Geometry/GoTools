//===========================================================================
//                                                                           
// File: AlgObj2DInt.h                                                       
//                                                                           
// Created: Mon Jan 24 11:03:43 2005                                         
//                                                                           
// Author: Sverre Briseid <Sverre.Briseid@sintef.no>
//                                                                           
// Revision: $Id: AlgObj2DInt.h,v 1.9 2006-03-08 09:31:19 jbt Exp $
//                                                                           
// Description:
//                                                                           
//===========================================================================


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

