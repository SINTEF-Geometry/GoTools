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

#ifndef _BERNSTEINTRIANGULARPOLY_H
#define _BERNSTEINTRIANGULARPOLY_H


#include "GoTools/utils/Array.h"
#include "GoTools/utils/ScratchVect.h"
#include "GoTools/utils/errormacros.h"
#include <vector>
#include <iostream>


namespace Go {


/**
 * Class that implements Bernstein polynomials on a triangle.
 * The formula for a triangular Bernstein polynomial of degree n is
 * \f[
 *    p(\beta_0, \beta_1, \beta_2) = \sum_{i + j + k = n}
 *       \frac{n!}{i!j!k!} \beta_0^i \beta_1^j \beta_2^k c_{(i,j,k)}.
 * \f]
 * The coefficients \f$c_{ijk}\f$ are stored in the following order:
 * \f$ c_{(n,0,0)}, c_{(n-1,1,0)}, c_{(n-2,2,0)}, \ldots, c_{(0,n,0)},
 * c_{(n-1,0,1)}, c_{(n-2, 1, 1)}, \ldots, c_{(0,0,n)} \f$.
 */

class BernsteinTriangularPoly {
public:
    /// Default constructor.
    BernsteinTriangularPoly() : deg_(-1) { }
    /// Constructor.
    /// \param deg degree of the polynomial
    /// \param coefs vector of Bernstein coefficients
    BernsteinTriangularPoly(int deg, const std::vector<double>& coefs)
	: deg_(deg), coefs_(coefs) { }
    /// Constructor.
    /// \param deg degree of the polynomial
    /// \param begin iterator to start of vector of Bernstein
    /// coefficients
    /// \param end iterator to end of vector of Bernstein coefficients
    template <typename ForwardIterator>
    BernsteinTriangularPoly(int deg,
			    ForwardIterator begin, ForwardIterator end)
	: deg_(deg), coefs_(begin, end) { }
    /// Constructor for a constant polynomial.
    /// \param coef the constant value of the polynomial
    explicit BernsteinTriangularPoly(double coef)
	: deg_(0), coefs_(1, coef) { }

    /// Get the degree
    /// \return the degree of the polynomial
    int degree() const
    { return deg_; }

    /// Access function for coefficients
    const double operator[] (int num) const
    { return coefs_[num]; }
    /// Access function for coefficients
    double& operator[] (int num)
    { return coefs_[num]; }

    /// Evaluation operator, based on the triangular de Casteljau
    /// algorithm. Adapted from 'tri_decas' in Farin: "Curves and
    /// Surfaces for CAGD".
    /// \param u parameter point given in barycentric coordinates
    /// \return the value of the polynomial at u
    template <typename T>
    T operator() (const Array<T, 3>& u) const
    {
        // Adapted from the de Casteljau algorithm 'tri_decas' in Farin,
        // "Curves and Surfaces for CAGD".

	ASSERT(deg_ >= 0);

        int sz = (int)coefs_.size();
	std::vector<T> tmp(sz);
	for (int i = 0; i < sz; ++i) {
	    tmp[i] = T(coefs_[i]);
	}

	for (int r = 1; r <= deg_; ++r) {
	    int m = -1;
	    for (int i = 0; i <= deg_-r; ++i) {
	        for (int l = 0; l <= i; ++l) {
		    m++;
		    tmp[m] = u[0] * tmp[m]
		        + u[1] * tmp[m + 1 + i]
		        + u[2] * tmp[m + 2 + i];
		}
	    }
	}

	return tmp[0];
    }

    /// Calculates the norm of the polynomial, defined as the sum of
    /// absolute values of the coefficients divided by the number of
    /// coefficients. This norm is not the same as the L1 norm of the
    /// polynomial over the triangle, but is an upper bound for it.
    double norm() const;
    /// Normalizes the polynomial by dividing all coefficients with
    /// norm(), which is related to (but not identical to) the L1 norm
    /// of the polynomial. \see norm()
    void normalize();

    /// Calculates the der'th derivative in the d-direction in terms of a new
    /// BernsteinTriangularPoly. The direction is given in barycentric
    /// coordinates, so it is a vector whose elements sum to zero.
    /// \param der order of derivative requested
    /// \param d direction of derivative
    /// \retval btp the derivative polynomial
    void deriv(int der, const Vector3D& d,
	       BernsteinTriangularPoly& btp) const;

    /// Evaluates the blossom of the polynomial. The blossom of the
    /// polynomial \f$p\f$ of degree n is a multi-affine, symmetric
    /// n-variate function \f$B_p\f$ that satisfies
    /// \f[
    ///   B_p(b, \ldots, b) = p(b).
    /// \f]
    /// \param uvec the vector of values at which the blossom is
    /// evaluated, each value is a domain point given in barycentric
    /// coordinates.
    template <typename T>
    T blossom(const std::vector<Array<T, 3> >& uvec) const
    {
        // This routine is a natural extension of operator().

	ASSERT(deg_ >= 0);

        int sz = (int)coefs_.size();
	std::vector<T> tmp(sz);
	for (int i = 0; i < sz; ++i)
	    tmp[i] = T(coefs_[i]);

	for (int r = 1; r <= deg_; ++r) {
	    const T& u = uvec[r-1][0];
	    const T& v = uvec[r-1][1];
	    const T& w = uvec[r-1][2];
	    int m = -1;
	    for (int i = 0; i <= deg_-r; ++i) {
	        for (int l = 0; l <= i; ++l) {
		    m++;
		    tmp[m] = u * tmp[m]
		        + v * tmp[m + 1 + i]
		        + w * tmp[m + 2 + i];
		}
	    }
	}

	return tmp[0];
    }

    /// Multiplication with another polynomial
    BernsteinTriangularPoly& operator*= (const BernsteinTriangularPoly& poly);
    /// Multiplication with a scalar
    BernsteinTriangularPoly& operator*= (double c);

    /// Addition with another polynomial
    BernsteinTriangularPoly& operator+= (const BernsteinTriangularPoly& poly);
    /// Addition with a scalar
    BernsteinTriangularPoly& operator+= (double c);

    /// Subtraction with another polynomial
    BernsteinTriangularPoly& operator-= (const BernsteinTriangularPoly& poly);
    /// Subtraction with a scalar
    BernsteinTriangularPoly& operator-= (double c);

    /// Division with a scalar
    BernsteinTriangularPoly& operator/= (double c);

    /// Read from an input stream
    void read(std::istream& is);
    /// Write to an output stream
    void write(std::ostream& os) const;

private:
    int deg_;
    //typedef ScratchVect<double, 21> VecType;
    typedef std::vector<double> VecType;
    VecType coefs_;
};


/// Multiplication of two polynomials
inline BernsteinTriangularPoly operator* (const BernsteinTriangularPoly& p1,
					  const BernsteinTriangularPoly& p2)
{
    BernsteinTriangularPoly tmp = p1;
    tmp *= p2;
    return tmp;
}


/// Multiplication of a polynomial with a scalar
inline BernsteinTriangularPoly operator* (const BernsteinTriangularPoly& p,
					  double c)
{
    BernsteinTriangularPoly tmp = p;
    tmp *= c;
    return tmp;
}


/// Multiplication of a scalar with a polynomial
inline BernsteinTriangularPoly operator* (double c,
					  const BernsteinTriangularPoly& p)
{
    BernsteinTriangularPoly tmp = p;
    tmp *= c;
    return tmp;
}


/// Addition of two polynomials
inline BernsteinTriangularPoly operator+ (const BernsteinTriangularPoly& p1,
					  const BernsteinTriangularPoly& p2)
{
    BernsteinTriangularPoly tmp = p1;
    tmp += p2;
    return tmp;
}


/// Addition of a scalar with a polynomial
inline BernsteinTriangularPoly operator+ (double c,
					  const BernsteinTriangularPoly& p)
{
    BernsteinTriangularPoly tmp = p;
    tmp += c;
    return tmp;
}


/// Addition of a polynomial with a scalar
inline BernsteinTriangularPoly operator+ (const BernsteinTriangularPoly& p,
					  double c)
{
    BernsteinTriangularPoly tmp = p;
    tmp += c;
    return tmp;
}


/// Subtraction of two polynomials
inline BernsteinTriangularPoly operator- (const BernsteinTriangularPoly& p1,
					  const BernsteinTriangularPoly& p2)
{
    BernsteinTriangularPoly tmp = p1;
    tmp -= p2;
    return tmp;
}


/// Subtraction of a polynomial from a scalar
inline BernsteinTriangularPoly operator- (double c,
					  const BernsteinTriangularPoly& p)
{
    BernsteinTriangularPoly tmp(c);
    tmp -= p;
    return tmp;
}


/// Subtraction of a scalar from a polynomial
inline BernsteinTriangularPoly operator- (const BernsteinTriangularPoly& p,
					  double c)
{
    BernsteinTriangularPoly tmp = p;
    tmp -= c;
    return tmp;
}


/// Division of a polynomial with a scalar
inline BernsteinTriangularPoly operator/ (const BernsteinTriangularPoly& p,
					  double c)
{
    BernsteinTriangularPoly tmp = p;
    tmp /= c;
    return tmp;
}


/// Read BernsteinTriangularPoly from input stream
inline std::istream& operator >> (std::istream& is,
				  Go::BernsteinTriangularPoly& p)
{
    p.read(is);
    return is;
}

/// Write BernsteinTriangularPoly to output stream
inline std::ostream& operator << (std::ostream& os,
				  const Go::BernsteinTriangularPoly& p)
{
    p.write(os);
    return os;
}


} // namespace Go



#endif // _BERNSTEINTRIANGULARPOLY_H

