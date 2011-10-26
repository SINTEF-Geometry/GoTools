//==========================================================================
//                                                                          
// File: BernsteinTetrahedralPoly.h                                          
//                                                                          
// Created: Tue Jan 21 09:11:36 2003                                         
//                                                                          
// Author: Jan B. Thomassen <jbt@sintef.no>
//                                                                          
// Revision: 
// $Id: BernsteinTetrahedralPoly.h,v 1.20 2007-03-13 14:43:19 jbt Exp $
//                                                                          
// Description:
//                                                                          
//==========================================================================

#ifndef _BERNSTEINTETRAHEDRALPOLY_H
#define _BERNSTEINTETRAHEDRALPOLY_H


#include "GoTools/utils/Array.h"
#include "GoTools/utils/errormacros.h"
#include <vector>
#include <iostream>


namespace Go {


class BernsteinPoly;


/**
 * Class that implements Bernstein polynomials on a tetrahedron.
 * The formula for a tetrahedral Bernstein polynomial of degree n is
 * \f[
 *    p(\beta_0, \beta_1, \beta_2, \beta_3) = \sum_{i + j + k + l = n}
 *       \frac{n!}{i!j!k!l!} \beta_0^i \beta_1^j \beta_2^k \beta_3^l
 *       c_{(i,j,k)}.
 * \f]
 * The coefficients \f$c_{ijk}\f$ are stored in the following order:
 * \f$ c_{(n,0,0,0)}, c_{(n-1,1,0,0)}, c_{(n-2,2,0,0)}, \ldots, c_{(0,n,0,0)},
 * c_{(n-1,0,1,0)}, c_{(n-2, 1, 1,0)}, \ldots, c_{(0,0,0,n)} \f$. That
 * is, a reverse lexicographical ordering with the most significant
 * index being the last one.
 */

class BernsteinTetrahedralPoly {
public:
    /// Default constructor.
    BernsteinTetrahedralPoly() : deg_(-1) { }
    /// Constructor.
    /// \param deg degree of the polynomial
    /// \param coefs vector of Bernstein coefficients
    BernsteinTetrahedralPoly(int deg, const std::vector<double>& coefs)
	: deg_(deg), coefs_(coefs) { }
    /// Constructor.
    /// \param deg degree of the polynomial
    /// \param begin iterator to start of vector of Bernstein
    /// coefficients
    /// \param end iterator to end of vector of Bernstein coefficients
    template <typename ForwardIterator>
    BernsteinTetrahedralPoly(int deg,
			     ForwardIterator begin, ForwardIterator end)
	: deg_(deg), coefs_(begin, end) { }
    /// Constructor for a constant polynomial.
    /// \param coef the constant value of the polynomial
    explicit BernsteinTetrahedralPoly(double coef)
	: deg_(0), coefs_(1, coef) { }

    /// Get the degree
    /// \return the degree of the polynomial
    int degree() const
    { return deg_; }

    /// Access function for coefficients
    const double& operator[] (int num) const
    { return coefs_[num]; }
    /// Access function for coefficients
    double& operator[] (int num)
    { return coefs_[num]; }

    /// Evaluation operator, based on the tetrahedral de Casteljau
    /// algorithm. Generalized from 'tri_decas' in Farin: "Curves and
    /// Surfaces for CAGD".
    /// \param u parameter point given in barycentric coordinates
    /// \return the value of the polynomial at u
    template <typename T>
    T operator() (const Array<T, 4>& u) const
    {
        // Generalized from the de Casteljau algorithm 'tri_decas' in
        // Farin, "Curves and Surfaces for CAGD".
        
	ASSERT(deg_ >= 0);

        int sz = (int)coefs_.size();
        std::vector<T> tmp(sz);
        for (int i = 0; i < sz; ++i)
        	tmp[i] = T(coefs_[i]);
        
        for (int r = 1; r <= deg_; ++r) {
        	int m = -1;
        	for (int i = 0; i <= deg_-r; ++i) {
        	    int k = (i+1) * (i+2) / 2;
        	    for (int j = 0; j <= i; ++j) {
        		for (int l = 0; l <= j; ++l) {
        		    m++;
        		    tmp[m] = u[0] * tmp[m]
        			+ u[1] * tmp[m + k]
        			+ u[2] * tmp[m + 1 + j + k]
        			+ u[3] * tmp[m + 2 + j + k];
        		}
        	    }
        	}
        }
        
        return tmp[0];
    }

    /// Calculates the norm of the polynomial, defined as the sum of
    /// absolute values of the coefficients divided by the number of
    /// coefficients. This norm is not the same as the L1 norm of the
    /// polynomial over the tetrahedron, but is an upper bound for it.
    double norm() const;
    /// Normalizes the polynomial by dividing all coefficients with
    /// norm(), which is related to (but not identical to) the L1 norm
    /// of the polynomial. \see norm()
    void normalize();

    /// Calculates the der'th derivative in the d-direction in terms of a new
    /// BernsteinTetrahedralPoly. The direction is given in barycentric
    /// coordinates, so it is a vector whose elements sum to zero.
    /// \param der order of derivative requested
    /// \param d direction of derivative
    /// \retval btp the derivative polynomial
    void deriv(int der, const Vector4D& d,
	       BernsteinTetrahedralPoly& btp) const;

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
    T blossom(const std::vector<Array<T, 4> >& uvec) const
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
	    const T& x = uvec[r-1][3];
	    int m = -1;
	    for (int i = 0; i <= deg_-r; ++i) {
	        int k = (i+1) * (i+2) / 2;
		for (int j = 0; j <= i; ++j) {
		    for (int l = 0; l <= j; ++l) {
		        m++;
			tmp[m] = u * tmp[m]
			    + v * tmp[m + k]
			    + w * tmp[m + 1 + j + k]
			    + x * tmp[m + 2 + j + k];
		    }
		}
	    }
	}

	return tmp[0];
    }

    /// Routine that picks out the line with endpoints a and
    /// b as a new BernsteinPoly. The new polynomial is defined
    /// on [0,1]. If p is the
    /// old bivariate polynomial and q is the new one, q is defined by
    /// \f[ q(t) = p((1-t)a_0 + tb_0, (1-t)a_1 + tv_1). \f]
    /// The function is implemented by blossoming.
    /// \param a the starting point of the line in barycentric
    /// coordinates
    /// \param b the endpoint of the line in barycentric coordinates
    BernsteinPoly pickLine(const Array<double, 4>& a,
			   const Array<double, 4>& b) const;

    /// Multiplication with another polynomial
    BernsteinTetrahedralPoly&
    operator*= (const BernsteinTetrahedralPoly& poly);
    /// Multiplication with a scalar
    BernsteinTetrahedralPoly& operator*= (double c);

    /// Addition with another polynomial
    BernsteinTetrahedralPoly&
    operator+= (const BernsteinTetrahedralPoly& poly);
    /// Addition with a scalar
    BernsteinTetrahedralPoly& operator+= (double c);

    /// Subtraction with another polynomial
    BernsteinTetrahedralPoly&
    operator-= (const BernsteinTetrahedralPoly& poly);
    /// Subtraction with a scalar
    BernsteinTetrahedralPoly& operator-= (double c);

    /// Division with a scalar
    BernsteinTetrahedralPoly& operator/= (double c);

    /// Read from an input stream
    void read(std::istream& is);
    /// Write to an output stream
    void write(std::ostream& os) const;

private:
    int deg_;
    std::vector<double> coefs_;
};


/// Multiplication of two polynomials
inline BernsteinTetrahedralPoly
operator* (const BernsteinTetrahedralPoly& p1,
	   const BernsteinTetrahedralPoly& p2)
{
    BernsteinTetrahedralPoly tmp = p1;
    tmp *= p2;
    return tmp;
}


/// Multiplication of a polynomial with a scalar
inline BernsteinTetrahedralPoly
operator* (const BernsteinTetrahedralPoly& p, double c)
{
    BernsteinTetrahedralPoly tmp = p;
    tmp *= c;
    return tmp;
}


/// Multiplication of a scalar with a polynomial
inline BernsteinTetrahedralPoly
operator* (double c, const BernsteinTetrahedralPoly& p)
{
    BernsteinTetrahedralPoly tmp = p;
    tmp *= c;
    return tmp;
}


/// Addition of two polynomials
inline BernsteinTetrahedralPoly
operator+ (const BernsteinTetrahedralPoly& p1,
	   const BernsteinTetrahedralPoly& p2)
{
    BernsteinTetrahedralPoly tmp = p1;
    tmp += p2;
    return tmp;
}


/// Addition of a scalar with a polynomial
inline BernsteinTetrahedralPoly
operator+ (double c, const BernsteinTetrahedralPoly& p)
{
    BernsteinTetrahedralPoly tmp = p;
    tmp += c;
    return tmp;
}


/// Addition of a polynomial with a scalar
inline BernsteinTetrahedralPoly
operator+ (const BernsteinTetrahedralPoly& p, double c)
{
    BernsteinTetrahedralPoly tmp = p;
    tmp += c;
    return tmp;
}


/// Subtraction of two polynomials
inline BernsteinTetrahedralPoly
operator- (const BernsteinTetrahedralPoly& p1,
	   const BernsteinTetrahedralPoly& p2)
{
    BernsteinTetrahedralPoly tmp = p1;
    tmp -= p2;
    return tmp;
}


/// Subtraction of a polynomial from a scalar
inline BernsteinTetrahedralPoly
operator- (double c, const BernsteinTetrahedralPoly& p)
{
    BernsteinTetrahedralPoly tmp(c);
    tmp -= p;
    return tmp;
}


/// Subtraction of a scalar from a polynomial
inline BernsteinTetrahedralPoly
operator- (const BernsteinTetrahedralPoly& p, double c)
{
    BernsteinTetrahedralPoly tmp = p;
    tmp -= c;
    return tmp;
}


/// Division of a polynomial with a scalar
inline BernsteinTetrahedralPoly
operator/ (const BernsteinTetrahedralPoly& p, double c)
{
    BernsteinTetrahedralPoly tmp = p;
    tmp /= c;
    return tmp;
}


/// Read BernsteinTetrahedralPoly from input stream
inline std::istream& operator >> (std::istream& is,
				  Go::BernsteinTetrahedralPoly& p)
{
    p.read(is);
    return is;
}


/// Write BernsteinTetrahedralPoly to output stream
inline std::ostream& operator << (std::ostream& os,
				  const Go::BernsteinTetrahedralPoly& p)
{
    p.write(os);
    return os;
}


} // namespace Go


#endif // _BERNSTEINTETRAHEDRALPOLY_H

