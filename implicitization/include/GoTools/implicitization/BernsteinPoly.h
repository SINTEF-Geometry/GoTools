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

#ifndef _BERNSTEINPOLY_H
#define _BERNSTEINPOLY_H

#include <vector>
#include <iostream>

namespace Go {


/**
 * Class that implements Bernstein polynomials on the interval
 * [0,1]. The formula for Bernstein polynomials is
 * \f[
 *    p(t) = \sum_{i=0}^n {n \choose i} t^i(1-t)^{n-i}c_i.
 * \f] 
 * The coefficients \f$c_i\f$ are stored in a vector of doubles.
 */

class BernsteinPoly {
public:
    /// Default constructor
    BernsteinPoly() { }
    /// Constructor
    /// \param coefs a vector of doubles containing the Bernstein
    /// coefficients
    explicit BernsteinPoly(const std::vector<double>& coefs)
	: coefs_(coefs) { }
    /// Constructor
    /// \param begin iterator to start of container containing the
    /// Bernstein coefficients
    /// \param end iterator to end of same container
    template <typename ForwardIterator>
    BernsteinPoly(ForwardIterator begin, ForwardIterator end)
	: coefs_(begin, end) { }
    /// Constructor for a constant polynomial
    /// \param coef the constant value of the polynomial
    explicit BernsteinPoly(double coef)
	: coefs_(std::vector<double>(1, coef)) { }

    /// Get the degree of the polynomial
    /// \return the degree of the polynomial
    int degree() const
    { return (int)coefs_.size() - 1; }

    /// Access function for coefficients
    const double& operator[] (int num) const
    { return coefs_[num]; }
    /// Access function for coefficients
    double& operator[] (int num)
    { return coefs_[num]; }
    /// Access function for coefficients
    std::vector<double>::iterator coefsBegin()
    { return coefs_.begin(); }
    /// Access function for coefficients
    std::vector<double>::const_iterator coefsBegin() const
    { return coefs_.begin(); }
    /// Access function for coefficients
    std::vector<double>::iterator coefsEnd()
    { return coefs_.end(); }
    /// Access function for coefficients
    std::vector<double>::const_iterator coefsEnd() const
    { return coefs_.end(); }

    /// Evaluation operator, based on the de Casteljau algorithm.
    /// \param t the parameter value at which the polynomial is
    /// evaluated
    /// \return the value of the polynomial at t
    double operator() (double t) const;

    /// Check if the polynomial is the zero function
    /// \param eps the threshold value with which every coefficient
    /// is compared
    /// \return true if for every coefficient c, \f$|c| <= eps\f$,
    /// false otherwise
    bool isZero(const double eps = 0.0) const;
    /// Check if the polynomial is strictly positive. NOTE: If the
    /// return value is 'true', then the polynomial is guaranteed to
    /// be strictly positive. If the return value is 'false', the
    /// situation is undetermined and the polynomial may still be
    /// strictly positive.
    /// \param eps the threshold value with which every coefficient
    /// is compared
    /// \return true if for every coefficient c, \f$c > eps\f$, false
    /// otherwise
    bool isStrictlyPositive(const double eps = 0.0) const;

    /// Calculates a norm, defined as the sum of the absolute values
    /// of the coefficients divided by the number of
    /// coefficients. This norm is not the same as the L1 norm of the
    /// polynomial (on the interval [0,1]), unless the 
    /// polynomial does not change sign. In any case this norm is
    /// greater than or equal to the L1 norm of the polynomial.
    double norm() const;
    /// Normalizes the polynomial by dividing all coefficients with
    /// norm(), which is close to (but not identical to) the L1 norm
    /// of the polynomial. \see norm()
    void normalize();

    /// Calculates the derivative in terms of a new BernsteinPoly.
    /// \param der the order of derivatives requested
    /// \returns a new polynomial q given by 
    /// \f[ q(t) = p^{(der)}(t). \f]
    BernsteinPoly deriv(int der) const;

    /// Calculates the integral of the polynomial over an
    /// interval. The default interval is [0,1].
    double integral(double a = 0.0, double b = 1.0) const;

    /// Evaluates the blossom of the polynomial. The blossom of the
    /// polynomial \f$p\f$ of degree n is a multi-affine n-variate
    /// symmetric
    /// function \f$B_p\f$ that satisfies \f$B_p(t,\ldots,t) =
    /// p(t)\f$.
    /// \param tvec the vector of values at which the blossom is
    /// evaluated.
    /// \returns the value of the blossom at tvec
    double blossom(const std::vector<double>& tvec) const;
    /// Finds a new BernsteinPoly representing the original on some
    /// interval. The new polynomial is defined on [0,1]. If p is the
    /// old polynomial and q is the new one, q is defined by
    /// \f[ q(t) = p(a(1-t) + tb). \f]
    /// \param a the start of the interval
    /// \param b the end of the interval
    /// \returns a polynomial q as specified above
    BernsteinPoly pickInterval(double a, double b) const;

    /// Degree elevation.
    /// \param d the polynomial degree in which you want to represent
    /// this polynomial.
    void degreeElevate(int d);

    /// Multiplication with another polynomial
    BernsteinPoly& operator*= (const BernsteinPoly& poly);
    /// Multiplication with a scalar
    BernsteinPoly& operator*= (double c);

    /// Addition with another polynomial
    BernsteinPoly& operator+= (const BernsteinPoly& poly);
    /// Addition with a scalar
    BernsteinPoly& operator+= (double c);

    /// Subtraction with another polynomial
    BernsteinPoly& operator-= (const BernsteinPoly& poly);
    /// Subtraction with a scalar
    BernsteinPoly& operator-= (double c);

    /// Division with a scalar
    BernsteinPoly& operator/= (double c);

    /// Read from input stream
    void read(std::istream& is);
    /// Write to output stream
    void write(std::ostream& os) const;

private:
    std::vector<double> coefs_;
};


/// Multiplication of two polynomials
inline BernsteinPoly operator* (const BernsteinPoly& p1,
				const BernsteinPoly& p2)
{
    BernsteinPoly tmp = p1;
    tmp *= p2;
    return tmp;
}


/// Multiplication of a polynomial with a scalar
inline BernsteinPoly operator* (const BernsteinPoly& p,
				double c)
{
    BernsteinPoly tmp = p;
    tmp *= c;
    return tmp;
}


/// Multiplication of a scalar with a polynomial
inline BernsteinPoly operator* (double c,
				const BernsteinPoly& p)
{
    BernsteinPoly tmp = p;
    tmp *= c;
    return tmp;
}


/// Addition of two polynomials
inline BernsteinPoly operator+ (const BernsteinPoly& p1,
				const BernsteinPoly& p2)
{
    BernsteinPoly tmp = p1;
    tmp += p2;
    return tmp;
}


/// Addition of a scalar with a polynomial
inline BernsteinPoly operator+ (double c,
				const BernsteinPoly& p)
{
    BernsteinPoly tmp = p;
    tmp += c;
    return tmp;
}


/// Addition of a polynomial with a scalar
inline BernsteinPoly operator+ (const BernsteinPoly& p,
				double c)
{
    BernsteinPoly tmp = p;
    tmp += c;
    return tmp;
}


/// Subtraction of two polynomials
inline BernsteinPoly operator- (const BernsteinPoly& p1,
				const BernsteinPoly& p2)
{
    BernsteinPoly tmp = p1;
    tmp -= p2;
    return tmp;
}


/// Subtraction of a polynomial from a scalar
inline BernsteinPoly operator- (double c,
				const BernsteinPoly& p)
{
    BernsteinPoly tmp(c);
    tmp -= p;
    return tmp;
}


/// Subtraction of a scalar from a polynomial
inline BernsteinPoly operator- (const BernsteinPoly& p,
				double c)
{
    BernsteinPoly tmp = p;
    tmp -= c;
    return tmp;
}


/// Division of a polynomial with a scalar
inline BernsteinPoly operator/ (const BernsteinPoly& p,
				double c)
{
    BernsteinPoly tmp = p;
    tmp /= c;
    return tmp;
}


} // namespace Go


namespace std {


/// Read BernsteinPoly from input stream
inline std::istream& operator >> (std::istream& is,
				  Go::BernsteinPoly& p)
{
    p.read(is);
    return is;
}


/// Write BernsteinPoly to output stream
inline std::ostream& operator << (std::ostream& os,
				  const Go::BernsteinPoly& p)
{
    p.write(os);
    return os;
}


} // namespace std


#endif // _BERNSTEINPOLY_H

