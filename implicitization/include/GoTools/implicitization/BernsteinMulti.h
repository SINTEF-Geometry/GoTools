//===========================================================================
//                                                                           
// File: BernsteinMulti.h                                                  
//                                                                           
// Created: Wed May 16 11:22:13 2001                                         
//                                                                           
// Author: Jan B. Thomassen <jbt@math.sintef.no>
//                                                                           
// Revision: $Id: BernsteinMulti.h,v 1.45 2006-03-29 12:12:36 jbt Exp $
//                                                                           
// Description:
//                                                                           
//===========================================================================

#ifndef _BERNSTEINMULTI_H
#define _BERNSTEINMULTI_H

#include "GoTools/utils/Array.h"
#include <vector>
#include <iostream>


namespace Go {


class BernsteinPoly;


/**
 * Class that implements bivariate tensor product Bernstein
 * polynomials on the domain [0,1]x[0,1].
 *
 * The formula for a degree (m,n) polynomial is
 * \f[
 *   p(u,v) = \sum_{i=0}^m \sum_{j=0}^n B_i(u) B_j(v) c_{ij}
 * \f]
 * where the basis functions in the first parameter direction are
 * defined by 
 * \f[
 *   B_i(u) = {m \choose i} u^i(1-u)^{n-i}
 * \f]
 * The coefficients \f$c_{ij}\f$ are stored in the following order:
 * \f$ c_{00}, c_{10}, c_{20}, \ldots c_{mn} \f$.
 */

class BernsteinMulti {
public:
    /// Default constructor
    BernsteinMulti() : degu_(-1), degv_(-1) { }
    /// Constructor.
    /// \param m degree in the first parameter direction
    /// \param n degree in the second parameter direction
    /// \param coefs vector of Bernstein coefficients
    BernsteinMulti(int m, int n,
		   const std::vector<double>& coefs)
	: degu_(m), degv_(n), coefs_(coefs)
    { 
	ALWAYS_ERROR_IF((degu_+1) * (degv_+1) != (int)coefs_.size(),
			"Number of coefficients does not match degrees");
    }
    /// Constructor
    /// \param m degree in the first parameter direction
    /// \param n degree in the second parameter direction
    /// \param begin iterator to start of container with Bernstein
    /// coefficients
    /// \param end iterator to end of container
    template <typename ForwardIterator>
    BernsteinMulti(int m, int n,
		   ForwardIterator begin, ForwardIterator end)
	: degu_(m), degv_(n), coefs_(begin, end) { }
    /// Constructs a constant polynomial.
    /// \param coef the constant value of the polynomial
    explicit BernsteinMulti(double coef)
	: degu_(0), degv_(0), coefs_(std::vector<double>(1, coef)) { }

    /// Get the degree in the first parameter direction.
    /// \return the degree in the first parameter direction
    int degreeU() const
    { return degu_; }
    /// Get the degree in the second parameter direction.
    /// \return the degree in the second parameter direction
    int degreeV() const
    { return degv_; }

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
    /// \param u first parameter
    /// \param v second parameter
    /// \return the value of the polynomial at (u,v)
    double operator() (double u, double v) const;
    /// Evaluation operator - as above, but taking the arguments in a
    /// Vector2D.
    /// \param pt parameter point
    /// \return the value of the polynomial at pt
    double operator() (const Vector2D& pt) const
    { return (*this)(pt[0], pt[1]); }
    /// Evaluation operator - as above, but taking the arguments in a
    /// std::vector.
    /// \param pt parameter point
    /// \return the value of the polynomial at pt
    double operator() (const std::vector<double>& pt) const
    { return (*this)(pt[0], pt[1]); }

    /// Check if the polynomial is the zero function
    /// \param eps the threshold value with which every coefficient
    /// is compared
    /// \return true if for every coefficient c, \f$|c| <= eps\f$,
    /// false otherwise
    bool isZero(const double eps = 0.0) const;
    /// Check if polynomial is strictly positive
    /// \param eps the threshold value with which every coefficient
    /// is compared
    /// \return true if for every coefficient c, \f$c > eps\f$, false
    /// otherwise
    bool isStrictlyPositive(const double eps = 0.0) const;
    /// Check if polynomial is strictly negative
    /// \param eps the threshold value with which every coefficient
    /// is compared
    /// \return true if for every coefficient c, \f$c < -eps\f$, false
    /// otherwise
    bool isStrictlyNegative(const double eps = 0.0) const;
    /// Check if polynomial is non-negative
    /// \param eps the threshold value with which every coefficient
    /// is compared
    /// \return true if for every coefficient c, \f$c >= -eps\f$,
    /// false otherwise
    bool isNonNegative(const double eps = 0.0) const;

    /// Calculates the norm of the polynomial, defined as the sum of
    /// absolute values of the coefficients divided by the number of
    /// coefficients. This norm is not the same as the L1 norm of the
    /// polynomial on the domain [0,1]x[0,1], unless the polynomial
    /// does not change sign. In any case this norm is greater than or
    /// equal to the L1 norm of the polynomial.
    double norm() const;
    /// Normalizes the polynomial by dividing all coefficients with
    /// norm(), which is close to (but not identical to) the L1 norm
    /// of the polynomial. \see norm()
    void normalize();

    /// Calculates the mean value of the coefficients.
    double mean() const;

    /// Finds the derivative of the polynomial in terms of a new
    /// BernsteinMulti.
    /// \param der1 the order of derivatives requested in the first
    /// parameter direction
    /// \param der2 the order of derivatives requested in the second
    /// parameter direction
    /// \returns a new polynomial q given by 
    /// \f[ q(u,v) = D^{(der1,der2)}p(u,v). \f]
    BernsteinMulti deriv(int der1, int der2) const;
    /// Calculates the determinant of the Hessian as a BernsteinMulti
    /// Describes Gaussian curvature of the graph of the polynomial
    BernsteinMulti detHess() const;
    /// Calculates the trace of the Hessian as a BernsteinMUlti
    /// Describes mean curvature (modulo a factor 2) of the graph of
    /// the polynomial
    BernsteinMulti traceHess() const;

    /// Evaluates the blossom of the polynomial. The blossom of the
    /// polynomial \f$p\f$ of degree (m,n) is a multi-affine
    /// (m+n)-variate function \f$B_p\f$ that satisfies
    /// \f[
    ///   B_p(\underbrace{u, \ldots, u}_m, \underbrace{v, \ldots, v}_n)
    ///   = p(u,v)
    /// \f]
    /// and is symmetric in the first m and last n arguments
    /// separately.
    /// \param uvec the vector of values at which the blossom is
    /// evaluated in the first parameter direction
    /// \param vvec the vector of values at which the blossom is
    /// evaluated in the second parameter direction
    /// \returns the value of the blossom at (uvec,vvec)
    double blossom(const std::vector<double>& uvec,
		   const std::vector<double>& vvec) const;
    /// Finds a new BernsteinMulti representing the original on some
    /// rectangular domain. The new polynomial is defined on
    /// [0,1]x[0,1]. If p is the
    /// old polynomial and q is the new one, q is defined by
    /// \f[
    ///   q(u,v) = p((1-u) \cdot u_0 + u \cdot u_1, (1-v) \cdot v_0 + v
    ///    \cdot v_1).
    /// \f]
    /// The function is implemented by blossoming.
    /// \param u0 the start of the interval in the first parameter
    /// direction 
    /// \param u1 the end of the interval in the first parameter
    /// direction
    /// \param v0 the start of the interval in the second parameter
    /// direction 
    /// \param v1 the end of the interval in the second parameter
    /// direction
    /// \returns a polynomial q as specified above
    BernsteinMulti pickDomain(double u0, double u1,
			      double v0, double v1) const;
    /// Routine that picks out the line with endpoints a and
    /// b as a new BernsteinPoly. The new polynomial is defined
    /// on [0,1]. If p is the
    /// old bivariate polynomial and q is the new one, q is defined by
    /// \f[ q(t) = p((1-t)a_0 + tb_0, (1-t)a_1 + tv_1). \f]
    /// The function is implemented by blossoming.
    BernsteinPoly pickLine(Array<double,2> a, 
			   Array<double,2> b) const;

    /// Degree elevation.
    /// \param du the polynomial degree in which you want to represent
    /// this polynomial in the first parameter direction
    /// \param dv the polynomial degree in which you want to represent
    /// this polynomial in the second parameter direction
    void degreeElevate(int du, int dv);

    /// Multiplication with another polynomial
    BernsteinMulti& operator*= (const BernsteinMulti& multi);
    /// Multiplication with a scalar
    BernsteinMulti& operator*= (double c);
    
    /// Addition with another polynomial
    BernsteinMulti& operator+= (const BernsteinMulti& multi);
    /// Addition with a scalar
    BernsteinMulti& operator+= (double c);
    
    /// Subtraction with another polynomial
    BernsteinMulti& operator-= (const BernsteinMulti& multi);
    /// Subtraction with a scalar
    BernsteinMulti& operator-= (double c);
    
    /// Division with a scalar
    BernsteinMulti& operator/= (double c);

    /// Function that "binds" the u-parameter.
    /// \return a BernsteinPoly
    /// in the v-direction for a constant value of u.
    BernsteinPoly bindU(double u) const;
    /// Function that "binds" the v-parameter.
    /// \return a BernsteinPoly
    /// in the u-direction for a constant value of v.
    BernsteinPoly bindV(double v) const;

    /// Read from an input stream
    void read(std::istream& is);
    /// Write to an output stream
    void write(std::ostream& os) const;

private:
    // Degrees in the u'th and v'th direction.
    int degu_, degv_;
    // Coefficients are stored in a vector with (degu_+1)(degv_+1)
    // entries. The "u-index" runs fastest; i.e. for control points
    // c_ij, the sequence is
    // c00, c10, c20, ..., c01, c11, c21, ..., cmn
    std::vector<double> coefs_;
};


/// Multiplication of two polynomials
inline BernsteinMulti operator* (const BernsteinMulti& m1,
				 const BernsteinMulti& m2)
{
    BernsteinMulti tmp = m1;
    tmp *= m2;
    return tmp;
}


/// Multiplication of a polynomial with a scalar
inline BernsteinMulti operator* (const BernsteinMulti& m1,
				 double c)
{
    BernsteinMulti tmp = m1;
    tmp *= c;
    return tmp;
}


/// Multiplication of a scalar with a polynomial
inline BernsteinMulti operator* (double c,
				 const BernsteinMulti& m1)
{
    BernsteinMulti tmp = m1;
    tmp *= c;
    return tmp;
}


/// Addition of two polynomials
inline BernsteinMulti operator+ (const BernsteinMulti& m1,
				 const BernsteinMulti& m2)
{
    BernsteinMulti tmp = m1;
    tmp += m2;
    return tmp;
}


/// Addition of a polynomial with a scalar
inline BernsteinMulti operator+ (const BernsteinMulti& m,
				 double c)
{
    BernsteinMulti tmp = m;
    tmp += c;
    return tmp;
}


/// Addition of a scalar with a polynomial
inline BernsteinMulti operator+ (double c,
				 const BernsteinMulti& m)
{
    BernsteinMulti tmp = m;
    tmp += c;
    return tmp;
}


/// Subtraction of two polynomials
inline BernsteinMulti operator- (const BernsteinMulti& m1,
				 const BernsteinMulti& m2)
{
    BernsteinMulti tmp = m1;
    tmp -= m2;
    return tmp;
}


/// Subtraction of a scalar from a polynomial
inline BernsteinMulti operator- (const BernsteinMulti& m,
				 double c)
{
    BernsteinMulti tmp = m;
    tmp -= c;
    return tmp;
}


/// Subtraction of a polynomial from a scalar
inline BernsteinMulti operator- (double c,
				 const BernsteinMulti& m)
{
    BernsteinMulti tmp(c);
    tmp -= m;
    return tmp;
}


/// Division of a polynomial with a scalar
inline BernsteinMulti operator/ (const BernsteinMulti& m,
				 double c)
{
    BernsteinMulti tmp = m;
    tmp /= c;
    return tmp;
}


} // namespace Go


namespace std {


/// Read BernsteinMulti from input stream
inline std::istream& operator >> (std::istream& is,
				  Go::BernsteinMulti& m)
{
    m.read(is);
    return is;
}


/// Write BernsteinMulti to output stream
inline std::ostream& operator << (std::ostream& os,
				  const Go::BernsteinMulti& m)
{
    m.write(os);
    return os;
}


} // namespace std


#endif // _BERNSTEINMULTI_H

