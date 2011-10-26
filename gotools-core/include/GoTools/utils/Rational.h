//===========================================================================
//                                                                           
// File: Rational.h                                                          
//                                                                           
// Created: Tue Dec 10 17:11:16 2002                                         
//                                                                           
// Author: Atgeirr F Rasmussen <atgeirr@sintef.no>
//                                                                           
// Revision: $Id: Rational.h,v 1.3 2005-06-09 07:29:40 oan Exp $
//                                                                           
// Description:
//                                                                           
//===========================================================================

#ifndef _RATIONAL_H
#define _RATIONAL_H

#include <iostream>
#include <cstdlib>

namespace Go {
    
    /** Class representing rational numbers
     *
     */


class Rational
{
public:
    
    /// Construct a rational number whose value is 0.
    Rational() : p_(0), q_(1) {}

    /// Construct a rational number whose value is 'p' (integer)
    Rational(int p) : p_(p), q_(1) {}

    /// Construct a rational number whose value is 'p'/'q' (fraction)
    Rational(int p, int q) : p_(p), q_(q) {}

    /// Add a different rational number to this rational number
    Rational& operator += (const Rational& other)
    {
	p_ = p_*other.q_ + q_*other.p_;
	q_ = q_*other.q_;
	simplify();
	return *this;
    }

    /// Subtract a different rational number from this rational number
    Rational& operator -= (const Rational& other)
    {
	Rational tmp = -other;
	(*this) += tmp;
	return *this;
    }

    /// Multiply this rational number with a different rational number
    Rational& operator *= (const Rational& other)
    {
	p_ = p_*other.p_;
	q_ = q_*other.q_;
	simplify();
	return *this;
    }

    /// Divide this rational number by a different rational number
    Rational& operator /= (const Rational& other)
    {
	p_ = p_*other.q_;
	q_ = q_*other.p_;
	simplify();
	return *this;
    }

    /// Return the additive inverse of this rational number
    Rational operator- () const
    {
	return Rational(-p_, q_);
    }

    /// Test this rational number for equality with another rational number
    bool operator == (const Rational r)
    {
	return (p_ == r.p_ && q_ == r.q_);
    }

    /// Test this rational number for difference with another rational number
    bool operator != (const Rational r)
    {
	return (p_ != r.p_ || q_ != r.q_);
    }

    /// Write this rational number to a stream
    void write(std::ostream& os) const
    {
	os << p_ << '/' << q_;
    }

    /// Simplify the internal fractional expression of this rational number
    void simplify()
    {
	int n = std::min(abs(p_), abs(q_));
	for (int i = 2; i <= n; ++i) {
	    while (p_%i==0 && q_%i==0) {
		p_ /= i;
		q_ /= i;
		n /= i;
	    }
	}
    }

private:
    int p_;
    int q_;
};

Rational operator + (const Rational& r1, const Rational r2)
{
    Rational res = r1;
    res += r2;
    return res;
}

Rational operator - (const Rational& r1, const Rational r2)
{
    Rational res = r1;
    res -= r2;
    return res;
}

Rational operator * (const Rational& r1, const Rational r2)
{
    Rational res = r1;
    res *= r2;
    return res;
}

Rational operator / (const Rational& r1, const Rational r2)
{
    Rational res = r1;
    res /= r2;
    return res;
}

std::ostream& operator << (std::ostream& os, const Rational& p)
{
    p.write(os);
    return os;
}

}; // end namespace Go
#endif // _RATIONAL_H

