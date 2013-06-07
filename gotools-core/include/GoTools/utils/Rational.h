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

