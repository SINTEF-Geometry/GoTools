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

#ifndef _BEZIERTRIANGLE_H
#define _BEZIERTRIANGLE_H


#include "GoTools/utils/Array.h"
#include "GoTools/implicitization/BernsteinTriangularPoly.h"
#include "newmat.h"
#include "GoTools/implicitization/Binomial.h"
#include <vector>

/// \namespace std 
/// We implement two extensions: \b identity_element and \b power.

namespace std
{


// identity_element (not part of the C++ standard).

template <class _Tp> inline _Tp identity_element(plus<_Tp>) {
  return _Tp(0);
}
template <class _Tp> inline _Tp identity_element(multiplies<_Tp>) {
  return _Tp(1);
}



// Returns __x ** __n, where __n >= 0.  _Note that "multiplication"
// is required to be associative, but not necessarily commutative.

template <class _Tp, class _Integer, class _MonoidOperation>
_Tp __power(_Tp __x, _Integer __n, _MonoidOperation __oper)
{
  if (__n == 0)
    return identity_element(__oper);
  else {
    while ((__n & 1) == 0) {
      __n >>= 1;
      __x = __oper(__x, __x);
    }

    _Tp __result = __x;
    __n >>= 1;
    while (__n != 0) {
      __x = __oper(__x, __x);
      if ((__n & 1) != 0)
        __result = __oper(__result, __x);
      __n >>= 1;
    }
    return __result;
  }
}


template <class _Tp, class _Integer>
inline _Tp __power(_Tp __x, _Integer __n)
{
  return __power(__x, __n, multiplies<_Tp>());
}

// Alias for the internal name __power.  Note that power is an extension,
// not part of the C++ standard.

template <class _Tp, class _Integer, class _MonoidOperation>
inline _Tp power(_Tp __x, _Integer __n, _MonoidOperation __oper)
{
  return __power(__x, __n, __oper);
}

template <class _Tp, class _Integer>
inline _Tp power(_Tp __x, _Integer __n)
{
  return __power(__x, __n);
}

} // namespace std

namespace Go
{


/// Not documented
template <int N>
class BezierTriangle
{
public:
    /// Default constructor
    BezierTriangle() : deg_(0), last_derivdir_(-1.0, -1.0, -1.0) { }

    /// Constructor
    BezierTriangle(int deg, const Array<double, N>* cp)
	: deg_(deg), last_derivdir_(-1.0, -1.0, -1.0)
    {
	int sz = (deg_+1)*(deg_+2)/2;
	std::vector<double> tmp(sz);
	for (int i = 0; i < N; ++i) {
	    for (int j = 0; j < sz; ++j) {
		tmp[j] = cp[j][i];
	    }
	    polys_[i] = BernsteinTriangularPoly(deg_, tmp);
	}
    }

    /// Not documented
    Array<double, N> eval(Vector3D beta)
    {
	Array<double, N> res;
	for (int i = 0; i < N; ++i) {
	    res[i] = polys_[i](beta);
	}
	return res;
    }

    /// Not documented
    Array<double, N> evalderiv(Vector3D beta, Vector3D dir)
    {
	if (last_derivdir_.dist2(dir) > 0.0) {
	    for (int i = 0; i < N; ++i) {
		polys_[i].deriv(1, dir, derivpolys_[i]);
	    }
	}
	last_derivdir_ = dir;
	Array<double, N> res;
	for (int i = 0; i < N; ++i) {
	    res[i] = derivpolys_[i](beta);
	}
	return res;
    }

    /// Not documented
    void interpolate(int deg,
		     const Vector3D* nodes,
		     const Array<double, N>* values)
    {
	deg_ = deg;
	// Setting up interpolation matrix
	int sz = (deg+1)*(deg+2)/2;
	Matrix A(sz, sz);
	Binomial bin;
	for (int row = 0; row < sz; ++row) {
	    Vector3D b = nodes[row];
	    int col = 0;
	    for (int i = deg; i >= 0; --i) {
		for (int j = deg - i; j >= 0; --j) {
		    int k = deg - i - j;
		    double rr = std::power(b[0], i)
			* std::power(b[1], j)
			* std::power(b[2], k);
		    rr *= bin.trinomial(deg, i, j);
		    A.element(row, col) = rr;
		    ++col;
		}
	    }
	}

	// Solving the system for N right hand sides
	CroutMatrix ALUfact = A;
	DEBUG_ERROR_IF(ALUfact.IsSingular(),
		 "Matrix is singular! This should never happen!");
	ColumnVector coef(sz);
	ColumnVector rhs(sz);
	for (int i = 0; i < N; ++i) {
	    for (int j = 0; j < sz; ++j) {
		rhs.element(j) = values[j][i];
	    }
	    coef = ALUfact.i()*rhs;
	    polys_[i] = BernsteinTriangularPoly(deg,
						coef.Store(),
						coef.Store() + sz);
	}
    }
   
private:
    int deg_;
    Array<BernsteinTriangularPoly, N> polys_;
    mutable Array<BernsteinTriangularPoly, N> derivpolys_;
    mutable Vector3D last_derivdir_;
};


} // namespace Go

#endif // _BEZIERTRIANGLE_H

