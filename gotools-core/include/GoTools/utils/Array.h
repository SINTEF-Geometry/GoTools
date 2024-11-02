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

#ifndef _ARRAY_H
#define _ARRAY_H


#include "GoTools/utils/errormacros.h"
#include "GoTools/utils/config.h"
#include <iostream>
#include <cmath>
#include <algorithm>


namespace Go
{

    /** Compile-time sized array.
     *  Encapsulates an array with full value semantics, i.e. copying
     *  and assignment works as expected. The class also has some
     *  vector algebra functionality, such as scalar product,
     *  multiplication by scalars etc. In such contexts the objects
     *  are talked about as 'vectors' and not 'Arrays'.
     */
    template <typename T, int Dim>
    class Array
    {
    private:
	T p_[Dim];

    public:
	/// Default constructor, does not initialize elements.
	Array()
	{}
      
	/// Constructor filling the 'Dim'-tuple with the same value. (J.O. 080211)
	/// I don't think this one should interfere with other constructors, particularly since there
	/// doesn't seem to be support for an M==1 case already...?!
	explicit Array(T x)
	{
	    for (int i=0; i<Dim; i++)
		p_[i] = x;
	}

	/// Constructor taking 2 arguments, only compiles if
	/// Array dimension argument is 2.
	Array(T x, T y)
	{
          static_assert(Dim == 2, "Expecting Dim == 2");
	    p_[0] = x;
	    p_[1] = y;
	}
	/// Constructor taking 3 arguments, only compiles if
	/// Array dimension argument is 3.
	Array(T x, T y, T z)
	{
          static_assert(Dim == 3, "Expecting Dim == 3");
	    p_[0] = x;
	    p_[1] = y;
	    p_[2] = z;
	}
	/// Constructor taking 4 arguments, only compiles if
	/// Array dimension argument is 4.
	Array(T x, T y, T z, T w)
	{
	    static_assert(Dim == 4, "Expecting Dim == 4");
	    p_[0] = x;
	    p_[1] = y;
	    p_[2] = z;
	    p_[3] = w;
	}
	/// Copy constructor.
	Array(const Array& v)
	{
	    setValue(v.p_);
	}
	/// Constructor taking a pointer or other
	/// random access iterator and copying all
	/// elements into the Array.
	template <typename RandomAccessIterator>
	explicit Array(RandomAccessIterator start)
	{
	    setValue(start);
	}
	/// Takes a pointer or other
	/// random access iterator and copies all
	/// elements into the Array.
	template <typename RandomAccessIterator>
	void setValue(RandomAccessIterator from)
	{
	    for (int i = 0; i < Dim; ++i)
		p_[i] = from[i];
	}
	/// Constructor that takes an Array of a
	/// compatible type.
	template <typename U>
	explicit Array(const Array<U, Dim>& v)
	{
	    setValueConvert(v.begin());
	}
	/// Takes a pointer or other
	/// random access iterator with compatible
	/// value type and copies all
	/// elements into the Array.
	template <typename RandomAccessIterator>
	void setValueConvert(RandomAccessIterator from)
	{
	    for (int i = 0; i < Dim; ++i)
		p_[i] = T(from[i]);
	}
	/// Assignment operator.
	Array& operator = (const Array& v)
	{
	    if (&v != this) {
		setValue(v.p_);
	    }
	    return *this;
	}

	/// Reads an Array elementwise from
	/// a standard istream.
	void read(std::istream& is)
	{
	    for (int i = 0; i < Dim; ++i)
		is >> p_[i];
	}
	/// Writes an Array elementwise to
	/// a standard ostream. Precision is
	/// set to 16 by this function.
	void write(std::ostream& os) const
	{
	    os.precision(16);
	    for (int i = 0; i < Dim-1; ++i)
		os << p_[i] << ' ';
	    os << p_[Dim-1];
	}

	/// Access to first element.  (Deprecated, kept for backward
	/// compatibility)
	const T& x() const { return p_[0]; }
	/// Access to first element.  (Deprecated, kept for backward
	/// compatibility)
	T& x()       { return p_[0]; }
	/// Access to first element.  (Deprecated, kept for backward
	/// compatibility)
	const T& y() const { return p_[1]; }
	/// Access to first element.  (Deprecated, kept for backward
	/// compatibility)
	T& y()       { return p_[1]; }
	/// Access to first element.  (Deprecated, kept for backward
	/// compatibility)
	const T& z() const { return p_[2]; }
	/// Access to first element.  (Deprecated, kept for backward
	/// compatibility)
	T& z()       { return p_[2]; }

	/// Read-only index access.
	const T& operator [] (int i) const { return p_[i]; }
	/// Index access.
	T& operator [] (int i)       { return p_[i]; }

	/// Get a read-only start iterator.
	const T* begin() const { return p_; }
	/// Get a start iterator.
	T* begin()       { return p_; }
	/// Get a read-only end iterator.
	const T* end() const { return p_+Dim; }
	/// Get a end iterator.
	T* end()       { return p_+Dim; }

	/// Get the dimension of the Array.
	int size() const { return Dim; }

	/// Get the square of the euclidian length
	/// of the vector.
	T length2() const
	{
	    T l2(0.0);
	    for (int i = 0; i < Dim; ++i)
		l2 += p_[i]*p_[i];
	    return l2;
	}
	/// Get the euclidian length
	/// of the vector.
	T length() const
	{
#ifdef __BORLANDC__
	    return std::sqrt(length2());
#else
	    return sqrt(length2());
#endif
	}
	/// Get the infinity-norm (or max-norm) length
	/// of the vector.
	T lengthInf() const
	{
	    T linf(0.0);
	    for (int i = 0; i < Dim; ++i) {
		linf = std::max(linf, std::abs(p_[i]));
	    }
	    return linf;
	}
	/// Get the square of the euclidian length
	/// of the difference between this vector
	/// and another vector.
	T dist2(const Array &v) const
	{
	    T l2(0.0);
	    T d;
	    for (int i = 0; i < Dim; ++i) {
		d = p_[i] - v.p_[i];
		l2 += d*d;
	    }
	    return l2;
	}
	/// Get the euclidian length
	/// of the difference between this vector
	/// and another vector.
	T dist(const Array &v) const
	{
#ifdef __BORLANDC__
	    return std::sqrt(dist2(v));
#else
	    return sqrt(dist2(v));
#endif
	}
	/// Get the infinity-norm (or max-norm) length
	/// of the difference between this vector
	/// and another vector.
	T distInf(const Array &v) const
	{
	    T linf(0.0);
	    T d;
	    for (int i = 0; i < Dim; ++i) {
		d = p_[i] - v.p_[i];
		linf = std::max(linf, std::abs(d));
	    }
	    return linf;
	}

	/// Normalize this vector, i.e. divide every element
	/// by length().
	void normalize()
	{
	    (*this) /= length();
	}

	void normalize_checked()
	{
	  double len = length();
	  if (len > 1.0e-12)
	    (*this) /= len;
	}

	/// The sum of two vectors.
	Array operator + (const Array &v) const
	{
	    Array res(*this);
	    res += v;
	    return res;
	}

	/// Add a vector to this vector.
	bool operator == (const Array &v) const
	{
	  for (int i = 0; i < Dim; ++i)
	    if (p_[i] != v.p_[i])
	      return false;
	  return true;
	}
	/// Add a vector to this vector.
	Array& operator += (const Array &v)
	{
	    for (int i = 0; i < Dim; ++i)
		p_[i] += v.p_[i];
	    return *this;
	}

	/// The difference between two vectors.
	Array operator - (const Array &v) const
	{
	    Array res(*this);
	    res -= v;
	    return res;
	}

	/// Subtract a vector from this vector.
	Array& operator -=(const Array &v)
	{
	    for (int i = 0; i < Dim; ++i)
		p_[i] -= v.p_[i];
	    return *this;
	}

	
	/// The product of a vector and a scalar.
	Array operator * (T d) const
	{
	    Array res(*this);
	    res *= d;
	    return res;
	}

	/// Multiply this vector by a scalar.
	Array& operator *= (T d)
	{
	    for (int i = 0; i < Dim; ++i)
		p_[i] *= d;
	    return *this;
	}

	/// A vector divided by a scalar.
	Array  operator / (double d) const
	{
	    Array res(*this);
	    res /= d;
	    return res;
	}
	/// Divide this vector with a scalar.
	Array& operator /= (double d)
	{
	    for (int i = 0; i < Dim; ++i)
		p_[i] /= d;
	    return *this;
	}

	/// The negation of a vector.
	Array operator - () const
	{
	    Array res(*this);
	    for (int i = 0; i < Dim; ++i)
		res.p_[i] = - p_[i];
	    return res;
	}

	/// The scalar product (or inner product, or dot
	/// product) of two vectors.
	T operator * (const Array &v) const
	{
	    T res(0.0);
	    for (int i = 0; i < Dim; ++i)
		res += p_[i]*v.p_[i];
	    return res;
	}

	/// The cross product of two vectors.
	/// Only compiles if dimension is 3.
	Array operator % (const Array &v) const
	{
	    static_assert(Dim == 3, "Expecting Dim == 3");
#ifdef __BORLANDC__
	    return Array<T, Dim>(p_[1]*v.p_[2] - p_[2]*v.p_[1],
			 p_[2]*v.p_[0] - p_[0]*v.p_[2],
			 p_[0]*v.p_[1] - p_[1]*v.p_[0]);
#else
	    return Array(p_[1]*v.p_[2] - p_[2]*v.p_[1],
			 p_[2]*v.p_[0] - p_[0]*v.p_[2],
			 p_[0]*v.p_[1] - p_[1]*v.p_[0]);
#endif
	}

	/// The cross product of two vectors.
	/// Only compiles if dimension is 3.
	Array cross(const Array &v) const
	{
	    return operator%(v);
	}

	/// The cosine of the angle between this and
	/// another vector.
	T cosAngle(const Array& v) const
	{
	    T cosang = ((*this)*v)/(length()*v.length());
	    // Roundoff may push us outside the [-1, 1] interval,
	    // which may lead to no end of trouble.
	    cosang = std::max(T(-1.0), cosang);
	    cosang = std::min(T(1.0), cosang);
	    return cosang;
	}

	/// The angle between this and
	/// another vector.
	T angle(const Array& v) const
	{
#ifdef __BORLANDC__
	    return std::acos(cosAngle(v));
#else
	    return acos(cosAngle(v));
#endif
	}

	/// Sets the vector to the zero vector.
	void zero()
	{
	    for (int i = 0; i < Dim; ++i) {
		p_[i] = T(0.0);
	    }
	}

    };

    /// Typedef for ease of use in frequently used case.
    typedef Array<double, 2> Vector2D;
    /// Typedef for ease of use in frequently used case.
    typedef Array<double, 3> Vector3D;
    /// Typedef for ease of use in frequently used case.
    typedef Array<double, 4> Vector4D;


    /// The product of a vector and a scalar.
    template<typename T, int Dim>
    inline Go::Array<T, Dim> operator * (T d, const Go::Array<T, Dim>& v)
    { return v*d; }

    /// Stream extraction for Array.
    template <typename T, int Dim>
    inline std::istream& operator >> (std::istream& is, Go::Array<T,Dim>& v)
    { v.read(is); return is; }

    /// Stream insertion for Array.
    template <typename T, int Dim>
    inline std::ostream& operator << (std::ostream& os, const Go::Array<T,Dim>& v)
    { v.write(os); return os; }


} // namespace Go


#endif // _ARRAY_H



