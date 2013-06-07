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

#ifndef _MATRIXXD_H
#define _MATRIXXD_H

#include "GoTools/utils/Array.h"
#include <algorithm>
#include <iostream>
#include <cmath>


namespace Go
{

    /** Square n-dimensional (compile time constant dim) matrix.
     *
     */
	
template <typename T, int Dim>
class MatrixXD
{
public:
    /// Default constructor
    MatrixXD() {}

    /// Default destructor
    ~MatrixXD() {}

    /// Element access operator
    T& operator()(int row, int col)
    {
	return m_[row][col];
    }

    /// Access data array directly 
    T** get() {return m_;}

    /// Const element access operator.
    const T& operator()(int row, int col) const
    {
	return m_[row][col];
    }

    /// Set all elements to zero
    void zero()
    {
	for (int i = 0; i < Dim; ++i) {
	    for (int j = 0; j < Dim; ++j) {
		m_[i][j] = 0;
	    }
	}
    }

    /// Identity matrix
    void identity()
    {
	zero();
	for (int i = 0; i < Dim; ++i) {
	    m_[i][i] = 1;
	}
    }

    /// Transpose matrix
    void transpose()
    {
	for (int i = 0; i < Dim; ++i) {
	    for (int j = 0; j < Dim; ++j) {
		std::swap(m_[i][j], m_[j][i]);
	    }
	}
    }

    /// Similar to glRotate(), sets the matrix to be a rotation of
    /// angle radians about the axis given by (x, y, z).
    /// (x, y, z) does not have to be normalized.
    /// Only makes sense for Dim == 3 or Dim == 4.
    void setToRotation(T angle, T x, T y, T z);

    /// Sets the matrix to be a rotation that 
    /// takes the point p to q.
    /// p and q are assumed to lie on the 2-sphere.
    /// Only makes sense for Dim == 3.
    void setToRotation(const Vector3D& p, const Vector3D& q);

    /// Matrix-matrix multiplication
    MatrixXD operator * (const MatrixXD& other) const
    {
	MatrixXD res;
	T inner;
	for (int i = 0; i < Dim; ++i) {
	    for (int j = 0; j < Dim; ++j) {
		inner = 0;
		for (int k = 0; k < Dim; ++k) {
		    inner += m_[i][k] * other.m_[k][j];
		}
		res.m_[i][j] = inner;
	    }
	}
	return res;
    }

    ///  Multiply this matrix with a matrix
    MatrixXD& operator *= (const MatrixXD& other)
    {
	(*this) = (*this) * other;
	return *this;
    }

    /// Multiplication by a scalar
    MatrixXD operator * (T scalar) const
    {
	MatrixXD dup(*this);
	dup *= scalar;
	return dup;
    }

    /// Multiplication by a scalar
    MatrixXD& operator *= (T scalar)
    {
	for (int i = 0; i < Dim; ++i) {
	    for (int j = 0; j < Dim; ++j) {
		m_[i][j] *= scalar;
	    }
	}
	return *this;
    }

    /// Addition
    MatrixXD operator + (const MatrixXD& other) const
    {
	MatrixXD dup(*this);
	dup += other;
	return dup;
    }

    /// Add a matrix to this matrix
    MatrixXD& operator += (const MatrixXD& other)
    {
	for (int i = 0; i < Dim; ++i) {
	    for (int j = 0; j < Dim; ++j) {
		m_[i][j] += other.m_[i][j];
	    }
	}
	return *this;
    }

    /// Negation
    MatrixXD operator - () const
    {
	MatrixXD dup(*this);
	for (int i = 0; i < Dim; ++i) {
	    for (int j = 0; j < Dim; ++j) {
		dup.m_[i][j] = -m_[i][j];
	    }
	}
	return dup;
    }

    /// Matrix-vector multiplication
    template <class VectorType>
    VectorType operator * (const VectorType& vec) const
    {
	VectorType res(vec); // Make a copy so that it gets the right size.
	T inner;
	for (int i = 0; i < Dim; ++i) {
	    inner = 0;
	    for (int k = 0; k < Dim; ++k) {
		inner += m_[i][k] * vec[k];
	    }
	    res[i] = inner;
	}
	return res;
    }

    /// Matrix-vector multiplication
    template <class FromConstIteratorType, class ToIteratorType>
    void mult (const FromConstIteratorType from,
	       const ToIteratorType to) const
    {
	T inner;
	for (int i = 0; i < Dim; ++i) {
	    inner = 0;
	    for (int k = 0; k < Dim; ++k) {
		inner += m_[i][k] * from[k];
	    }
	    to[i] = inner;
	}
    }

    /// Determinant
    T det() const;

    /// Trace
    T trace() const
    {
	T tr(0.0);
	for (int i = 0; i < Dim ; ++i) {
	    for (int j = 0; j < Dim ; ++j) {
		tr +=m_[i][j];
	    }
	}
	return tr;
    }

    /// Frobenius norm
    T frobeniusNorm() const
    {
	T fn(0.0);
	for (int i = 0; i < Dim ; ++i) {
	    for (int j = 0; j < Dim ; ++j) {
		fn +=m_[i][j]*m_[i][j];
	    }
	}
	fn = std::sqrt(fn);
	return fn;
    }

    /// Submatrix with given row and column removed.
    MatrixXD<T, Dim-1> submatrix(int r, int c) const
    {
	MatrixXD<T, Dim-1> subm;
	for (int j = 0; j < Dim; ++j) {
	    if (j == r) continue;
	    int joff = 0;
	    if (j > r) joff = -1;
	    for (int k = 0; k < Dim; ++k) {
		if (k == c) continue;
		int koff = 0;
		if (k > c) koff = -1;
		subm(j + joff, k + koff) = m_[j][k];
	    }
	}
	return subm;
    }

private:
    T m_[Dim][Dim];
};


    template <typename T, int Dim>
    inline void MatrixXD<T,Dim>::setToRotation(T angle, T x, T y, T z)
    {
	static_assert(Dim == 3 || Dim == 4, "Expected Dim == 3 or 4");
	THROW("This code should never be entered!");
    }

    template <>
    inline void MatrixXD<double,3>::setToRotation(double angle, double x, double y, double z)
    {
	Array<double, 3> u(x, y, z);
	u.normalize();
	MatrixXD<double,3> S;
	S(0,0) = S(1,1) = S(2,2) = 0.0;
	S(0,1) = -u[2];
	S(1,0) = u[2];
	S(0,2) = u[1];
	S(2,0) = -u[1];
	S(1,2) = -u[0];
	S(2,1) = u[0];
	MatrixXD<double,3> uut;
	uut(0,0) = u[0]*u[0];
	uut(0,1) = u[0]*u[1];
	uut(0,2) = u[0]*u[2];
	uut(1,0) = u[1]*u[0];
	uut(1,1) = u[1]*u[1];
	uut(1,2) = u[1]*u[2];
	uut(2,0) = u[2]*u[0];
	uut(2,1) = u[2]*u[1];
	uut(2,2) = u[2]*u[2];
	// Now, make this matrix into
	// uut + cos(angle)*(I-uut) + sin(angle)*S;
	double cosang = std::cos(angle);
	double sinang = std::sin(angle);
	uut *= (double(1.0) - cosang);
	S *= sinang;
	identity();
	(*this) *= cosang;
	(*this) += uut;
	(*this) += S;
    }


    template <>
    inline void MatrixXD<double,4>::setToRotation(double angle, double x, double y, double z)
    {
	identity();
	MatrixXD<double,3> r;
	r.setToRotation(angle, x, y, z);
	m_[0][0] = r(0,0);
	m_[0][1] = r(0,1);
	m_[0][2] = r(0,2);
	m_[1][0] = r(1,0);
	m_[1][1] = r(1,1);
	m_[1][2] = r(1,2);
	m_[2][0] = r(2,0);
	m_[2][1] = r(2,1);
	m_[2][2] = r(2,2);
    }

    template <typename T, int Dim>
    inline void MatrixXD<T,Dim>::setToRotation(const Vector3D& p,
					       const Vector3D& q)
    {
	static_assert(Dim == 3, "Expected Dim == 3");
	THROW("This code should never be entered!");
    }

    template <>
    inline void MatrixXD<double, 3>::setToRotation(const Vector3D& p,
						   const Vector3D& q)
    {
	Vector3D v = p % q;

// 	double alpha = p.angle(q);
// 	setToRotation(alpha, v[0], v[1], v[2]);

	MatrixXD<double,3> S;
	S(0,0) = S(1,1) = S(2,2) = 0.0;
	S(0,1) = -v[2];
	S(1,0) = v[2];
	S(0,2) = v[1];
	S(2,0) = -v[1];
	S(1,2) = -v[0];
	S(2,1) = v[0];
	MatrixXD<double,3> vvt;
	vvt(0,0) = v[0]*v[0];
	vvt(0,1) = v[0]*v[1];
	vvt(0,2) = v[0]*v[2];
	vvt(1,0) = v[1]*v[0];
	vvt(1,1) = v[1]*v[1];
	vvt(1,2) = v[1]*v[2];
	vvt(2,0) = v[2]*v[0];
	vvt(2,1) = v[2]*v[1];
	vvt(2,2) = v[2]*v[2];
	// Now, make this matrix into
	// vvt + cos(angle)*(I+vvt) + S;
	double cosang = p*q;
	if (cosang+1.0 > -1.0e-13 && cosang+1.0 < 1.0e-13)
	    ;
	else
	    vvt *= 1.0/(1.0 + cosang);
	identity();
	(*this) *= cosang;
	(*this) += vvt;
	(*this) += S;
    }

    /// \cond
    /// Compute determinant of matrix
    template <typename T, int Dim>
    class DetComp
    {
    public:
	T det(const T m[Dim][Dim])
	{
	    // @@ Slow implementation...
	    // Developing along the first coordinate (rows).
	    T result(0);
	    for (int i = 0; i < Dim; ++i) {
		// Make the submatrix
		MatrixXD<T, Dim-1> subm;
		for (int j = 1; j < Dim; ++j) {
		    for (int k = 0; k < Dim; ++k) {
			if (k == i) continue;
			int koff = 0;
			if (k > i) koff = -1;
			subm(j - 1, k + koff) = m[j][k];
		    }
		}
		// Add or subtract the sub determinant
		if (i/2 == (i+1)/2) {
		    result += subm.det()*m[0][i];
		} else {
		    result -= subm.det()*m[0][i];
		}
	    }
	    return result;
	}
    };
    /// \endcond

    /// Specialization of the determinant function for the 1x1 case.
    /// Terminates the det() template recursion.
    //    template<> removed by JAM
    /// \cond
    template<typename T>
    class DetComp<T,1>
    {
    public:
	T det(const T m[1][1])
	{
	    return m[0][0];
	}
    };

    /// Determinant
    template <typename T, int Dim>
    inline T MatrixXD<T,Dim>::det() const
    {
	DetComp<T,Dim> dc;
	return dc.det(m_);
    }

    /// input operator
    template <typename T, int Dim>
    inline std::istream& operator>> (std::istream& is, MatrixXD<T, Dim>& m)
    {
	for (int i = 0; i < Dim; ++i) {
	    for (int j = 0; j < Dim; ++j) {
		is >> m(i,j);
	    }
	}
	return is;
    }
    /// \endcond

    /// output operator
    template <typename T, int Dim>
    inline std::ostream& operator<< (std::ostream& os, const MatrixXD<T, Dim>& m)
    {
	for (int i = 0; i < Dim; ++i) {
	    for (int j = 0; j < Dim; ++j) {
		os <<  m(i,j) << ' ';
	    }
	    os << '\n';
	}
	return os;
    }


} // namespace Go


#endif // _MATRIXXD_H

