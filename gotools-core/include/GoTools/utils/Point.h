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

#ifndef _POINT_H
#define _POINT_H


#include "GoTools/utils/Array.h"
#include "GoTools/utils/errormacros.h"
#include <algorithm>
#include <math.h>
#include "GoTools/utils/config.h"
#include "GoTools/utils/Values.h"


namespace Go {


/** Run-time sized point class.
 *  Encapsulates a point or vector with full value semantics, i.e. copying
 *  and assignment works as expected. The class has some
 *  vector algebra functionality, such as scalar product,
 *  multiplication by scalars etc, and objects will sometimes be
 *  called 'vectors' in the following. Based on double precision floating
 *  point numbers.
 */
class GO_API Point
{
private:
    double* pstart_;
    int n_;
    bool owns_;

public:
    /// Default constructor, does not initialize elements.
    /// The only functions you are allowed to use with a
    /// default constructed (0-dim) Point are the
    /// assignment operator, resize and setValue(...). This is not enforced.
    Point()
	: pstart_(0), n_(0), owns_(true)
    {}
    /// Constructor taking a dimension argument.
    /// Resulting point is of the specified dimension,
    /// but not initialized.
    explicit Point(int dim)
	: pstart_(new double[dim]), n_(dim), owns_(true)
    {}
    /// Constructor taking 2 arguments, makes the
    /// 2D-point (x,y).
    Point(double x, double y)
	: pstart_(new double[2]), n_(2), owns_(true)
    {
	pstart_[0] = x;
	pstart_[1] = y;
    }
    /// Constructor taking 3 arguments, makes the
    /// 3D-point (x,y,z).
    Point(double x, double y, double z)
	: pstart_(new double[3]), n_(3), owns_(true)
    {
	pstart_[0] = x;
	pstart_[1] = y;
	pstart_[2] = z;
    }

    /// Constructor making a Point as a copy of an Array.
    /// The Array must have a value type convertible to double.
    // @@@ 'typename T' must be present in order to compile on MS
    // v6.0. (Using just 'int Dim' and 'Array<double, Dim>' will cause
    // "internal compiler error" !) @jbt
    template <typename T, int Dim>
    explicit Point(const Array<T, Dim>& v)
	: pstart_(0), n_(Dim), owns_(true)
    {
	pstart_ = new double[n_];
#if (!defined (_MSC_VER)  || _MSC_VER > 1599) // Getting rid of warning C4996 on Windows
	std::copy(v.begin(), v.end(), pstart_);
#else
	stdext::unchecked_copy(v.begin(), v.end(), pstart_);
#endif // _MSC_VER
    }

// #ifndef _MSC_VER
    /// Constructor making a Point from an iterator range.
    template <typename RandomAccessIterator>
    Point(RandomAccessIterator first, RandomAccessIterator last)
	: pstart_(0), n_((int)(last - first)), owns_(true)
    {
	pstart_ = new double[n_];
#if (!defined (_MSC_VER)  || _MSC_VER > 1599) // Getting rid of warning C4996 on Windows
	std::copy(first, last, pstart_);
#else
	stdext::unchecked_copy(first, last, pstart_);
#endif // _MSC_VER
    }
// #endif
// #ifdef _MSC_VER
// #if _MSC_VER > 1200
//     template <typename ForwardIterator>
//     Point(ForwardIterator first, ForwardIterator last)
// 	: pstart_(0), n_(last - first), owns_(true)
//     {
// 	pstart_ = new double[n_];
// 	std::copy(first, last, pstart_);
//     }
// #else // Old version
//     Point(const double* first, const double* last)
// 	: pstart_(0), n_(last - first), owns_(true)
//     {
// 	pstart_ = new double[n_];
// 	std::copy(first, last, pstart_);
//     }
//     Point(const float* first, const float* last)
// 	: pstart_(0), n_(last - first), owns_(true)
//     {
// 	pstart_ = new double[n_];
// 	std::copy(first, last, pstart_);
//     }
// #endif // _MSC_VER > 1200
// #endif

    /// Make a point from an existing range of doubles.
    /// If the own parmeter is false, the Point will refer
    /// to the input data, and not own them. If it is true,
    /// it will act as a regular point, owning its data.
    Point(double* begin, double* end, bool own)
	: pstart_(0), n_((int)(end-begin)), owns_(own)
    {
	if (owns_) {
	    pstart_ = new double[n_];
#if (!defined (_MSC_VER)  || _MSC_VER > 1599) // Getting rid of warning C4996 on Windows
	    std::copy(begin, end, pstart_);
#else
	    stdext::unchecked_copy(begin, end, pstart_);
#endif // _MSC_VER
	} else {
	    pstart_ = begin;
	}
    }

    /// Copy constructor.
    Point(const Point& v)
	: pstart_(0), n_(v.n_), owns_(true)
    {
	pstart_ = new double[n_];
#if (!defined (_MSC_VER)  || _MSC_VER > 1599) // Getting rid of warning C4996 on Windows
	std::copy(v.pstart_, v.pstart_ + n_, pstart_);
#else
	stdext::unchecked_copy(v.pstart_, v.pstart_ + n_, pstart_);
#endif // _MSC_VER
    }

    /// Assignment operator.
    Point& operator = (const Point &v)
    {
	Point temp(v);
	swap(temp);
	return *this;
    }

    /// Destructor.
    ~Point()
    {
	if (owns_) delete [] pstart_;
    }

    /// Swaps two Point instances. Never throws.
    void swap(Point& other)
    {
	std::swap(pstart_, other.pstart_);
	std::swap(n_, other.n_);
	std::swap(owns_, other.owns_);
    }

    /// Reads a Point elementwise from
    /// a standard istream. The Point must already be
    /// initialized with the correct dimension.
    void read(std::istream& is);

    /// Writes a Point elementwise to
    /// a standard ostream. Precision is
    /// set to 16 by this function. Dimension
    /// is not stored.
    void write(std::ostream& os) const;


    /// Read-only index access.
    const double& operator [] (int i) const  { return pstart_[i]; }
    /// Index access.
    double& operator [] (int i)        { return pstart_[i]; }
    /// Get a read-only start iterator.
    const double* begin() const { return pstart_; }
    /// Get a start iterator.
    double* begin()       { return pstart_; }
    /// Get a read-only end iterator.
    const double* end() const { return pstart_ + n_; }
    /// Get a end iterator.
    double* end()       { return pstart_ + n_; }
    /// Get the dimension of the Point.
    int size() const { return n_; }
    /// Get the dimension of the Point. Same as size().
    int dimension() const { return n_; }

    /// Changing dimension. This loses all info in the point.
    void resize(int d)
    {
	if (n_ < d) {
	    Point temp(d);
	    swap(temp);
	} else {
	    n_ = d;
	}
    }

    /// Set function for 2D.
    void  setValue(double x, double y)
    {
	resize(2);
	pstart_[0] = x;
	pstart_[1] = y;
    }
    /// Set function for 3D.
    void  setValue(double x, double y, double z)
    {
	resize(3);
	pstart_[0] = x;
	pstart_[1] = y;
	pstart_[2] = z;
    }
    /// Set function that copies values from
    /// an input array. Dimension must already
    /// be initialized!
    void  setValue(const double* array)
    {
	for (int i = 0; i < n_; ++i) {
	    pstart_[i] = array[i];
	}
    }
    /// Set function that sets all elements equal to
    /// the input. Dimension must already
    /// be initialized!
    void  setValue(double val)
    {
	for (int i = 0; i < n_; ++i) {
	    pstart_[i] = val;
	}
    }

    void resetValue(int idx, double val)
	{
	    DEBUG_ERROR_IF(idx < 0 || idx >= n_,
		 "Dimension mismatch.");
	    pstart_[idx] = val;
	}

    /// Get the square of the euclidian length
    /// of the vector.
    double length2() const
    {
	double l2 = 0;
	for (int i = 0; i < n_; ++i)
	    l2 += pstart_[i]*pstart_[i];
	return l2;
    }
    /// Get the euclidian length
    /// of the vector.
    double length()  const
    {
	return sqrt(length2());
    }
    /// Get the infinity-norm (or max-norm) length
    /// of the vector.
    double lengthInf() const
    {
	double linf = 0;
	for (int i = 0; i < n_; ++i) {
	    linf = std::max(linf, fabs(pstart_[i]));
	}
	return linf;
    }
    /// Get the square of the euclidian length
    /// of the difference between this vector
    /// and another vector.
    double dist2(const Point &v) const
    {
	double l2 = 0;
	double d;
	for (int i = 0; i < n_; ++i) {
	    d = pstart_[i] - v.pstart_[i];
	    l2 += d*d;
	}
	return l2;
    }
    /// Get the euclidian length
    /// of the difference between this vector
    /// and another vector.
    double dist(const Point &v) const
    {
	return sqrt(dist2(v));
    }
    /// Get the infinity-norm (or max-norm) length
    /// of the difference between this vector
    /// and another vector.
    double distInf(const Point &v) const
    {
	double linf = 0;
	double d;
	for (int i = 0; i < n_; ++i) {
	    d = pstart_[i] - v.pstart_[i];
	    linf = std::max(linf, fabs(d));
	}
	return linf;
    }
    /// Normalize this vector, i.e. divide every element
    /// by length().
    void normalize()
    {
	double tl = length();
	DEBUG_ERROR_IF(tl == 0.0, "Cannot normalize vector of zero length");
	(*this) /= tl;
    }

    double normalize_checked() 
    {
	double tl = length(); 
	if (tl<1.0e-12)
	    return 0;

	(*this) /= tl; 
	return tl;
    }

    /// The sum of two vectors.
    Point  operator + (const Point     &v) const
    {
	Point res(*this);
	res += v;
	return res;
    }
    /// Add a vector to this vector.
    void     operator +=(const Point     &v)
    {
	DEBUG_ERROR_IF(n_!=v.n_,
		 "Dimension mismatch.");
	for (int i = 0; i < n_; ++i)
	    pstart_[i] += v.pstart_[i];
    }

    /// The difference between two vectors.
    Point  operator - (const Point     &v) const
    {
	Point res(*this);
	res -= v;
	return res;
    }
    /// Subtract a vector from this vector.
    void     operator -=(const Point     &v)
    {
	DEBUG_ERROR_IF(n_!=v.n_,
		 "Dimension mismatch.");
	for (int i = 0; i < n_; ++i)
	    pstart_[i] -= v.pstart_[i];
    }
    /// The product of a vector and a scalar.
    Point  operator * (double d) const
    {
	Point res(*this);
	res *= d;
	return res;
    }
    /// Multiply this vector by a scalar.
    void     operator *=(double d)
    {
	for (int i = 0; i < n_; ++i)
	    pstart_[i] *= d;
    }

    /// A vector divided by a scalar.
    Point  operator / (double d) const
    {
	Point res(*this);
	res /= d;
	return res;
    }
    /// Divide this vector with a scalar.
    void     operator /=(double d)
    {
	for (int i = 0; i < n_; ++i)
	    pstart_[i] /= d;
    }

    /// The negation of a vector.
    Point  operator - () const
    {
	Point res(*this);
	for (int i = 0; i < n_; ++i)
	    res.pstart_[i] = - pstart_[i];
	return res;
    }

    /// The scalar product (or inner product, or dot
    /// product) of two vectors.
    double   operator * (const Point     &v) const
    {
	DEBUG_ERROR_IF(n_!=v.n_,
		 "Dimension mismatch.");
	double res = 0;
	for (int i = 0; i < n_; ++i)
	    res += pstart_[i]*v.pstart_[i];
	return res;
    }

    /// The cross product of two vectors.
    /// Dimensions should be 3.
    Point  operator % (const Point     &v) const
    {
	DEBUG_ERROR_IF(n_ != 3 || v.n_ != 3,
		 "Dimension mismatch.");
	return Point(pstart_[1]*v.pstart_[2] - pstart_[2]*v.pstart_[1],
		     pstart_[2]*v.pstart_[0] - pstart_[0]*v.pstart_[2],
		     pstart_[0]*v.pstart_[1] - pstart_[1]*v.pstart_[0]);
    }

    /// The cross product of two vectors.
    /// Throws if dimensions are not 3.
    Point cross(const Point     &v) const
    {
	return operator%(v);
    }

    /// Set this Point to the cross product of two other Points.
    void setToCrossProd(const Point &u, const Point &v)
    {
	DEBUG_ERROR_IF(u.n_!=v.n_,
		 "Dimension mismatch.");
	DEBUG_ERROR_IF(u.n_!=3,
		 "Dimension must be 3.");

	bool have_already = owns_ && (n_ >= v.n_);
	if (!have_already) {
	    Point temp(3);
	    swap(temp);
	} else {
	    n_ = 3;
	}
	pstart_[0]=u.pstart_[1]*v.pstart_[2] - u.pstart_[2]*v.pstart_[1];
	pstart_[1]=u.pstart_[2]*v.pstart_[0] - u.pstart_[0]*v.pstart_[2];
	pstart_[2]=u.pstart_[0]*v.pstart_[1] - u.pstart_[1]*v.pstart_[0];
    }

    /// The cosine of the angle between this and
    /// another vector.
    double cosAngle(const Point& v) const
    {
	double tl1 = length();
	double tl2 = v.length();
	//(DEBUG_ERROR_IF(tl1*tl2 == 0.0, "Vector of zero length");
	if (tl1*tl2 == 0.0) 
	{
	    //MESSAGE("Vector of zero length in angle compuation");
	    return 0.0;
	}
	double res = ((*this)*v)/(tl1*tl2);
	res = std::min(1.0, std::max(-1.0, res));
	return res;
    }

    /// The angle between this and another vector. Range is [0.0, M_PI].
    double angle(const Point& v) const
    {
	return acos(cosAngle(v));
    }

    /// The angle between this and another vector of dimension 2. Range is [0.0, 2*M_PI].
    double angle2(const Point& v) const
    {
	DEBUG_ERROR_IF(n_!=v.n_,
		 "Dimension mismatch.");
	DEBUG_ERROR_IF(v.n_!=2,
		 "Dimension must be 2.");
	double ang = acos(cosAngle(v)); // Range is [0, M_PI).
	double cross_prod = pstart_[0]*v.pstart_[1] - pstart_[1]*v.pstart_[0];
	if (cross_prod < 0.0)
	    ang = 2*M_PI - ang;
	return ang;
    }

     /// The smallest angle between this and another vector regardless
    /// of orientation. Range is [0.0, 0.5*M_PI].
    double angle_smallest(const Point& v) const
    {
	double tl1 = length();
	double tl2 = v.length();
	//DEBUG_ERROR_IF(tl1*tl2 == 0.0, "Vector of zero length");
	double tcos = fabs(((*this)*v)/(tl1*tl2));
	tcos = std::min(1.0, tcos);
	return acos(tcos);
    }

};

    /// The product of a vector and a scalar.
    inline Point operator * (double d, const Point& p)
    { return p*d; }

    /// Stream extraction for Point.
    inline std::istream& operator>>(std::istream& is, Go::Point& v)
    { v.read(is); return is; }

    /// Stream insertion for Point.
    inline std::ostream& operator<<(std::ostream& os, const Go::Point& v)
    { v.write(os); return os; }

    /// Less than operator
    inline bool operator<(const Point& p1, const Point& p2)
    {
	const int dim = p1.dimension();
	DEBUG_ERROR_IF(p2.dimension() != dim, "Dimension Mismatch");
	for (int i = dim-1; i >= 0; --i) {
	    if (p1[i] > p2[i]) return false;
	    if (p1[i] < p2[i]) return true;
	}
	return false;
    }

    /// Equal operator
    inline bool operator==(const Point& p1, const Point& p2)
    {
	const int dim = p1.dimension();
	DEBUG_ERROR_IF(p2.dimension() != dim, "Dimension Mismatch");
	for (int i = 0; i < dim; ++i) {
	    if (p1[i] != p2[i]) return false;
	}
	return true;
    }

} // End of namespace Go



#endif // _POINT_H





