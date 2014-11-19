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

#ifndef _POINTCLOUD_H
#define _POINTCLOUD_H


#include "GoTools/utils/errormacros.h"
#include "GoTools/geometry/GeomObject.h"
#include "GoTools/utils/Array.h"
#include <vector>
#include <assert.h>

namespace Go {

// NB! This class is now a template. The dimension is the template
// parameter. Typedefs are provided for ease of use. Thus, what used
// to be a "PointCloud" is now a "PointCloud3D". @@@jbt

/** 
 * Represent a cloud of points in 'Dim'-dimensional space.
 */
template <int Dim>
class PointCloud : public GeomObject {
public:
    /// Constructor for an empty PointCloud that can only be assigned to
    /// or read(...) into.
    PointCloud()
    {};

    /// start is supposed to point to the start of the array to be copied,
    /// its valuetype should be convertible to double

    /// Constructor specifying a PointCloud from a number of predefined points.
    /// \param start iterator to the beginning of an array where the defining 
    ///        points are stored (these will be copied to the objects internal
    ///        datastructure).
    /// \param numpoints number of points to be read into the PointCloud
    template <typename ForwardIterator>
    PointCloud(ForwardIterator start, int numpoints)
	: points_(numpoints)
    {
#if (!defined (_MSC_VER)  || _MSC_VER > 1599) // Getting rid of warning C4996 on Windows
	std::copy(start, start + Dim * numpoints, points_[0].begin());
#else
	stdext::unchecked_copy(start, start + Dim * numpoints, points_[0].begin());
#endif // _MSC_VER
    }

    /// Constructor specifying a PointCloud from an Array of points.
    PointCloud(std::vector<Array<double, Dim> >& points)
        : points_(points)
    {}

    /// Virtual destructor, enables safe inheritance.
    virtual ~PointCloud()
    {}

    // inherited from GeomObject
    virtual BoundingBox boundingBox() const
    {
	BoundingBox box;
	box.setFromArray(points_[0].begin(), points_.back().end(), Dim);
	return box;
    }

    /// Query the dimension of the space where the PointCloud lies.
    virtual int dimension() const
    { return Dim; }

    // inherited from GeomObject
    virtual ClassType instanceType() const
    { return PointCloud::classType(); }
    
    // inherited from GeomObject
    static ClassType classType()
    { return Class_PointCloud; }

    // inherited from GeomObject
// #ifdef _MSC_VER
// #if _MSC_VER < 1300
//     virtual GeomObject* clone() const
//     { return new PointCloud(*this); }
// #else
//     virtual PointCloud* clone() const
//     { return new PointCloud(*this); }
// #endif // _MSC_VER < 1300
// #else
    virtual PointCloud* clone() const
    { return new PointCloud(*this); }
// #endif

    /// Query the number of points in the PointCloud
    /// \return the number of points in the PointCloud
    int numPoints() const
    { return (int)points_.size(); }
    
    /// Get the coordinates of a point in the point cloud
    /// \param i the index of the requested point
    /// \return a reference to the Array containing the 
    ///         coordinates of the requested point
    Array<double, Dim>& point(int i) 
    { return points_[i]; }

    /// Get the coordinates of a point in the point cloud (const-version)
    /// \param i the index of the requested point
    /// \return a const reference to the Array containing the 
    ///         coordinates of the requested point
    const Array<double, Dim>& point(int i) const
    { return points_[i]; }

    /// Get a pointer to the memory area where the point coordinates are 
    /// (consecutively) stored.
    /// \return a pointer to the coordinate storage area.
    double* rawData()
    { return points_[0].begin(); }

    /// Get a const pointer to the memory area where the point coordinates
    /// are (consecutively) stored.
    /// \return a constant pointer to the coordinate storage area
    const double* rawData() const
    { return points_[0].begin(); }

    /// Get a reference to the vector where the points are stored
    /// \return a reference to the point coordinate vector
    const std::vector<Array<double, Dim> >& pointVector()
    { return points_; }

    void translate(Array<double, Dim> vec)
    {
	for (size_t ki=0; ki<points_.size(); ++ki)
	{
	    points_[ki] += vec;
	}
    }

    // inherited from Streamable
    virtual void read (std::istream& is)
    {
	bool is_good = is.good();
	if (!is_good) {
// The THROW macro results in linker error for some reason (step_reader module).
	    // THROW("Invalid geometry file!");
	    throw std::runtime_error("Invalid geometry file!");
	}
	int nump;
	is >> nump;
	// ALWAYS_ERROR_IF(nump < 1, "Less than one point in cloud.");
	if (nump < 1) {
	    throw std::runtime_error("Less than one point in cloud.");
	}
	is_good = is.good();
	if (!is_good) {
// The THROW macro results in linker error for some reason (step_reader module).
	    // THROW("Invalid geometry file!");
	    throw std::runtime_error("Invalid geometry file!");
	}
	points_.resize(nump);
	for (int i = 0; i < nump; ++i) {
	    for (int d = 0; d < Dim; ++d) {
		is >> points_[i][d];
	    }
	}

	is_good = is.good();
	if (!is_good) {
// The THROW macro results in linker error for some reason (step_reader module).
	    // THROW("Invalid geometry file!");
	    throw std::runtime_error("Invalid geometry file!");
	}
    }

    // inherited from Streamable
    virtual void write (std::ostream& os) const
    {
	os << std::setprecision(15);

	int nump = (int)points_.size();
	os << nump << '\n';
	for(int i = 0; i < nump; ++i) {
	    os << points_[i][0];
	    for (int d = 1; d < Dim; ++d) {
		os << ' ' << points_[i][d];
	    }
	    os << '\n';
	}
	os << std::endl;
    }

    // void translate(const std::vector<double>& dir)
    // {
    // 	int nump = (int)points_.size();
    // 	assert(dir.size() == Dim);
    // 	for(int i = 0; i < nump; ++i) {
    // 	    for (int d = 0; d < Dim; ++d) {
    // 		points_[i][d] += dir[d];
    // 	    }
    // 	}
    // }

private:
    std::vector<Array<double, Dim> > points_;

};


typedef PointCloud<3> PointCloud3D;
typedef PointCloud<4> PointCloud4D;


} // namespace Go




#endif // _POINTCLOUD_H

