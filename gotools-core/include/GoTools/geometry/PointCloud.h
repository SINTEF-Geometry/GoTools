//===========================================================================
//                                                                           
// File: PointCloud.h                                                      
//                                                                           
// Created: Wed Mar 13 15:11:11 2002                                         
//                                                                           
// Author: Atgeirr F Rasmussen <atgeirr@sintef.no>
//                                                                           
// Revision: $Id: PointCloud.h,v 1.19 2007-12-04 16:12:01 jbt Exp $
//                                                                           
// Description:
//                                                                           
//===========================================================================

#ifndef _POINTCLOUD_H
#define _POINTCLOUD_H


#include "GoTools/utils/errormacros.h"
#include "GoTools/geometry/GeomObject.h"
#include "GoTools/utils/Array.h"
#include <vector>


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

    // inherited from Streamable
    virtual void read (std::istream& is)
    {
	bool is_good = is.good();
	if (!is_good) {
	    THROW("Invalid geometry file!");
	}
	int nump;
	is >> nump;
	ALWAYS_ERROR_IF(nump < 1, "Less than one point in cloud.");
	is_good = is.good();
	if (!is_good) {
	    THROW("Invalid geometry file!");
	}
	points_.resize(nump);
	for (int i = 0; i < nump; ++i) {
	    for (int d = 0; d < Dim; ++d) {
		is >> points_[i][d];
	    }
	}

	is_good = is.good();
	if (!is_good) {
	    THROW("Invalid geometry file!");
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

private:
    std::vector<Array<double, Dim> > points_;

};


typedef PointCloud<3> PointCloud3D;
typedef PointCloud<4> PointCloud4D;


} // namespace Go




#endif // _POINTCLOUD_H

