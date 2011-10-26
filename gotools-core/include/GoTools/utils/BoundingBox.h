//===========================================================================
//                                                                           
// File: BoundingBox.h                                                      
//                                                                           
// Created: August 2000
//                                                                           
// Author: Atgeirr F Rasmussen <atgeirr@sintef.no>
//                                                                           
// Revision: $Id: BoundingBox.h,v 1.17 2009-05-13 07:33:08 vsk Exp $
//                                                                           
//===========================================================================

#ifndef _BOUNDINGBOX_H
#define _BOUNDINGBOX_H

#include "GoTools/utils/Point.h"
#include <vector>
#include "GoTools/utils/config.h"

namespace Go
{


    /** Axis-aligned bounding box.
     *  A BoundingBox object can be an axis-aligned box in any
     *  number of dimensions.
     */

class GO_API BoundingBox
{
public:
    /// The default constructor makes an
    /// uninitialized object, with dimension zero.
    BoundingBox() : valid_(false) {}
    /// Creates a BoundingBox of the specified dimension,
    /// but apart from that still uninitialized.
    explicit BoundingBox(int dim)
	: low_(dim), high_(dim), valid_(false) {}
    /// Creates a BoundingBox with the Point low specifying
    /// the lower bound in all dimensions and high specifying
    /// the upper bound.
    BoundingBox(const Point& low, const Point& high)
	: low_(low), high_(high), valid_(false) { check(); }
    /// Do not inherit from this class -- nonvirtual destructor.
    ~BoundingBox();

    /// Makes the bounding box have lower bounds as specified in low
    /// and upper bounds as specified in high.
    void setFromPoints(const Point& low, const Point& high);

    /// Given an array of dim-dimensional points stored as doubles
    /// or floats, makes the smallest bounding box containing all
    /// points in the array.
    template <typename FloatType>
    void setFromArray(const FloatType* start, const FloatType* end, int dim)
    {
	low_ = Point(start, start+dim);
	high_ = Point(start, start+dim);
	start += dim;
	while (start != end) {
	    for (int d = 0; d < dim; ++d) {
		if (start[d] < low_[d]) {
		    low_[d] = start[d];
		} else if (start[d] > high_[d]) {
		    high_[d] = start[d];
		}
	    }
	    start += dim;
	}

	check();
    }
    /// Given an array of dim-dimensional points stored as doubles
    /// or floats, makes the smallest bounding box containing all
    /// points in the array.
    template <typename ForwardIterator>
    void setFromArray(ForwardIterator start,
		      ForwardIterator end, int dim)
    {
	low_ = Point(start, start+dim);
	high_ = Point(start, start+dim);
	start += dim;
	while (start != end) {
	    for (int d = 0; d < dim; ++d) {
		if (start[d] < low_[d]) {
		    low_[d] = start[d];
		} else if (start[d] > high_[d]) {
		    high_[d] = start[d];
		}
	    }
	    start += dim;
	}

	check();
    }
    /// Given a vector of dim-dimensional points, makes the smallest
    /// bounding box containing all points in the array.
    void setFromPoints(const std::vector<Point>& points)
    {
	int dim = points[0].dimension();
	low_ = points[0];
	high_ = points[0];
	for (size_t i = 1; i < points.size(); ++i) {
	    for (int d = 0; d < dim; ++d) {
		if (points[i][d] < low_[d]) {
		    low_[d] = points[i][d];
		} else if (points[i][d] > high_[d]) {
		    high_[d] = points[i][d];
		}
	    }
	}

	check();
    }

    /// Read a bounding box from a standard istream.
    void read(std::istream& is);
    /// Write a bounding box to a standard ostream.
    void write(std::ostream& os) const;

    /// The dimension of the bounding box.
    int  dimension()  const { return low_.size(); }

    /// The lower bound of the bounding box.
    const Point& low() const { return low_; }
    /// The upper bound of the bounding box.
    const Point& high() const { return high_; }

    /// Returns true if the point pt is inside the
    /// box, or within tol of the boundary.
    bool containsPoint(const Point& pt, double tol = 0.0) const;

    /// Returns true if the two boxes overlap, or are a
    /// distance less than tol apart.
    bool overlaps(const BoundingBox& box, double tol = 0.0) const;
    bool getOverlap(const BoundingBox& box, double& overlap, double tol = 0.0) const;
    /// Returns true if this box contain the box passed as a parameter,
    /// if enlarged by tol in all directions.
    bool containsBox(const BoundingBox& box, double tol = 0.0) const;

    /// After the call, the bounding box will contain
    /// both this box and the given point.
    void addUnionWith(const Point& pt);
    /// After the call, the bounding box will contain
    /// both initial boxes.
    void addUnionWith(const BoundingBox& box);

    /// Is the bounding box initialized?
    bool valid() const { return valid_; }

    /// Check that box validity.
    /// Call valid() to find out if check succeeded.
    void check() const;

private:
    // Data members
    Point low_;
    Point high_;
    mutable bool valid_;
};


} // namespace Go


namespace std {


/// Read BoundingBox from input stream
inline std::istream& operator >> (std::istream& is,
				  Go::BoundingBox& bbox)
{
    bbox.read(is);
    return is;
}


/// Write BoundingBox to output stream
inline std::ostream& operator << (std::ostream& os,
				  const Go::BoundingBox& bbox)
{
    bbox.write(os);
    return os;
}


} // namespace std


#endif // _BOUNDINGBOX_H
