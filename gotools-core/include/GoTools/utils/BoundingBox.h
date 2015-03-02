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

    /// Compute the intersections with this box and a line
    std::vector<Point> lineIntersect(const Point& p1, const Point& dir) const;

    /// Is the bounding box initialized?
    bool valid() const { return valid_; }

    /// Unset bounding box
    void unset()
      {
	valid_ = false;
	low_.resize(0);
	high_.resize(0);
      }

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
