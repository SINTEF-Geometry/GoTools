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

#ifndef _ROTATEDBOX_H
#define _ROTATEDBOX_H

#include "GoTools/utils/errormacros.h"
#include "GoTools/utils/CompositeBox.h"
#include "GoTools/utils/MatrixXD.h"
#include "GoTools/utils/config.h"
#include <memory>

namespace Go
{


    /** A rotated version of CompositeBox.
     *  It works in the same way, except that the boxes are
     *  aligned with an arbitrary (given) coordinate system.
     */

class RotatedBox
{
public:
    /// Given an array of dim-dimensional points and a
    /// coordinate system given by (axis[0], axis[1], axis[0] x axis[1]),
    /// construct a rotated box containing all points.
    /// If in 2D, coordinate system is (axis[0], Rot(Pi/2)*axis[0])
    template <typename RandomAccessIterator>
    RotatedBox(RandomAccessIterator start,
	       int dim,
	       int num_u,
	       int num_v,
	       const Point* axis)
    {
	// Make coordinate system.
	setCs(axis, dim);
	// Make box.
	setFromArray(start, dim, num_u, num_v);
    }

    /// Creates a RotatedBox with the Point low specifying
    /// the lower bound in all dimensions and high specifying
    /// the upper bound. The inner and edge boxes are equal.
    RotatedBox(const Point& low, const Point& high,
	       const Point* axis)
    {
	// Make coordinate system.
	setCs(axis, low.size());
	// Make box.
	setFromPoints(low, high);
    }

    /// Do not inherit from this class -- nonvirtual destructor.
    ~RotatedBox()
    {
    }

    /// Given an array of dim-dimensional points stored as doubles
    /// or floats, makes the smallest composite box containing all
    /// points in the array. The array must be like a control point
    /// grid, with num_u points in the fastest running direction,
    /// and num_v points in the other direction. For curves, use
    /// num_v = 1.
    template <typename RandomAccessIterator>
    void setFromArray(RandomAccessIterator start,
		      int dim,
		      int num_u,
		      int num_v)
    {
	// Make a temporary array of points
	std::vector<double> pts(dim*num_u*num_v);
	std::copy(start, start + dim*num_u*num_v, &pts[0]);
	// Transform the points.
	if (dim == 2) {
	    for (int i = 0; i < num_u*num_v; ++i) {
		Point p(&pts[0] + i*2,
			&pts[0] + (i+1)*2, false);
		p = cs2_*p;
		pts[i*2] = p[0];
		pts[i*2+1] = p[1];
	    }
	} else if (dim == 3) {
	    for (int i = 0; i < num_u*num_v; ++i) {
		Point p(&pts[0] + i*3,
			&pts[0] + (i+1)*3, false);
		p = cs3_*p;
		pts[i*3] = p[0];
		pts[i*3+1] = p[1];
		pts[i*3+2] = p[2];
	    }
	} else {
	    THROW("Only supports 2 and 3 dimensions.");
	}
	// Make the composite box.
	box_.reset(new CompositeBox(&pts[0], dim, num_u, num_v));
    }

    /// Makes the bounding box have lower bounds as specified in low
    /// and upper bounds as specified in high.
    void setFromPoints(const Point& low, const Point& high)
    {
	int dim = low.size();
	// Transform the points.
	if (dim == 2) {
	    Point p1 = cs2_*low;
	    Point p2 = cs2_*high;
	    // Make the composite box.
	    box_.reset(new CompositeBox(p1, p2));
	} else if (dim == 3) {
	    Point p1 = cs3_*low;
	    Point p2 = cs3_*high;
	    // Make the composite box.
	    box_.reset(new CompositeBox(p1, p2));
	} else {
	    THROW("Only supports 2 and 3 dimensions.");
	}
    }

    /// The dimension of the rotated box.
    int dimension()  const
    {
	return box_->dimension();
    }

    /// The lower bound of the bounding box, in coordinate system of
    /// rotated box.
    const Point& low() const
    {
	return box_->edge().low();
    }

    /// The upper bound of the bounding box, in coordinate system of
    /// rotated box.
    const Point& high() const
    {
	return box_->edge().high();
    }

    /// The lower bound of the bounding box, in standard coordinate
    /// system.
    const Point low_rot() const
    {
	Point low = box_->edge().low();
	Point rot_low;
	if (dimension() == 2) {
	    rot_low = cs2_i_*low;
	} else if (dimension() == 3) {
	    rot_low = cs3_i_*low;
	} else {
	    THROW("Only supports 2 and 3 dimensions.");
	}

	return rot_low;
    }

    /// The upper bound of the bounding box, in standard coordinate
    /// system.
    const Point high_rot() const
    {
	Point high = box_->edge().high();
	Point rot_high;
	if (dimension() == 2) {
	    rot_high = cs2_i_*high; // We need the inverste rotation ...
	} else if (dimension() == 3) {
	    rot_high = cs3_i_*high;
	} else {
	    THROW("Only supports 2 and 3 dimensions.");
	}

	return rot_high;
    }

    /// The composite box. WARNING: The coordinates of this box must
    /// be interpreted in the coordinate system given by coordsystem().
    const CompositeBox& box() const
    {
	return *box_;
    }

    /// Returns true if the point pt is inside the
    /// box, up to tolerances. Tolerances may be specified
    /// separately for inner and edge boxes.
    bool containsPoint(const Point& pt,
		       double toli = 0.0,
		       double tole = 0.0) const
    {
	int dim = dimension();
	Point rotp;
	if (dim == 2) {
	    rotp = cs2_*pt;
	} else if (dim == 3) {
	    rotp = cs2_*pt;
	} else {
	    THROW("Only supports 2 and 3 dimensions.");
	}
	return box_->containsPoint(rotp, toli, tole);
    }

    /// Returns true if the two boxes overlap, up to 
    /// tolerances. Tolerances may be specified
    /// separately for inner and edge boxes.
    bool overlaps(const RotatedBox& box,
		  double toli = 0.0,
		  double tole = 0.0) const
    {
	MatrixXD<double, 2> m2 = cs2_;
	m2 += -box.cs2_;
	if (m2.frobeniusNorm() > 1e-12) {
	    THROW("The two RotatedBox objects have different coordinate systems.");
	}
	MatrixXD<double, 3> m3 = cs3_;
	m3 += -box.cs3_;
	if (m3.frobeniusNorm() > 1e-12) {
	    THROW("The two RotatedBox objects have different coordinate systems.");
	}
	return box_->overlaps(box.box(), toli, tole);
    }

    /// Returns true if this box contain the box passed 
    /// as a parameter,up to tolerances. Tolerances may be
    /// specified separately for inner and edge boxes.
    bool containsBox(const RotatedBox& box,
		     double toli = 0.0,
		     double tole = 0.0) const
    {
	MatrixXD<double, 2> m2 = cs2_;
	m2 += -box.cs2_;
	if (m2.frobeniusNorm() > 1e-12) {
	    THROW("The two RotatedBox objects have different coordinate systems.");
	}
	MatrixXD<double, 3> m3 = cs3_;
	m3 += -box.cs3_;
	if (m3.frobeniusNorm() > 1e-12) {
	    THROW("The two RotatedBox objects have different coordinate systems.");
	}
	return box_->containsBox(box.box(), toli, tole);
    }


private:
    void setCs(const Point* axis, int dim)
    {
	// What we actually compute is the inverse coordinate system.
	// This way when we multiply a vector with cs we get its
	// coordinates in the system given by (axis[0], Rot(Pi/2)*axis[0]).
	// Initiate the other array to zero to avoid an exception
	if (dim == 2) {
	    cs2_(0,0) = axis[0][0];
	    cs2_(0,1) = axis[0][1];
	    cs2_(1,0) = -axis[0][1];
	    cs2_(1,1) = axis[0][0];

	    // Setting the inverse matrix.
	    double det = cs2_(0,0)*cs2_(1,1) - cs2_(0,1)*cs2_(1,0);
	    double det_inv = 1.0/det;
	    cs2_i_(0,0) = det_inv*cs2_(1,1);
	    cs2_i_(0,1) = -det_inv*cs2_(0,1);
	    cs2_i_(1,0) = -det_inv*cs2_(1,0);
	    cs2_i_(1,1) = det_inv*cs2_(0,0);

	    cs3_(0,0) = cs3_(0,1) = cs3_(0,2) = 0.0;
	    cs3_(1,0) = cs3_(1,1) = cs3_(1,2) = 0.0;
	    cs3_(2,0) = cs3_(2,1) = cs3_(2,2) = 0.0;
	    cs3_i_(0,0) = cs3_i_(0,1) = cs3_i_(0,2) = 0.0;
	    cs3_i_(1,0) = cs3_i_(1,1) = cs3_i_(1,2) = 0.0;
	    cs3_i_(2,0) = cs3_i_(2,1) = cs3_i_(2,2) = 0.0;
	} else if (dim == 3) {
	    Point zaxis = axis[0] % axis[1];
	    zaxis.normalize();
	    cs3_(0,0) = axis[0][0];
	    cs3_(0,1) = axis[0][1];
	    cs3_(0,2) = axis[0][2];
	    cs3_(1,0) = axis[1][0];
	    cs3_(1,1) = axis[1][1];
	    cs3_(1,2) = axis[1][2];
	    cs3_(2,0) = zaxis[0];
	    cs3_(2,1) = zaxis[1];
	    cs3_(2,2) = zaxis[2];

	    // Setting the inverse matrix.
	    double det = cs3_(0,0)*(cs3_(2,2)*cs3_(1,1) - cs3_(2,1)*cs3_(1,2)) -
		cs3_(1,0)*(cs3_(2,2)*cs3_(0,1) - cs3_(2,1)*cs3_(0,2)) +
		cs3_(2,0)*(cs3_(1,2)*cs3_(0,1) - cs3_(1,1)*cs3_(0,2));
	    double det_inv = 1.0/det;
	    cs3_i_(0,0) = det_inv*(cs3_(2,2)*cs3_(1,1) - cs3_(2,1)*cs3_(1,2));
	    cs3_i_(0,1) = -det_inv*(cs3_(2,2)*cs3_(0,1) - cs3_(2,1)*cs3_(0,2));
	    cs3_i_(0,2) = det_inv*(cs3_(1,2)*cs3_(0,1) - cs3_(1,1)*cs3_(0,2));
	    cs3_i_(1,0) = -det_inv*(cs3_(2,2)*cs3_(1,0) - cs3_(2,0)*cs3_(1,2));
	    cs3_i_(1,1) = det_inv*(cs3_(2,2)*cs3_(0,0) - cs3_(2,0)*cs3_(0,2));
	    cs3_i_(1,2) = -det_inv*(cs3_(1,2)*cs3_(0,0) - cs3_(1,0)*cs3_(0,2));
	    cs3_i_(2,0) = det_inv*(cs3_(2,1)*cs3_(1,0) - cs3_(2,0)*cs3_(1,1));
	    cs3_i_(2,1) = -det_inv*(cs3_(2,1)*cs3_(0,0) - cs3_(2,0)*cs3_(0,1));
	    cs3_i_(2,2) = det_inv*(cs3_(1,1)*cs3_(0,0) - cs3_(1,0)*cs3_(0,1));

	    cs2_(0,0) = cs2_(0,1) = 0.0;
	    cs2_(1,0) = cs2_(1,1) = 0.0;
	    cs2_i_(0,0) = cs2_i_(0,1) = 0.0;
	    cs2_i_(1,0) = cs2_i_(1,1) = 0.0;
	} else {
	    THROW("Only supports 2 and 3 dimensions.");
	}
    }

    // Data members
    shared_ptr<CompositeBox> box_;
    MatrixXD<double, 2> cs2_;
    MatrixXD<double, 3> cs3_;
    MatrixXD<double, 2> cs2_i_; // We also store the inverse matrices.
    MatrixXD<double, 3> cs3_i_;
};

} // namespace Go



#endif // _ROTATEDBOX_H

