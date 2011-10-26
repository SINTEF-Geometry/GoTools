//===========================================================================
//                                                                           
// File: CompositeBox.h                                                      
//                                                                           
// Created: Wed Dec  8 10:57:06 2004                                         
//                                                                           
// Author: Atgeirr F Rasmussen <atgeirr@sintef.no>
//                                                                           
// Revision: $Id: CompositeBox.h,v 1.8 2009-05-13 07:33:08 vsk Exp $
//                                                                           
//===========================================================================

#ifndef _COMPOSITEBOX_H
#define _COMPOSITEBOX_H

#include "GoTools/utils/BoundingBox.h"

namespace Go
{


    /** Composite of two bounding boxes. The inner box is supposed to 
     *  be a bounding box for the interior of an object, while the
     *  edge box bounds the outer edges. When requests are made, the
     *  CompositeBox is treated as a single axis-aligned box, but
     *  tolerances can be given for both sub-boxes separately.
     *  A CompositeBox object can have any number of dimensions.
     */

class GO_API CompositeBox
{
public:
    /// Constructor using setFromArray()
    template <typename RandomAccessIterator>
    CompositeBox(RandomAccessIterator start,
		 int dim,
		 int num_u,
		 int num_v)
    {
	obsolete_inner_ = 0;
	setFromArray(start, dim, num_u, num_v);
    }
    /// Creates a CompositeBox with the Point low specifying
    /// the lower bound in all dimensions and high specifying
    /// the upper bound. The inner and edge boxes are equal.
    CompositeBox(const Point& low, const Point& high)
	: inner_(low, high), edge_(low, high), obsolete_inner_(0)
    {
    }
    /// Do not inherit from this class -- nonvirtual destructor.
    ~CompositeBox();

    /// Makes the bounding box have lower bounds as specified in low
    /// and upper bounds as specified in high.
    void setFromPoints(const Point& low, const Point& high)
    {
	inner_.setFromPoints(low, high);
	edge_ = inner_;
	obsolete_inner_ = 0;
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
	if (num_u < 3 || num_v == 2) {
	    // We cannot make a separate inner box. Just make
	    // two identical boxes. The num_v == 1 case is the curve
	    // special case, which is handled.
	    edge_.setFromArray(start, start+dim*num_u*num_v, dim);
	    Point mid = 0.5*(edge_.low() + edge_.high());
	    inner_.setFromPoints(mid, mid);
	    obsolete_inner_ = 0;
	    return;
	}
	// First do the inner box.
	int vstart = (num_v == 1) ? 0 : 1;
	int vend = (num_v == 1) ? num_v : num_v - 1;
	Point il(start + dim*(1 + num_u*vstart),
		 start + dim*(2 + num_u*vstart));
	Point ih(il);
	for (int i = vstart; i < vend; ++i) {
	    for (int j = 1; j < num_u-1; ++j) {
		RandomAccessIterator pnt = start + dim*(j + num_u*i);
		for (int d = 0; d < dim; ++d) {
		    if (pnt[d] < il[d]) {
			il[d] = pnt[d];
		    } else if (pnt[d] > ih[d]) {
			ih[d] = pnt[d];
		    }
		}
	    }
	}
	inner_.setFromPoints(il, ih);

	// Then the edge box.
	Point el(start, start+dim);
	Point eh(el);
	if (num_v == 1) {
	    // We have a curve. Only check the endpoint
	    RandomAccessIterator pnt = start + dim*(num_u - 1);
	    for (int d = 0; d < dim; ++d) {
		if (pnt[d] < el[d]) {
		    el[d] = pnt[d];
		} else if (pnt[d] > eh[d]) {
		    eh[d] = pnt[d];
		}
	    }
	} else {
	    // We have a 2D-array, check all boundary points
	    for (int i = 0; i < num_v; ++i) {
		for (int j = 0; j < num_u; j += num_u-1) {
		    RandomAccessIterator pnt = start + dim*(j + num_u*i);
		    for (int d = 0; d < dim; ++d) {
			if (pnt[d] < el[d]) {
			    el[d] = pnt[d];
			} else if (pnt[d] > eh[d]) {
			    eh[d] = pnt[d];
			}
		    }
		}
	    }
	    for (int j = 0; j < num_u; ++j) {
		for (int i = 0; i < num_v; i += num_v-1) {
		    RandomAccessIterator pnt = start + dim*(j + num_u*i);
		    for (int d = 0; d < dim; ++d) {
			if (pnt[d] < el[d]) {
			    el[d] = pnt[d];
			} else if (pnt[d] > eh[d]) {
			    eh[d] = pnt[d];
			}
		    }
		}
	    }
	}
	edge_.setFromPoints(el, eh);
    }

    /// Read a composite box from a standard istream.
    void read(std::istream& is);
    /// Write a composite box to a standard ostream.
    void write(std::ostream& os) const;

    /// The dimension of the composite box.
    int  dimension()  const
    {
	return inner_.dimension();
    }

    /// The lower bound of the composite box.
    Point low(double toli = 0.0,
	      double tole = 0.0) const;


    /// The upper bound of the composite box.
    Point high(double toli = 0.0,
	       double tole = 0.0) const;


    /// The inner box.
    const BoundingBox& inner() const
    {
	return inner_;
    }
    /// The edge box.
    const BoundingBox& edge() const
    {
	return edge_;
    }

    /// Returns true if the point pt is inside the box, up to
    /// tolerances. Tolerances may be specified separately for inner
    /// and edge boxes.
    bool containsPoint(const Point& pt,
		       double toli = 0.0,
		       double tole = 0.0) const;

    /// Returns true if the two boxes overlap, up to
    /// tolerances. Tolerances may be specified separately for inner
    /// and edge boxes.
    bool overlaps(const CompositeBox& box,
		  double toli = 0.0,
		  double tole = 0.0) const;

    /// Get the overlap between the two boxes. Tolerances may be
    /// specified separately for inner and edge boxes.
    bool getOverlap(const CompositeBox& box,
		  double& overlap,
		  double toli = 0.0,
		  double tole = 0.0) const;

    /// Returns true if this box contain the box passed as a
    /// parameter,up to tolerances. Tolerances may be specified
    /// separately for inner and edge boxes.
    bool containsBox(const CompositeBox& box,
		     double toli = 0.0,
		     double tole = 0.0) const;


private:
    // Data members
    BoundingBox inner_;
    BoundingBox edge_;
    int obsolete_inner_;
};

} // namespace Go


#endif // _COMPOSITEBOX_H

