//===========================================================================
//                                                                           
// File: BoundingBox.C                                                     
//                                                                           
// Created: Wed Jun  7 19:45:05 2000                                         
//                                                                           
// Author: Atgeirr F Rasmussen <atgeirr@sintef.no>
//                                                                           
// Revision: $Id: BoundingBox.C,v 1.15 2006-02-13 13:24:00 jbt Exp $
//                                                                           
// Description:
//                                                                           
//===========================================================================


#include "GoTools/utils/BoundingBox.h"
#include <iostream>

using namespace Go;
using namespace std;

//===========================================================================
BoundingBox::~BoundingBox()
//===========================================================================
{
}

//===========================================================================
void BoundingBox::setFromPoints(const Point& low, const Point& high)
//===========================================================================
{
    low_ = low;
    high_ = high;
    check();
}

//===========================================================================
void BoundingBox::read(std::istream& is)
//===========================================================================
{
    ALWAYS_ERROR_IF(low_.dimension() == 0,
		    "Boundingbox has no set dimension yet - cannot read.");

    is >> low_ >> high_;
    check();
}


//===========================================================================
void BoundingBox::write(std::ostream& os) const
//===========================================================================
{
    ALWAYS_ERROR_IF(!valid_, "Not initialized - cannot write.");
    streamsize prev = os.precision(15);
    os << low_ << endl << high_;
    os.precision(prev);   // Reset precision to it's previous value
}


//===========================================================================
bool BoundingBox::containsPoint(const Point& pt, double tol) const
//===========================================================================
{
    ALWAYS_ERROR_IF(!valid_, "Not initialized - cannot call.");
    for (int d = 0; d < low_.dimension(); ++d) {
	if (pt[d] < low_[d] - tol)
	    return false;
	if (pt[d] > high_[d] + tol)
	    return false;
    }
    return true;
}


//===========================================================================
bool BoundingBox::overlaps(const BoundingBox& box, double tol) const
//===========================================================================
{
    double overlap;
    return getOverlap(box, overlap, tol);
}

//===========================================================================
bool BoundingBox::getOverlap(const BoundingBox& box, double& overlap,
			     double tol) const
//===========================================================================
{
  int kd;
  double t2, t3;
  overlap = 1.0e10;  // Just a large number
  Point other_high = box.high();
  Point other_low = box.low();
  for (kd=0; kd<low_.dimension(); kd++)
    {
      if (high_[kd] > other_high[kd])
	{
	  t2 = low_[kd];
	  t3 = other_high[kd];
	}
      else
	{
	  t2 = other_low[kd];
	  t3 = high_[kd];
	}
      overlap = std::min(overlap, t3-t2);
      if (t3 < t2 - tol)
	return false;
    }
  return true;
}



//===========================================================================
bool BoundingBox::containsBox(const BoundingBox& box, double tol) const
//===========================================================================
{
    return (containsPoint(box.low(), tol) && containsPoint(box.high(), tol));
}



//===========================================================================
void BoundingBox::addUnionWith(const Point& pt)
//===========================================================================
{
    ALWAYS_ERROR_IF (low_.dimension() != pt.dimension(),
		     "Dimension mismatch.");

    if (!valid_) {
	low_ = pt;
	high_ = pt;
	valid_ = true;
	return;
    }
    for (int d = 0; d < low_.dimension(); ++d) {
	if (pt[d] < low_[d])
	    low_[d] = pt[d];
	else if (pt[d] > high_[d])
	    high_[d] = pt[d];
    }
}


//===========================================================================
void BoundingBox::addUnionWith(const BoundingBox& box)
//===========================================================================
{
    addUnionWith(box.low_);
    addUnionWith(box.high_);
}


//===========================================================================
void BoundingBox::check() const
//===========================================================================
{
    int n = low_.dimension();
    if (high_.dimension() != n) {
	THROW("Dimension mismatch.");
    }
    for (int i = 0; i < n; ++i) {
	if (low_[i] > high_[i]) {
	    THROW("Low point is higher than high point.");
	}
    }
    valid_ = true;
}
