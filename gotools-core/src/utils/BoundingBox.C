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
vector<Point> BoundingBox::lineIntersect(const Point& p1, const Point& dir) const
//===========================================================================
{
  // For each box side, represent the side as an infite plane and intersect with 
  // the line. Check if the found point lies inside the box
  vector<Point> result;
  double tol = 1.0e-6;

  // Bottom
  Point norm(0.0, 0.0, -1.0);
  double div = dir*norm;
  double t;
  if (fabs(div) > tol)
    {
      t = (low_ - p1)*norm/div;
      if (t >= 0.0)
	{
	  Point pos = p1 + t*dir;
	  if (containsPoint(pos, tol))
	    result.push_back(pos);
	}
    }

  // Top
  norm.setValue(0.0, 0.0, 1.0);
  div = dir*norm;
  if (fabs(div) > tol)
    {
      t = (high_ - p1)*norm/div;
      if (t >= 0.0)
	{
	  Point pos = p1 + t*dir;
	  if (containsPoint(pos, tol))
	    result.push_back(pos);
	}
    }

  // Left side
  norm.setValue(-1.0, 0.0, 0.0);
  div = dir*norm;
  if (fabs(div) > tol)
    {
      t = (low_ - p1)*norm/div;
      if (t >= 0.0)
	{
	  Point pos = p1 + t*dir;
	  if (containsPoint(pos, tol))
	    result.push_back(pos);
	}
    }

  // Right side
  norm.setValue(1.0, 0.0, 0.0);
  div = dir*norm;
  if (fabs(div) > tol)
    {
      t = (high_ - p1)*norm/div;
      if (t >= 0.0)
	{
	  Point pos = p1 + t*dir;
	  if (containsPoint(pos, tol))
	    result.push_back(pos);
	}
    }

  // Front side
  norm.setValue(0.0, -1.0, 0.0);
  div = dir*norm;
  if (fabs(div) > tol)
    {
      t = (low_ - p1)*norm/div;
      if (t >= 0.0)
	{
	  Point pos = p1 + t*dir;
	  if (containsPoint(pos, tol))
	    result.push_back(pos);
	}
    }

  // Back side
  norm.setValue(0.0, 1.0, 0.0);
  div = dir*norm;
  if (fabs(div) > tol)
    {
      t = (high_ - p1)*norm/div;
      if (t >= 0.0)
	{
	  Point pos = p1 + t*dir;
	  if (containsPoint(pos, tol))
	    result.push_back(pos);
	}
    }

  return result;
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
