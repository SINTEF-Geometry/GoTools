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
  if (valid_)
    {
      addUnionWith(box.low_);
      addUnionWith(box.high_);
    }
  else
    {
      low_ = box.low_;
      high_ = box.high_;
      valid_ = true;
    }
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
