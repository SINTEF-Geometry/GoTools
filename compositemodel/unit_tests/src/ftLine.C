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

#include "GoTools/compositemodel/ftLine.h"
#include "GoTools/utils/BoundingBox.h"



namespace Go
{

//===========================================================================
ftLine::~ftLine()
{
}

//===========================================================================
void ftLine::getTwoPlanes(ftPlane& plane1, ftPlane& plane2) const

// Return two planes intersecting in the current line
//
//===========================================================================
{
    Point vec = dir_;
    vec.normalize();

    Point dum(0.0, 0.0, 0.0);
    if (vec[0] < vec[1] && vec[0] < vec[2])
	dum[0] = 1.0;
    else if (vec[1] < vec[2])
	dum[1] = 1.0;
    else 
	dum[2] = 1.0;

    Point norm1, norm2;
    norm1 = vec.cross(dum);
    norm1.normalize();

    norm2 = vec.cross(norm1);
    norm2.normalize();

    plane1 = ftPlane(norm1, point_);
    plane2 = ftPlane(norm2, point_);

}


// Determine if this line comes close enough to a given point pt after
// a rescaling of pt's coordinates

//===========================================================================
bool ftLine::closeToScaledPoint(const Point& pt,      // Input point before scaling
				const Point& scale,    // Scale factor, all directions
				double dist2) const   // Square of maximum distance
//===========================================================================
{

  // Vector from point_ to the scaled point
  int dim = point_.dimension();
  Point pToS(dim);
  for (int ki=0; ki<dim; ++ki)
    pToS[ki] = pt[ki]*scale[ki] - point_[ki];

  double scProd = pToS * dir_;
  return pToS*pToS - (scProd * scProd) / (dir_ * dir_) <= dist2;
}



//===========================================================================
ftLine ftLine::scaled(const Point& scale) const
//===========================================================================
{
  int dim = point_.dimension();
  Point pnt(dim), dir(dim);
  for (int ki=0; ki<dim; ++ki)
    {
      pnt[ki] = point_[ki]*scale[ki];
      dir[ki] = dir_[ki]*scale[ki];
    }
      
  return ftLine (dir, pnt);
}



//===========================================================================
bool ftLine::intersectsSphereOfBox(const BoundingBox& box) const
//===========================================================================
{
  Point diagonal = box.high() - box.low();
  int dim = point_.dimension();
  Point scale(dim);
  for (int ki=0; ki<dim; ++ki)
    scale[ki] = 1.0/diagonal[ki];
  ftLine scaleLine = scaled(scale);
  return scaleLine.closeToScaledPoint( (box.high() + box.low()) / 2,
				       scale,
				       0.75);    // 0.75 is square of distance from center to corner in unit box
}



//===========================================================================
bool ftLine::planesIntersectBox(const BoundingBox& box) const
//===========================================================================
{
  ftPlane pl1(Point(0.0, 0.0, 0.0), Point(0.0, 0.0, 0.0));   // No empty constructur exists for ftPlane
  ftPlane pl2(pl1);
  getTwoPlanes(pl1, pl2);
  return pl1.intersectsBox(box) && pl2.intersectsBox(box);
}



//===========================================================================
bool ftLine::intersectsBox(const BoundingBox& box) const
//===========================================================================
{
  // Checking is done by projection down to xy, xz or yz-plane. Choose the
  // plane where the line direction is 'steepest'. This is done to avoid zero
  // division in special cases.

  // Suppose projection is done on xy-plane (similar arguments apply for xy
  // and xz), and let p0 (resp. p1) be the projection of the intersection
  // between the line and the plane containing the bottom (resp. top) face of
  // the box. Then the line intersects the box if and only if the line
  // segment between p0 and p1 intersects the projection of the box

  int dim = point_.dimension();
  Point low = box.low();
  Point high = box.high();
  Point diagonal = high-low;

  // Get best projection direction
  int proj_dir;
  double abs_dir[3];
  for (int i = 0; i < 3; ++i)
    abs_dir[i] = dir_[i] >= 0 ? dir_[i] : -dir_[i];
  if (abs_dir[0] > abs_dir[1] && abs_dir[0] > abs_dir[2])
    proj_dir = 0;
  else if (abs_dir[1] > abs_dir[2])
    proj_dir = 1;
  else
    proj_dir = 2;

  int plane_idx[2];
  plane_idx[0] = proj_dir == 0 ? 1 : 0;
  plane_idx[1] = proj_dir == 2 ? 1 : 2;
  if (dim == 2)
    {
      plane_idx[0] = 0;
      plane_idx[1] = 1;
    }

  // Project some points
  Point low2D(low[plane_idx[0]], low[plane_idx[1]]);
  Point high2D(high[plane_idx[0]], high[plane_idx[1]]);
  Point point2D(point_[plane_idx[0]], point_[plane_idx[1]]);
  Point dir2D(dir_[plane_idx[0]], dir_[plane_idx[1]]);

  // Get intersection points
  double scalar = (dim == 2) ? 0 : (low[proj_dir] - point_[proj_dir]) / dir_[proj_dir];
  Point p0 = point2D + (dir2D * scalar);
  scalar = (dim == 2) ? 0 : (high[proj_dir] - point_[proj_dir]) / dir_[proj_dir];
  Point p1 = point2D + (dir2D * scalar);

  // Check some special cases where the line segment dose not intersect,
  // but the extended line might do
  // This is when both p0 and p1 are on the same side of the box
  if (dim == 3)
    {
      if (p0[0] < low2D[0] && p1[0] < low2D[0]) return false;
      if (p0[0] > high2D[0] && p1[0] > high2D[0]) return false;
      if (p0[1] < low2D[1] && p1[1] < low2D[1]) return false;
      if (p0[1] > high2D[1] && p1[1] > high2D[1]) return false;
    }

  // Now the line segment hits the box if and only if its extended line does.
  // This is the same as if the corner points of the box are separated by
  // the line

  Point dir2Dnorm(dir_[plane_idx[1]], -dir_[plane_idx[0]]);
  double scProd = dir2Dnorm * (p0 - low2D);
  if (scProd * (dir2Dnorm * (p0 - high2D)) <= 0) return true;
  if (scProd * (dir2Dnorm * (p0 - Point(low2D[0], high2D[1]))) <= 0) return true;
  if (scProd * (dir2Dnorm * (p0 - Point(high2D[0], low2D[1]))) <= 0) return true;

  return false;

}

} // namespace Go
