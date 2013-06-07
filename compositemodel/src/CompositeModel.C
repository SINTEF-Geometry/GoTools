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

#include "GoTools/compositemodel/CompositeModel.h"

namespace Go
{

  //===========================================================================
  CompositeModel::CompositeModel(double gap, double neighbour, double kink, double bend)
    : toptol_(tpTolerances(gap, neighbour, kink, bend))
      //===========================================================================
  {
      closest_idx_ = -1;
  }



  //===========================================================================
  CompositeModel::~CompositeModel()
  //===========================================================================
  {
  }


  //===========================================================================
  void CompositeModel::setTolerances(double gap, double neighbour,
				     double kink, double bend)
  //===========================================================================
  {
    toptol_ = tpTolerances(gap, neighbour, kink, bend);
  }


  //===========================================================================
  double CompositeModel::boxVecDist(const BoundingBox& box, 
				  const Point& vec) const
  //===========================================================================
  {
    Point b = box.low() - vec;
    Point t = vec - box.high();
    double dist2 = 0;
    double d = std::max(0.0, std::max(b[0], t[0]));
    dist2 += d*d;
    d = std::max(0.0, std::max(b[1], t[1]));
    dist2 += d*d;
    d = std::max(0.0, std::max(b[2], t[2]));
    dist2 += d*d;
    return sqrt(dist2);
  }

//===========================================================================
bool
CompositeModel::boxExtreme(const BoundingBox& box, const Point& dir, 
			 const Point& curr_pnt) const
//===========================================================================
{
  // Fetch box corners
  Point corners[8];
  corners[0] = corners[1] = corners[2] = corners[3] = box.low();
  corners[4] = corners[5] = corners[6] = corners[7] = box.high();
  corners[1][0] = corners[7][0];
  corners[2][1] = corners[7][1];
  corners[3][0] = corners[7][0];
  corners[3][1] = corners[7][1];
  corners[4][0] = corners[0][0];
  corners[4][1] = corners[0][1];
  corners[5][1] = corners[0][1];
  corners[6][0] = corners[0][0];

  // For all corners, check if it lies further in the given
  // direction than the input point
  for (int ki=0; ki<8; ki++)
    {
      if (corners[ki]*dir > curr_pnt*dir)
	return true;
    }

  return false;
}

} // namespace Go
