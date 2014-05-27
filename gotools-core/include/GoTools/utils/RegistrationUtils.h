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

#ifndef _REGISTRATIONUTILS_H
#define _REGISTRATIONUTILS_H


#include <vector>
#include "GoTools/utils/Point.h"


namespace Go
{

  /// Given two sequences of points in 3D, get the rotation, rescaling and translation that sends the second point set
  /// as close as possible to the first (i.e. that minimizes the sum of the square distances). The sequences must be of
  /// same length (at least three), and the points will be matched in the order they come in the vectors, i.e.
  /// points_transform[i] should, after the transformation, be as close as possible to points_fixed[i].
  /// The rotation matrix (together with the rescaling) is stored in 'rotation_matrix', the tronslation (to be performed
  /// after the roation) is stored in 'translate'
  bool registration(const std::vector<Point>& points_fixed, const std::vector<Point>& points_transform,
		    std::vector<std::vector<double> >& rotation_matrix, Point& translate);

} // namespace Go


#endif // _REGISTRATIONUTILS_H

