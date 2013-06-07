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

#ifndef _POINTONEDGE
#define _POINTONEDGE


#include "GoTools/utils/Point.h"

namespace Go
{

class ftEdge;

//===========================================================================
/** PointOnEdge - represents a point lying on an edge. Storage.
 * 
 */
//===========================================================================
class PointOnEdge
{

 public:
  /// Constructor
  /// \param the edge on which the point lies
  /// \param par associated edge parameter
  PointOnEdge(ftEdge* edge, double par);

  /// Destructor
  ~PointOnEdge();
	
  /// Point position
  const Point& position() const { return pt_; }
  /// Associated edge
  ftEdge* edge() const { return edge_; }
  /// Associated edge parameter
  double par() const { return par_; }
    
 private:
  ftEdge* edge_;
  double par_;
  Point pt_;
};

} // namespace Go


#endif // _POINTONEDGE_H

