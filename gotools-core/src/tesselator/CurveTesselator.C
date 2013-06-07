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

#include "GoTools/tesselator/CurveTesselator.h"

namespace Go
{


//===========================================================================
CurveTesselator::~CurveTesselator()
//===========================================================================
{
}

//===========================================================================
void CurveTesselator::tesselate()
//===========================================================================
{
  int dim = curve_.dimension();
  Point pt(3);
  int n = mesh_->numVertices();
  for (int i = 0; i < n; ++i) {
    double rt = double(i)/double(n-1);
    double ta = curve_.startparam();
    double tb = curve_.endparam();

    // Ensure a finite extension in case of an unbounded curve
    // 	double limit = 1.0e3;
    // 	ta = std::max(ta, -1.0*limit);
    // 	tb = std::min(tb, 1.0*limit);

    double t = ta*(1.0 - rt) + rt*tb;
    curve_.point(pt, t);
    int j;
    for (j=0; j<dim; ++j)
      mesh_->vertexArray()[i*3+j] = pt[j];
    for (; j<3; ++j)
      mesh_->vertexArray()[i*3+j] = 0.0;
  }
}

} // namespace Go

