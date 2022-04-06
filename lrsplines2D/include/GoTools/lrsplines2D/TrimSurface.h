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

#ifndef TRIMSURFACE_H
#define TRIMSURFACE_H

#include "GoTools/geometry/BoundedSurface.h"
#include <vector>
#include <math.h>

namespace Go
{
  /// Trimming of LR B-spline surfaces with
  /// respect to a corresponding point cloud
  namespace TrimSurface
  {
    /// Create bounded surface from an LR B-spline surface and associated
    /// point cloud. The domain of the bounded surface will approximately
    /// correspond to the domain of the point cloud while the domain of the
    /// initial surface exceeds the point domain.
    /// \param surf surface approximating a point cloud
    /// \param isotrim wether or not to search for a correspondance between
    /// the initial surface boundaries and the point cloud domain
    /// \param points point cloud corresponding to the given surface
    /// Note that the sequence of the points will be changed
    /// \param tightness indicates how close the trimming curve should go the
    /// domain of the point cloud, dense point clouds allow for a tighter bound than
    /// more sparse or unevenly distributed points. A number between 1 and 7.
    /// \param trim_surf the resulting bounded surface
    /// \param only_outer indicates wether only outer trimming should be performed
    /// or also trimming of holes
    bool makeBoundedSurface(shared_ptr<ParamSurface>& surf,
			    bool isotrim[], 
			    std::vector<double>& points,
			    int tightness,
			    shared_ptr<BoundedSurface>& trim_surf,
			    bool only_outer=true);

    /// Given a set of polygons defining the trimming curves (output from TrimUtils), 
    /// define the bounded surface
    bool defineBdSurface(shared_ptr<ParamSurface>& surf,
			 double domain[], bool isotrim[], double eps,
			 std::vector< std::vector<std::vector<double> > >& seqs,
			 shared_ptr<BoundedSurface>& trim_surf);
  };
};

#endif

