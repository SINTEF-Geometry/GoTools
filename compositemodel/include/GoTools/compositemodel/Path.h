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

#ifndef _PATH_H
#define _PATH_H

#include "GoTools/compositemodel/Vertex.h"
#include "GoTools/compositemodel/ftEdge.h"
#include <vector>


namespace Go
{
  /// Functions related to sequence of edges
  namespace Path
  {
      /// Estimate mid point, normal and radius defined by an edge
    /// sequence
    bool estimateHoleInfo(const std::vector<ftEdge*>& edges, Point& centre, 
			  Point& axis, double& radius, double& angle);

    /// Classify vertices in a path of edges as corner or non-corner
    /// depending on a given tolerance
    void classifyCorners(const std::vector<ftEdge*>& edges, double tol,
			 std::vector<shared_ptr<Vertex> >& corner,
			 std::vector<shared_ptr<Vertex> >& non_corner);

    /// Identify a loops starting and ending in a given vertex in an ordered
    /// sequence of edges
    std::vector<ftEdge*> identifyLoop(std::vector<ftEdge*> edges, 
				      shared_ptr<Vertex> vx);

    void closestPoint(std::vector<ftEdge*> edges, const Point& pt, 
		      int& clo_ind, double& clo_par, 
		      Point& clo_pt, double& clo_dist);  

    /// Combine edges into 4 curves, preferably with splits in real
    /// corners
    void getEdgeCurves(std::vector<ftEdge*>& loop, 
		       std::vector<shared_ptr<ParamCurve> >& space_cvs,
		       std::vector<Point>& joint_points,
		       double eps, double tol,
		       bool corner_in_Tjoint = true);

    /// Extract edge chain with no joints between more than two edges
    /// or corners
    std::vector<ftEdge*> edgeChain(ftEdge *edg, double angtol, 
				   shared_ptr<Vertex>& v1, shared_ptr<Vertex>& v2);

    /// Extract edge chain between two vertices
    /// face - The edges must belong to this face
    /// vx1 - Start vertice in chain
    /// vx2 - End vertice in chain
    /// out_vx - Face vertex not belonging to the chain
    std::vector<ftEdge*> edgeChainFace(ftSurface *face, shared_ptr<Vertex> vx1,
				       shared_ptr<Vertex> vx2, 
				       shared_ptr<Vertex> out_vx);

    /// Identify the face common to a set of edges, if any
    ftSurface* identifyCommonFace(std::vector<ftEdge*>& edgs);

}  // namespace Path

}  // namespace Go


#endif // _PATH_H
