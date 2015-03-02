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

#ifndef LRSURFSTITCH_H
#define LRSURFSTITCH_H

#include "GoTools/geometry/ParamSurface.h"
#include "GoTools/lrsplines2D/LRSplineSurface.h"
#include "GoTools/lrsplines2D/Direction2D.h"
#include "GoTools/lrsplines2D/Mesh2D.h"
#include "GoTools/lrsplines2D/LRBSpline2D.h"

namespace Go
{
  namespace LRSurfStitch
  {
    int averageCorner(std::vector<std::pair<shared_ptr<ParamSurface>,int> >& sfs,
		      double tol);

    bool averageEdge(shared_ptr<ParamSurface> surf1, int edge1,
		     shared_ptr<ParamSurface> surf2, int edge2, double tol);

    bool averageEdge(shared_ptr<LRSplineSurface> surf1, int edge1,
		     shared_ptr<LRSplineSurface> surf2, int edge2, double tol);

    void fetchEdgeCorners(shared_ptr<LRSplineSurface> surf, int edge,
			  double& u1, double& v1, double& u2, double& v2);

    void extractMissingKnots(std::vector<double>& union_vec, 
			     std::vector<double>& vec,
			     double tol, int order,
			     std::vector<double>& resvec);

    void defineRefinements(const Mesh2D& mesh, Direction2D dir,
			   int edge, int ix, std::vector<double>& knot_vals, 
			   int element_width,
			   std::vector<LRSplineSurface::Refinement2D>& refs); 

    void extractBoundaryBsplines(shared_ptr<LRSplineSurface> surf,
				 int edge,
				 std::vector<LRBSpline2D*>& bsplines);

};
};

#endif
