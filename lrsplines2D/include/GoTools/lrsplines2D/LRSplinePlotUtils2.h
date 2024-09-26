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

#ifndef LRSPLINEPLOT_UTILS2_H
#define LRSPLINEPLOT_UTILS2_H

#include "GoTools/lrsplines2D/Mesh2D.h"
//#include "BSplineFunction.h"

namespace Go
{

// NB: Plot functions below require locale to be properly set.
//     To do this, include <locale> in the main program 
//     compilation unit, and execute:
//     setlocale(LC_ALL, "en_US.UTF.8");
void plot_mesh(const Mesh2D& m, int thick_threshold = 2); // plots mesh to standard output
void plot_largest_tensorgrid(const Mesh2D& m); // extract/plot largest possible tensor grid within m.
void plot_supports_at_corner(const Mesh2D& m, int xpos, int ypos, int xdeg, int ydeg);
void plot_support(const Mesh2D& m, const int* const kvec1, const int* const kvec2, int len1, int len2);
void plot_all_supports(const Mesh2D& m, int x_deg, int y_deg);
void plot_rect_domain(const Mesh2D& m, int xmin, int ymin, int xmax, int ymax);
void plot_history(const Mesh2D& m); 
  //void plot_bspline_function(const Mesh2D& m, const BSplineFunction& b);

// Identify the largest tensorgrid found within a general Mesh (when run twice, once in each direction)
std::vector<int> orig_tensorgrid_knotpositions(const Mesh2D& m, Direction2D d);

}; // end namespace Go

#endif
