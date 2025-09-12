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

#ifndef _EXAMPLES_LRSPLINES2D_DOXYMAIN_H
#define _EXAMPLES_LRSPLINES2D_DOXYMAIN_H
/**
\page examples_LRSplines2D Example programs related to lrsplines2D

A number of example programs related to lrsplines2D.

refine_lrsurf.C

This program demonstrates definition of an LR spline surface from a
tensor-product spline surface and how to refine this surface with
specified knot line segments.

evaluateLRSurface.C

This program demonstrates the various evaluation possibilities for an 
LR spline surface.

investigate_LRSplineSurface.C

This program demonstrates a set of enquire functionalities for an 
LR spline surface.

investigate_Element2D.C

This program reads an LR spline surface from file.
It iterates trough all elements in the surface and demonstrate available 
enquiries.

approximateWithLRFunc.C

Given a point cloud that can be parameterized on its x- and y-coordinated,
set appropriate parameters and approximate the points with an LR surface
in 1D (function).

approximateParPointsWithLRSurf.C

Given a parameterized point cloud, set a ppropriate parameters and 
approximate the points with a 3D LR surface.

isoContoursLRFunc.C

Compute equally spaced contour curves from LR spline function.

comparePointsLRSurf3D.C

Compute distance between a points in a parameterized point cloud and
the corresponding 3D surface. Collect points into groups according the
this distance.

identify_and_resolve_linear_dependence.C

Identify B-splines involved in a linear dependence situation.
Remove the dependence by applying structured mesh refinement to the
main B-splines in the dependence relation.

\example refine_lrsurf refine_lrsurf.C
\verbatim
\endverbatim

This program demonstrates definition of an LR spline surface from a
tensor-product spline surface and how to refine this surface with
specified knot line segments.

\example evaluateLRSurface evaluateLRSurface.C
\verbatim
\endverbatim

This program demonstrates the various evaluation possibilities for an 
LR spline surface.

\example investigate_LRSplineSurface investigate_LRSplineSurface.C
\verbatim
\endverbatim


This program demonstrates a set of enquire functionalities for an 
LR spline surface.

\example investigate_Element2D investigate_Element2D.C
\verbatim
\endverbatim

This program reads an LR spline surface from file.
It iterates trough all elements in the surface and demonstrate available 
enquiries.

\example approximateWithLRFunc approximateWithLRFunc.C
\verbatim
\endverbatim

Given a point cloud that can be parameterized on its x- and y-coordinated,
set appropriate parameters and approximate the points with an LR surface
in 1D (function).

\example approximateParPointsWithLRSurf approximateParPointsWithLRSurf.C
\verbatim
\endverbatim

Given a parameterized point cloud, set a ppropriate parameters and 
approximate the points with a 3D LR surface.

\example isoContoursLRFunc isoContoursLRFunc.C
\verbatim
\endverbatim

Compute equally spaced contour curves from LR spline function.

\example comparePointsLRSurf3D comparePointsLRSurf3D.C
\verbatim
\endverbatim

Compute distance between a points in a parameterized point cloud and
the corresponding 3D surface. Collect points into groups according the
this distance.

\example identify_and_resolve_linear_dependence identify_and_resolve_linear_dependence.C
\verbatim
\endverbatim

Identify B-splines involved in a linear dependence situation.
Remove the dependence by applying structured mesh refinement to the
main B-splines in the dependence relation.

*/
#endif
