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

#ifndef _LRSPLINES2D_DOXYMAIN_H
#define _LRSPLINES2D_DOXYMAIN_H
/**
 \page lrsplines2d GoTools LRsplines2D

 This module contains classes and tools for working with LR spline surfaces.

\subsection Summary

The module lrsplines2D represents LR-spline surfaces.

\subsection The problem with tensor product grids

The simplest way to generate parametric surface patches is by applying the tensor-product construction to univariate parametric space curves. Examples of these are B-spline surfaces, B-spline volumes, NURBS surfaces and NURBS volumes, which are implemented in GoTools.

The number of control points of these tensor-product patches grows exponentially with the dimension of the geometric object. For instance, for curves with 10 control points, corresponding tensor-product surfaces will have 100 control points and tensor-product volumes will have 1000 control points. These control points form a rectangular gridded structure, making it impossible to refine the model locally without refining across the entire domain. Hence, in practice, it becomes unfeasible to add sufficient detail to the model where it is needed most.

\subsection LR-splines

Locally Refined B-splines (LR B-splines) aim to solve this problem by providing a mathematical framework, generalizing the tensor-product construction, for refining the model locally. This framework extends to any dimension.

\subsection Sources

    T. Dokken, T. Lyche, K.F. Pettersen. Polynomial splines over locally refined box-partitions, Computer Aided Geometric Design, Volume 30, Issue 3, March 2013, Pages 331--356
**/
#endif // _LRSPLINES2D_DOXYMAIN_H
