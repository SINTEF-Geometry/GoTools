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

#ifndef _EXAMPLES_LRSPLINES3D_DOXYMAIN_H
#define _EXAMPLES_LRSPLINES3D_DOXYMAIN_H
/**
\page examples_LRSplines3D Example programs related to lrsplines3D

A number of example programs related to lrsplines3D.

refine_lrvol.C

This program demonstrates definition of an LR spline volume from a
tensor-product spline volume and how to refine this volume with
specified mesh rectangles.

investigate_LRSplineVolume.C

This program demonstrates a set of enquire functionalities for an 
LR spline volume.

investigate_Element3D.C

This program reads an LR spline volume from file.
It iterates trough all elements in the volume and demonstrate available 
enquiries.

Bezier_extraction.C

This program reads an LR spline volume from file and compute the Bezier
coefficients for all patches. The patches are written to file as
tensor-product spline volumes.

evaluateLRVolume.C

This program demonstrates the various evaluation possibilities for an LR spline volume

\example refine_lrvol refine_lrvol.C
\verbatim
\endverbatim

This program demonstrates definition of an LR spline volume from a
tensor-product spline volume and how to refine this volume with
specified mesh rectangles.

\example investigate_LRSplineVolume investigate_LRSplineVolume.C
\verbatim
\endverbatim


This program demonstrates a set of enquire functionalities for an 
LR spline volume.

\example investigate_Element3D investigate_Element3D.C
\verbatim
\endverbatim

This program reads an LR spline volume from file.
It iterates trough all elements in the volume and demonstrate available 
enquiries.

\example Bezier_extraction Bezier_extraction.C
\verbatim
\endverbatim

This program reads an LR spline volume from file and compute the Bezier
coefficients for all patches. The patches are written to file as
tensor-product spline volumes.

\example evaluateLRVolume evaluateLRVolume.C
\verbatim
\endverbatim

This program demonstrates the various evaluation possibilities for an LR spline volume

*/
#endif
