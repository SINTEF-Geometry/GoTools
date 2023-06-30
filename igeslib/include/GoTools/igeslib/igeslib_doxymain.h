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

#ifndef _IGESLIB-DOXYMAIN_H
#define _IGESLIB-DOXYMAIN_H


/**
\page igeslib GoTools Igeslib

The module depends on the GoTools core library which again depend on SISL.

The \link Go::IGESconverter IGES converter \endlink
read an IGES file and represents its entites in the 
internal data structure of GoTools. It can also write a model represented in
GoTools to an IGES file or convert between an IGES file and the 
\link streamable_doc internal file format \endlink of GoTools.

GoTools core represent only geometric entities. Thus, IGES entities like annotation,
structure, property, associativity, view, drawing 
and figure will be neglected. Neither are constructive solid geometry or
finite element modelling entites handled. If such entities exist in a file
read by the IGES converter, warning messages will be issued.

The topological entities specified in IGES 5.3 is not handled by the current
version of the IGES converter. Thus, the entities vertex, edge, edge list,
loop, face and shell is not handled. However, the geometric entities 
corresponding to these topological entities will be read. Colour information
is read.

The IGES converter reads the following entities:
\li 100 - Circular segment
\li 102 - Composite curve
\li 104 - Conic arc
\li 106, form factor 12 - Copious data, linear path
\li 108 - Plane
\li 110 - Line
\li 116 - Point 
\li 118 - Ruled surface
\li 120 - Surface of revolution
\li 122 - Tabulated cylinder
\li 123 - Direction
\li 124 - Transformation matrix
\li 126 - NURBS curve
\li 128 - NURBS surface
\li 141 - Boundary
\li 142 - Curve on parametric surface
\li 143 - Bounded surface
\li 144 - Trimmed surface
\li 190 - Plane surface
\li 192 - Right cicular surface
\li 314 - Colour description


The content of an IGES file is transferred to the application as a vector
of \link Go::GeomObject GeomObjects\endlink.
By checking the type of each 
object and acting thereafter, the model can be stored and handled in 
the GoTools environment.

The following GoTools entities can be exported using IGES:
\li SplineCurve
\li Spline Surface
\li BoundedSurface

To write an IGES file, the file entities are added one by one to the IGES
convertor using the function addGeom which takes a GeomObject as parameter.
The actual file is written by the command writeIGES. An example showing how
data is read from and written to an IGES file can be found in 
\link intersect_with_plane \endlink.

*/

#endif // _IGESLIB-DOXYMAIN_H
