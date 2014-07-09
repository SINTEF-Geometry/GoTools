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

#ifndef _EXAMPLES_COMPOSITEMODEL_DOXYMAIN_H
#define _EXAMPLES_COMPOSITEMODEL_DOXYMAIN_H

/**

\example createSplitDisc createSplitDisc.C 
\verbatim
\endverbatim
This programs creates a face set representing a disc as two trimmed
surfaces: The trimmed disc and a rectangular surface lying inside this
disc. This is a starting point for the creation of a disc as a multi patch
model with spline surfaces and no degeneracies.
The construction uses planar, rectangular surfaces and a truncated sylinder,
but the operations performed using these surfaces, do not depend on that
level of regularity.

\example createBlockStructuredDisc createBlockStructuredDisc.C
\verbatim
\endverbatim
This example file creates a block structured set of spline surfaces 
from a face set
consisting of possibly trimmed surfaces with arbitrary topology
(no corner-to-corner conditions). 

\example createVolumeBoundaries createVolumeBoundaries.C
\verbatim
\endverbatim
This programs creates a set of B-spline surfaces intended as the boundary
surfaces for a spline volume. A number of different methods are used in the
surface construction. Thus, this expample program illustrates some of the
possibilities for surface construction.

\example face2splineset face2splineSet.C
\verbatim
\endverbatim
This program demonstrates how to create a set of spline surfaces,
meeting in a corner-to-corner configuration and with corresponding
coefficients at common boundaries, from one possibly trimmed face
///
The program reads a bounded surface from a file, splits this surface
into several bounded surfaces where each surface has (at most) 4 boundary
curves. Finally, each bounded surface is approximated by a spline surface
within a given tolerance and C0 continuities at common boundaries is
ensured.
*/

#endif // _EXAMPLES_COMPOSITEMODEL_DOXYMAIN_H
