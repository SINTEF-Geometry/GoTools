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

#ifndef _EXAMPLES_TRIVARIATEMODEL_DOXYMAIN_H
#define _EXAMPLES_TRIVARIATEMODEL_DOXYMAIN_H

/**
\example createMidShip createMidShip.C
\verbatim
\endverbatim
             Build a volume model representing a simplified mid ship,
             hull with stiffeneres in both directions and deck.
             The model is represented as a block structured volume
             model.

\example mirrorAndLoft example_mirrorAndLoft.C
\verbatim
\endverbatim
This program demonstrates how to create a volume model from a set
of spline surfaces. The surfaces must be connected. The surface set
is mirrored around a given plane, and lofting between corresponding 
surfaces is performed.

\example multiPatchSweep multiPatchSweep.C
\verbatim
\endverbatim
The idea of this program is to sweep a set surfaces to create a multi patch 
volume model.
The surface set will be created by the example programs in compositemodel:
createSplitDisc and createBlockStructuredDisc

\example BrepToTrivariate BrepToTrivariate.C 
\verbatim
\endverbatim
The idea of this program is to read a Brep model in g2-format and create
trivariate spline model.
Note that the functionality works only for some classes of Brep models.

\example createLinSweptVol createLinSweptVol.C 
\verbatim
\endverbatim
Create a block structured volume model from a face set
describing a boundary represented solid that may be created by
sweeping a planar face set along a stright line

\example createRotationalVol createRotationalVol.C 
\verbatim
\endverbatim
Create a block structured volume model from a face set describing a  
boundary represented solid representing a rotational object.

\example createVolumeBlocks createVolumeBlocks.C
\verbatim
\endverbatim
Create a block structured volume model from a face set
consisting of possibly trimmed surfaces with arbitrary topology
(no corner-to-corner conditions). The face set represents a boundary
represented solid and is described in a g2-file.

\example linearBrepToTrivariate linearBrepToTrivariate.C
\verbatim
\endverbatim
The idea of this program is to read a Brep model in g2-format and check
if it can be regenerated as a linear sweep. In that case a multi block
trivariate spline model will be created.

 */

#endif // _EXAMPLES_TRIVARIATEMODEL_DOXYMAIN_H
