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

#ifndef _EXAMPLES_TRIVARIATE_DOXYMAIN_H
#define _EXAMPLES_TRIVARIATE_DOXYMAIN_H

/**
\example coons_patch_volume_gen coons_patch_volume_gen.C
\verbatim
\endverbatim
This program demonstrates the use of the static function 
\em createCoonsPatch
in namespace \em CoonsPatchVolumeGen'
The function can create a new \em SplineVolume representing the coons patch of
six SplineSurfaces, the six faces of the volume.

\example createCoonsVolume createCoonsVolume.C
\verbatim
\endverbatim
The program creates a Coons volume and performs smoothing of this volume
keeping the boundary surfaces fixed.

\example linear_swept_volume linear_swept_volume.C
\verbatim
\endverbatim
This program demonstrates the use of the static function \em linearSweptVolume
in the class \em SweepVolumeCreator.
The function can generate a SplineVolume by sweeping a surface along a curve
or sweeping a curve along a surface. 
A sweeping point on the curve or the surface must be specified.
If the point lies on the the surface, the surface will be swept along the
curve. If the point lies on the the curve, the curve will be swept along the
surface. The curve and the surface must be such that it doesn't lead to
self-intersection.

\example loft_volume_creator loft_volume_creator.C
\verbatim
\endverbatim
This program demonstrates the use of a static function \em loftVolume
in namespace \em LoftVolumeCreator.
The function use lofting to create a new \em SplineVolume based on a set of
surfaces.  The surfaces are not changed during the lofting process.
The surfaces must lie in the same space.

\example rotational_swept_volume rotational_swept_volume.C
\verbatim
\endverbatim
This program demonstrates the use of the static function
\em rotationalSweptVolume in the class \em SweepVolumeCreator.
The function can generate a SplineVolume by rotating a surface around an axis.
The surface must be such that it doesn't lead to self-intersection.

 */

#endif // _EXAMPLES_TRIVARIATE_DOXYMAIN_H
