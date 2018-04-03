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

#ifndef _SISLCONVERSION_H
#define _SISLCONVERSION_H

#include "sislP.h"
#include "GoTools/geometry/SplineCurve.h"
#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/utils/config.h"

/// \file SISLconversion.h
/// Declaration file for a set of free conversion functions
/// between SISL and Spline curves and surfaces.

struct SISLCurve;
struct SISLSurf;

namespace Go
{
/// Convert a SplineCurve to a SISLCurve
/// \param cv the SplineCurve to convert
/// \param copy if 'true', then the generated SISLCurve will have its own
///             copy of the coefficient information contained in 'cv'.
///             Otherwise, it will share this information with 'cv' (ie. 
///             it will only contain a pointer into the corresponding
///             storage array in 'cv'.
/// \return A newly generated SISLCurve that describes the same curve as 
///         'cv'.  The user assumes ownership and is responsible for 
///         cleaning up (which means calling the SISL function freeCurve(...)
///         on the pointer when it should be destroyed).
SISLCurve GO_API *Curve2SISL( const SplineCurve& cv, bool copy = true);

/// Convert a SplineCurve to a rational SISLCurve
/// Arrays are copied
SISLCurve GO_API *Curve2SISL_rat( const SplineCurve& cv);

/// Convert a SISLCurve to a SplineCurve
/// \param cv the SISLcurve to convert
/// \return A newly generated SplineCurve that describes the same curve
///         as 'cv'.  The user assumes ownership and is responsible for 
///         cleaning up by calling the \c delete function.
SplineCurve GO_API *SISLCurve2Go( const SISLCurve* const cv);

/// Convert a SplineSurface to a SISLSurface
/// \param sf the SplineSurface to convert
/// \param copy if 'true', then the generated SISLSurf will have its own
///             copy of the coefficient information contained in 'sf'. 
///             Otherwise, it will share this information with 'sf' (ie.
///             it will only contain a pointer into the corresponding 
///             storage array in 'sf'.  
/// \return a newly generated SISLSurf that describes the same surface
///         as 'sf'.  The user assumes ownership and is responsible for 
///         cleaning up (which means calling the SISL function freeSurf(...)
///         on the pointer when it should be destroyed).    
SISLSurf GO_API *GoSurf2SISL( const SplineSurface& sf, bool copy = true);

/// Convert a SISLSurface to a SplineSurface
/// \param sf the SISLSurf to convert
/// \return A newly generated SplineSurface that describes the same surface
///         as 'sf'.  The user assumes ownership and is responsible for 
///         cleaning up by calling the \c delete function.
SplineSurface GO_API *SISLSurf2Go( SISLSurf* sf);

} // namespace Go

#endif // _SISLCONVERSION_H

