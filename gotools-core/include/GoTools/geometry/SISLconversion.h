//===========================================================================
//                                                                           
// File: SISLconversion.h                                                  
//                                                                           
// Created: Wed Nov  8 16:38:00 2000                                         
//                                                                           
// Author: Atgeirr F Rasmussen <atgeirr@sintef.no>
//                                                                           
// Revision: $Id: SISLconversion.h,v 1.4 2009-05-13 07:32:30 vsk Exp $
//                                                                           
// Description:
//                                                                           
//===========================================================================

#ifndef _SISLCONVERSION_H
#define _SISLCONVERSION_H

//#include "sisl.h"
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

