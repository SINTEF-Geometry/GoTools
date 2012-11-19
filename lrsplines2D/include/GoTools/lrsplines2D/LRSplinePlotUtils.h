//===========================================================================
//                                                                           
// File: LRSplinePlotUtils.h                                                 
//                                                                           
// Created: Fri Nov 16 14:13:35 2012                                         
//                                                                           
// Author: Sverre Briseid <Sverre.Briseid@sintef.no>
//                                                                           
// Revision: $Id$
//                                                                           
// Description:
//                                                                           
//===========================================================================

#ifndef _LRSPLINEPLOTUTILS_H
#define _LRSPLINEPLOTUTILS_H


#include "GoTools/lrsplines2D/LRSplineSurface.h"
#include <iostream>


namespace Go
{

    // Write to file, on PostScript-format, the parametric mesh.
    void writePostscriptMesh(Go::LRSplineSurface& lr_spline_sf, std::ostream &out);

}; // End namespace Go


#endif // _LRSPLINEPLOTUTILS_H

