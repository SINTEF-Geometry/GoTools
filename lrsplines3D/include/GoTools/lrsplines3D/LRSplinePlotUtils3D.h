//===========================================================================
//                                                                           
// File: LRSplinePlotUtils3D.h                                                 
//                                                                           
// Created: Wed 17 2013                                         
//                                                                           
// Author: Sverre Briseid <Sverre.Briseid@sintef.no>
//                                                                           
// Revision: $Id$
//                                                                           
// Description:
//                                                                           
//===========================================================================

#ifndef _LRSPLINEPLOTUTILS3D_H
#define _LRSPLINEPLOTUTILS3D_H


#include "GoTools/lrsplines3D/LRSplineVolume.h"
#include <iostream>


namespace Go
{

    // Write to file all element grid lines, in the parameter domain.
    void writeElementLineCloud(Go::LRSplineVolume& lr_spline_vol, std::ostream &out);

    // Write to file, on PostScript-format, the parametric mesh.
    void writePostscriptMesh(Go::LRSplineVolume& lr_spline_vol, std::ostream &out);

}; // End namespace Go


#endif // _LRSPLINEPLOTUTILS3D_H

