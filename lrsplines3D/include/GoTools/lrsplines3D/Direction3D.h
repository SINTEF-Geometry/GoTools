//===========================================================================
//                                                                           
// File: Direction3D.h                                                       
//                                                                           
// Created: Tue Feb 26 12:41:30 2013                                         
//                                                                           
// Author: Sverre Briseid <Sverre.Briseid@sintef.no>
//                                                                           
// Revision: $Id$
//                                                                           
// Description:
//                                                                           
//===========================================================================

#ifndef _DIRECTION3D_H
#define _DIRECTION3D_H


namespace Go
{

  /// Specifies parameter direction of volume (along first (x) parameter: d = XDIR;
  /// along second (y) parameter: YDIR; along third (z) parameter: ZDIR)
    enum Direction3D {XDIR=0, YDIR=1, ZDIR=2};

}; // end namespace Go


#endif // _DIRECTION3D_H

