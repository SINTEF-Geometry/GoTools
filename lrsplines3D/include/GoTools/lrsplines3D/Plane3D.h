//===========================================================================
//                                                                           
// File: Plane3D.h                                                       
//                                                                           
// Created: Mon Feb 25 11:07:45 2013                                         
//                                                                           
// Author: Sverre Briseid <Sverre.Briseid@sintef.no>
//                                                                           
// Revision: $Id$
//                                                                           
// Description:
//                                                                           
//===========================================================================

#ifndef _PLANE3D_H
#define _PLANE3D_H


// @@sbr201302 Axis-aligned: We should include this feature in the name I guess.
// What about LockedDir.h?

namespace Go
{

    enum Plane3D {XCONST=0, YCONST=1, ZCONST=2};

}; // end namespace Go


#endif // _PLANE3D_H

