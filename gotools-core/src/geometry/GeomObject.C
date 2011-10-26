//===========================================================================
//                                                                           
// File: GeomObject.C                                                      
//                                                                           
// Created: Fri Sep 15 14:57:12 2000                                         
//                                                                           
// Author: Atgeirr F Rasmussen <atgeirr@sintef.no>
//                                                                           
// Revision: $Id: GeomObject.C,v 1.4 2003-05-08 15:37:13 afr Exp $
//                                                                           
// Description:
//                                                                           
//===========================================================================


#include "GoTools/geometry/GeomObject.h"

namespace Go
{


//===========================================================================
GeomObject::~GeomObject()
//===========================================================================
{
}


//===========================================================================
void GeomObject::writeStandardHeader(std::ostream& os) const
//===========================================================================
{
  if (this->instanceType() == Class_SurfaceOnVolume)
    os << "200 " << ' ' << MAJOR_VERSION << ' '
       << MINOR_VERSION << " 0\n";
  else
    os << this->instanceType() << ' ' << MAJOR_VERSION << ' '
       << MINOR_VERSION << " 0\n";
}


} // namespace Go
