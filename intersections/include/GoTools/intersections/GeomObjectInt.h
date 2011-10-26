//===========================================================================
//                                                                           
// File: GeomObjectInt.h  
//                                                                           
// Created: 
//                                                                           
// Author: vsk
//                                                                           
// Revision: $Id: GeomObjectInt.h,v 1.5 2006-03-08 15:25:40 jbt Exp $
//                                                                           
// Description:
//                                                                           
//===========================================================================

#ifndef _GEOMOBJECTINT_H
#define _GEOMOBJECTINT_H


namespace Go {


/// This class is an abstract base class providing an interface to the
/// "intersection objects". An intersection object is a wrapper around
/// a corresponding geometric object, containing additional
/// functionality, and that is used in the intersection algorithms.

class GeomObjectInt {
public:
    /// Destructor
    virtual ~GeomObjectInt(){};
 
};


} // namespace Go


#endif // _GEOMOBJECTINT_H
