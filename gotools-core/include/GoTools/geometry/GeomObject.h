//===========================================================================
//                                                                           
// File: GeomObject.h                                                      
//                                                                           
// Created: Thu Jun 29 21:31:32 2000                                         
//                                                                           
// Author: Atgeirr F Rasmussen <atgeirr@sintef.no>
//                                                                           
// Revision: $Id: GeomObject.h,v 1.18 2009-05-13 07:30:48 vsk Exp $
//                                                                           
// Description:
//                                                                           
//===========================================================================

#ifndef _GEOMOBJECT_H
#define _GEOMOBJECT_H

#include "GoTools/utils/BoundingBox.h"
#include "GoTools/geometry/Streamable.h"
#include "GoTools/geometry/ClassType.h"
#include "GoTools/utils/config.h"

namespace Go
{

const int MAJOR_VERSION = 1;
const int MINOR_VERSION = 0;

    /** 
     *  Base class for geometrical objects (curves, surfaces, etc.) regrouping
     *  all the properties that they have in common.
     */

class GO_API GeomObject : public Streamable
{
public:
    virtual ~GeomObject();
    
    /// Return the object's bounding box
    virtual BoundingBox boundingBox() const = 0;
    
    /// Return the dimension of the space in which the object lies (usually 2 or 3)
    virtual int dimension() const = 0;

    /// Return the class type identifier of a given, derived instance of GeomObject
    virtual ClassType instanceType() const = 0;

    /// Return the class type identifier of a given class derived from GeomObject
    static ClassType classType();

    /// Clone the GeomObject and return a pointer to the clone.
    virtual GeomObject* clone() const = 0;

    /// Write header information of the GeomObject to stream.  This typically precedes
    /// the act of writing the object itself to a stream, to signal to the receiver
    /// what object is streamed.
    void writeStandardHeader(std::ostream& os) const;
};

} // namespace Go

#endif // _GEOMOBJECT_H




