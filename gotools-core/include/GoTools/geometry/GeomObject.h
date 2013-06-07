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




