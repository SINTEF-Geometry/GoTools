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

#ifndef _OBJECTHEADER_H
#define _OBJECTHEADER_H


#include <vector>
#include "GoTools/geometry/ClassType.h"
#include "GoTools/geometry/Streamable.h"

namespace Go
{

    /** An object representing the "header" usually preceeding a GeomObject 
     *  in a stream.  This header contains information about the GeomObject,
     *  and this information can be read by ObjectHeader and accessed by its 
     *  member functions.
     */ 


class GO_API ObjectHeader : public Streamable
{
public:
    /// Default constructor (uninitialized header)
    ObjectHeader()
	: class_type_(Class_Unknown),
	  major_version_(0),
	  minor_version_(0)
    {}

    /// Constructor creating an ObjectHeader that is initialized
    /// with a given ClassType, major and minor version.
    /// \param t the ClassType of the GeomObject that this ObjectHeader
    ///          shall represent.
    /// \param major major version number
    /// \param minor minor version number
    ObjectHeader(ClassType t, int major, int minor)
	: class_type_(t),
	  major_version_(major),
	  minor_version_(minor)
    {}

    /// Constructor creating an ObjectHeader that is initialized with
    /// a given ClassType, major and minor version and auxiliary, class-
    /// specific data.
    /// \param t the ClassType of the GeomObject that this ObjectHeader
    ///          shall represent.
    /// \param major major version number
    /// \param minor minor version number
    /// \param auxdata auxiliary data, specific to the ClassType.
    ObjectHeader(ClassType t, int major, int minor,
		   const std::vector<int>& auxdata)
	: class_type_(t),
	  major_version_(major),
	  minor_version_(minor),
	  auxillary_data_(auxdata)
    {}

    /// Virtual destructor, allowing safe destruction of derived objects.
    virtual ~ObjectHeader();

    /// Read the ObjectHeader from an input stream
    /// \param is the input stream from which the ObjectHeader is read
    virtual void read (std::istream& is);

    /// Write the ObjectHeader to an output stream
    /// \param os the output stream to which the ObjectHeader is written
    virtual void write (std::ostream& os) const;

    /// Get the ClassType stored in this ObjectHeader
    ClassType classType() const { return class_type_; }

    /// Get the major version number stored in this ObjectHeader
    int majorVersion() const { return major_version_; }

    /// Get the minor version number stored in this ObjectHeader
    int minorVersion() const { return minor_version_; }

    /// Get the size of the auxiliary data stored in this ObjectHeader 
    /// (size measured in number of ints).
    int auxdataSize() const { return (int)auxillary_data_.size(); }

    /// Get a certain piece of auxiliary data (an integer).
    /// \param i the requested integer's position in the auxiliary data vector
    int auxdata(int i) const { return auxillary_data_[i]; }

private:
    ClassType class_type_;
    int major_version_;
    int minor_version_;
    std::vector<int> auxillary_data_;
};


} // namespace Go


#endif // _OBJECTHEADER_H

