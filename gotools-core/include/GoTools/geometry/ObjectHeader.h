//===========================================================================
//                                                                           
// File: ObjectHeader.h                                                    
//                                                                           
// Created: Wed Nov  8 14:13:28 2000                                         
//                                                                           
// Author: Atgeirr F Rasmussen <atgeirr@sintef.no>
//                                                                           
// Revision: $Id: ObjectHeader.h,v 1.8 2009-05-13 07:30:49 vsk Exp $
//                                                                           
// Description:
//                                                                           
//===========================================================================

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
    ClassType classType() { return class_type_; }

    /// Get the major version number stored in this ObjectHeader
    int majorVersion() { return major_version_; }

    /// Get the minor version number stored in this ObjectHeader
    int minorVersion() { return minor_version_; }

    /// Get the size of the auxiliary data stored in this ObjectHeader 
    /// (size measured in number of ints).
    int auxdataSize() { return (int)auxillary_data_.size(); }

    /// Get a certain piece of auxiliary data (an integer).
    /// \param i the requested integer's position in the auxiliary data vector
    int auxdata(int i) { return auxillary_data_[i]; }

private:
    ClassType class_type_;
    int major_version_;
    int minor_version_;
    std::vector<int> auxillary_data_;
};


} // namespace Go


#endif // _OBJECTHEADER_H

