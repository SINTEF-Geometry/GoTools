//===========================================================================
//                                                                           
// File: Streamable.h                                                      
//                                                                           
// Created: Thu Aug 24 02:55:53 2000                                         
//                                                                           
// Author: Atgeirr F Rasmussen <atgeirr@sintef.no>
//                                                                           
// Revision: $Id: Streamable.h,v 1.12 2009-05-13 07:30:49 vsk Exp $
//                                                                           
// Description:
//                                                                           
//===========================================================================

#ifndef _STREAMABLE_H
#define _STREAMABLE_H

#include <iostream>
#include <iomanip>
#include "GoTools/utils/errormacros.h"
#include "GoTools/utils/config.h"

namespace Go
{

    /** 
     *  Base class for streamable objects, ie. objects which can be read from 
     *  and written to a stream.
     */

class GO_API Streamable
{
public:
    virtual ~Streamable();

    /// read object from stream
    /// \param is stream from which object is read
    virtual void read (std::istream& is) = 0;
    /// write object to stream
    /// \param os stream to which object is written
    virtual void write (std::ostream& os) const = 0;

    // Exception class
    class EofException{};
};

inline std::istream& operator >> (std::istream& is, Go::Streamable& obj)
{
    ALWAYS_ERROR_IF(is.eof(), "End of file reached. Cannot read.");
    obj.read(is);
    return is;
}

inline std::ostream& operator << (std::ostream& os,
				  const Go::Streamable& obj)
{
    obj.write(os);
    return os;
}

} // namespace Go



#endif // _STREAMABLE_H

