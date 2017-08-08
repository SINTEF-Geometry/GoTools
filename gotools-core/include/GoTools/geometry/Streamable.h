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

#ifndef _STREAMABLE_H
#define _STREAMABLE_H

#include <iostream>
#include <iomanip>
#include "GoTools/utils/errormacros.h"
#include "GoTools/utils/config.h"

namespace Go
{

    /** 
     *  Base class for streamable objects, i.e., objects which can be read from
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

