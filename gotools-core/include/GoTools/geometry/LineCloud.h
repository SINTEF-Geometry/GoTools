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

#ifndef _LINECLOUD_H
#define _LINECLOUD_H


#include "GoTools/geometry/GeomObject.h"
#include "GoTools/utils/Array.h"
#include <vector>
#include "GoTools/utils/config.h"


namespace Go
{

    /** GeomObject representing a collection of line segments in space.
     *
     */

class GO_API LineCloud : public GeomObject
{
public:

    /// Makes an unitialized LineCloud that can later be assigned or read() into.
    LineCloud()
    {};

    /** start is supposed to point to the start of the array to be copied,
	its valuetype should be convertible to double */

    /// Generate a LineCloud based on values stored in memory.  The lines should
    /// be stored as a sequence of point pairs indicating the start and end position
    /// of each line in the line cloud.  Each point is stored as (x_coord, y_coord, z_coord).
    /// The coordinates should be convertible to 'double'.  Ex: p0_start_x, p0_start_y, 
    /// p0_start_z, p0_end_x, p0_end_y, p0_end_z, p1_start_x, p1_start_y, p1_start_z,....
    /// \param start pointer to the start of the memory area from which point information
    ///              is to be copied
    /// \param numlines number of lines in the LineCloud
    template <typename ForwardIterator>
    LineCloud(ForwardIterator start, int numlines)
	: points_(numlines*2)
    {
	std::copy(start, start + 6*numlines, points_[0].begin());
    }

    /// Virtual destructor, enables safe inheritance.
    virtual ~LineCloud();

    /// Read a LineCloud from an input stream
    /// \param is the stream from which we will read the LineCloud
    virtual void read (std::istream& is);

    /// Write the LineCloud to stream
    /// \param os the stream that we will write the LineCloud to
    virtual void write (std::ostream& os) const;

    /// Get the BoundingBox of the LineCloud
    /// \return the BoundingBox enclosing the LineCloud
    virtual BoundingBox boundingBox() const;

    /// Query the dimension of the space in which the LineCloud is embedded 
    /// (currently, only dimension=3 is allowed...)
    /// \return the dimension of the space (always 3 for now)
    virtual int dimension() const
    { return 3; }

    /// Get the ClassType of this GeomObject (which is of course the ClassType identified
    /// LineCloud).
    /// \return the ClassType identifier of LineCloud
    virtual ClassType instanceType() const
    { return LineCloud::classType(); }

    /// Get the ClassType identifier of LineCloud.
    /// \return the ClassType identifier of LineCloud
    static ClassType classType()
    { return Class_LineCloud; }

    /// Clone this object
    /// \return a pointer to a cloned LineCloud.
// #ifdef _MSC_VER
// #if _MSC_VER < 1300
//     virtual GeomObject* clone() const
//     { return new LineCloud(*this); }
// #else
//     virtual LineCloud* clone() const
//     { return new LineCloud(*this); }
// #endif // _MSC_VER < 1300
// #else
    virtual LineCloud* clone() const
    { return new LineCloud(*this); }
// #endif

    /// Fill a line cloud with information read from memory.  The layout of the
    /// read information should be as for the LineCloud constructor: LineCloud(ForwardIterator 
    /// start, int numlines)
    /// \param points pointer to the memory area where the coordinates of the elements
    ///               in the LineCloud can be found. (This information will be copied).
    /// \param numlines number of lines in the LineCloud.
    void setCloud(const double* points, int numlines);

    /// Query the number of lines in the LineCloud
    /// \return the LineCloud's number of lines.
    int numLines() const { return (int)points_.size()/2; }

    /// Get a start or end point from a line in the LineCloud (non-const version)
    /// \param i the index of the start/end point.  If 'i' is pair, then the returned
    ///          point is a start point, else it is an end point.
    /// \return a reference to the requested start/end point
    Vector3D& point(int i) { return points_[i]; }

    /// Get a start or end point from a line in the LineCloud (const version)
    /// \param i the index of the start/end point.  If 'i' is even, then the returned
    ///          point is a start point, else it is an end point.
    /// \return a const-reference to the requested start/end point
    const Vector3D& point(int i) const { return points_[i]; }

    /// Get a pointer to the start of the internal memory area where line information
    /// is stored.
    /// \return a pointer to the beginning of the array where line information is stored.
    double* rawData() { return points_[0].begin(); }

private:
    std::vector<Vector3D> points_;
};


} // namespace Go



#endif // _LINECLOUD_H

