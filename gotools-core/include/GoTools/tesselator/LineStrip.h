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

#ifndef _LINESTRIP_H
#define _LINESTRIP_H

#include "GoTools/tesselator/GeneralMesh.h"
#include <vector>


namespace Go
{

/** LineStrip: Structure for storing values for a line strip, i.e.
result of curve tesselation
 */

class GO_API LineStrip : public GeneralMesh
{
public:
  /// Constructor given size of mesh
    LineStrip(int n = 200);
    /// Destrcutor
    virtual ~LineStrip();

    /// Change mesh size
    void resize(int n);

    /// Number of nodes
    virtual int numVertices()
    { return (int)vert_.size()/3; }

    virtual double* vertexArray() { return &vert_[0]; }
    virtual double* paramArray() { return &param_[0]; }
     virtual int atBoundary(int idx);
     /// Indices for each line strip
    unsigned int* stripArray() { return &strip_[0]; }
    /// Indices for triangles. Not used
    virtual unsigned int* triangleIndexArray() { return &strip_[0]; }

    /// Casting. Return object as line strip.
    virtual LineStrip* asLineStrip();

     // /// Translate all vertices by vert_translation, wrt the geometry.
     // void translate(const std::vector<double>& vert_translation);

private:
    std::vector<double> vert_;
    std::vector<double> param_;
    std::vector<unsigned int> strip_;

    // /// The 3D-translation wrt the geometry.
    // std::vector<double> vert_translation_;
};



} // namespace Go





#endif // end of #ifdef _LINESTRIP_H_
