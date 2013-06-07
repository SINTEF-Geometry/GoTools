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

#ifndef CURVETESSELATOR_H
#define CURVETESSELATOR_H

#include "GoTools/tesselator/Tesselator.h"
#include "GoTools/tesselator/LineStrip.h"
#include "GoTools/geometry/ParamCurve.h"

namespace Go
{

/** Tesselate curve to produce a linear approximation of the curve for
visualization purposed
 */
class GO_API CurveTesselator : public Tesselator
{
public:
  /// Constructor
    CurveTesselator(const ParamCurve& curve)
	: curve_(curve) 
	{
	    mesh_ = shared_ptr<LineStrip>(new LineStrip(500));
	}

  /// Destructor
    virtual ~CurveTesselator();
  
    /// Perform tesselation. The density depends on the resolution
    virtual void tesselate();

    /// Fetch the linear approximation represented as a LineStrip
    shared_ptr<LineStrip> getMesh()
    {
	return mesh_;
    }

    //
    // 010430: I'm not sure if we're going to need these functions, an
    //         alternative is to trigger these actions when asking for
    //         pointers to the discretizations, which will be done by
    //         the 'painter' when that one is requested to redraw the scene...
    //         (jon)
    //

    /// Set tesselation resolution (number of points to evaluate)
    void changeRes(int n)
    {
	mesh_->resize(n);
	tesselate();
    }
    /// Fetch info about resolution
    void getRes(int& n)
    {
	n = mesh_->numVertices();
    }

private:
    const ParamCurve& curve_;
    shared_ptr<LineStrip> mesh_;
};

} // namespace Go

#endif // CURVETESSELATOR_H

