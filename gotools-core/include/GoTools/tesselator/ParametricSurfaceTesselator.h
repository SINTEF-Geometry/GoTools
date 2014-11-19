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

#ifndef PARAMETRICSURFACETESSELATOR_H
#define PARAMETRICSURFACETESSELATOR_H

#include "GoTools/tesselator/Tesselator.h"
#include "GoTools/tesselator/RegularMesh.h"
#include "GoTools/geometry/ParamSurface.h"
#include "GoTools/tesselator/GenericTriMesh.h"
#include <memory>
#include "GoTools/utils/config.h"

namespace Go
{

/** ParametricSurfaceTesselator: create a mesh for a possibly trimmed surface with a suitable
    triangulation.
*/

class GO_API ParametricSurfaceTesselator : public Tesselator
{
public:
  /// Constructor. Surface and mesh size are given. The mesh size relates to 
  /// the underlying surface in the case of bounded surfaces.
    ParametricSurfaceTesselator(const ParamSurface& surf)
	: surf_(surf), m_(20), n_(20)
    {
 	mesh_ = shared_ptr<GenericTriMesh>(new GenericTriMesh(0,0,true,true));
    }

    virtual ~ParametricSurfaceTesselator();

    virtual void tesselate();

    // virtual GeneralMesh* getMesh()
    // {
    // 	return mesh_.get();
    // }

    /// Fetch the resulting mesh
    shared_ptr<GenericTriMesh> getMesh()
    {
	return mesh_;
    }

    /// Change mesh size
    void changeRes(int n, int m);

    /// Fetch info about mesh size
    void getRes(int& n, int& m)
    {
	m = m_;
	n = n_;
    }

private:
    const ParamSurface& surf_;
    shared_ptr<GenericTriMesh> mesh_;
    int m_;
    int n_;

};

} // namespace Go




#endif //  PARAMETRICSURFACETESSELATOR_H

