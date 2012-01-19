//===========================================================================
//                                                                           
// File: ParametricSurfaceTesselator.h                                               
//                                                                           
// Created:
//                                                                           
// Author: Atgeirr F Rasmussen <atgeirr@sintef.no>
//                                                                           
// Revision: $Id: ParametricSurfaceTesselator.h,v 1.5 2009-05-13 07:32:51 vsk Exp $
//                                                                           
// Description:
//                                                                           
//===========================================================================

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

