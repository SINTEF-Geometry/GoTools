//===========================================================================
//                                                                           
// File: RectangularVolumeTesselator.h                                       
//                                                                           
// Created: Thu Jul  5 15:07:09 2012                                         
//                                                                           
// Author: Sverre Briseid <Sverre.Briseid@sintef.no>
//                                                                           
// Revision: $Id$
//                                                                           
// Description:
//                                                                           
//===========================================================================

#ifndef _RECTANGULARVOLUMETESSELATOR_H
#define _RECTANGULARVOLUMETESSELATOR_H


#include "GoTools/tesselator/Tesselator.h"
#include "GoTools/tesselator/RegularMesh.h"
#include "GoTools/tesselator/LineStrip.h"
#include "GoTools/geometry/ParamSurface.h"
#include "GoTools/trivariate/ParamVolume.h"
#include "GoTools/trivariate/RegularVolMesh.h"
#include "GoTools/utils/config.h"

namespace Go
{

/** RectangularSurfaceTesselator: create a mesh for a boundary trimmed parametric
    surface. Visualization purposes.
*/

class GO_API RectangularVolumeTesselator : public Tesselator
{
public:
  /// Constructor. Volume and mesh size are given. The tesselator can be set
  /// to compute also iso parametric curves with a specified mesh size.
 RectangularVolumeTesselator(const ParamVolume& vol,
			     int res = 20)//,
			     // int vres = 20,
			     // int wres = 20)// ,
				   // bool iso = false,
				   // int uiso = 15,
				   // int viso = 15,
				   // int isores = 300)
	: vol_(vol)// ,
	// isolines_(iso), uiso_(uiso), viso_(viso), isores_(isores)
    {
	mesh_ = shared_ptr<RegularVolMesh>(new RegularVolMesh(res, true, true));
    }

  /// Destructor
    virtual ~RectangularVolumeTesselator();
  
    // @@sbr Note: Edge vertices are constructed for each face,
    // meaning that the same vertice may occur more than once in the
    // triangulation, possibly resulting in some small artifacts along
    // edges.
    virtual void tesselate();

    /// Fetch the computed mesh
    shared_ptr<RegularVolMesh> getMesh()
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

    /// Change mesh size
    void changeRes(int m)
    {
	mesh_->resize(m);
	tesselateVolume();
    }

    /// Fetch mesh size
    void getRes(int& m)
    {
	m = (mesh_->numStrips()/6) + 1;
	// m = mesh_->numVertices()/n;
    }

private:
    // We store triangle strips for each of the 6 sides.
    void tesselateVolume();

    const ParamVolume& vol_;
    shared_ptr<RegularVolMesh> mesh_;
};

} // namespace Go


#endif // _RECTANGULARVOLUMETESSELATOR_H

