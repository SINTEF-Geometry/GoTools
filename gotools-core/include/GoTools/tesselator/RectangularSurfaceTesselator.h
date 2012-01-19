//===========================================================================
//                                                                           
// File: RectangularSurfaceTesselator.h                                               
//                                                                           
// Created: Wed Nov 28 16:38:31 2001                                         
//                                                                           
// Author: Atgeirr F Rasmussen <atgeirr@sintef.no>
//                                                                           
// Revision: $Id: RectangularSurfaceTesselator.h,v 1.3 2009-05-13 07:32:51 vsk Exp $
//                                                                           
// Description:
//                                                                           
//===========================================================================

#ifndef RECTANGULARSURFACETESSELATOR_H
#define RECTANGULARSURFACETESSELATOR_H

#include "GoTools/tesselator/Tesselator.h"
#include "GoTools/tesselator/RegularMesh.h"
#include "GoTools/tesselator/LineStrip.h"
#include "GoTools/geometry/ParamSurface.h"
#include "GoTools/utils/config.h"

namespace Go
{

/** RectangularSurfaceTesselator: create a mesh for a boundary trimmed parametric
    surface. Visualization purposes.
*/

class GO_API RectangularSurfaceTesselator : public Tesselator
{
public:
  /// Constructor. Surface and mesh size are given. The tesselator can be set
  /// to compute also iso parametric curves with a specified mesh size.
 RectangularSurfaceTesselator(const ParamSurface& surf,
				   int ures = 20,
				   int vres = 20,
				   bool iso = false,
				   int uiso = 15,
				   int viso = 15,
				   int isores = 300)
	: surf_(surf),
	isolines_(iso), uiso_(uiso), viso_(viso), isores_(isores)
    {
	mesh_ = shared_ptr<RegularMesh>(new RegularMesh(ures, vres, true, true));
    }

  /// Destructor
    virtual ~RectangularSurfaceTesselator();
  
    virtual void tesselate();

    /// Fetch the computed mesh
    shared_ptr<RegularMesh> getMesh()
    {
	return mesh_;
    }

    /// Fetch tesselation of iso parametric curves
    std::vector<LineStrip>& getIsolineStrips()
    {
	return isolinestrips_;
    }

    //
    // 010430: I'm not sure if we're going to need these functions, an
    //         alternative is to trigger these actions when asking for
    //         pointers to the discretizations, which will be done by
    //         the 'painter' when that one is requested to redraw the scene...
    //         (jon)
    //

    /// Change mesh size
    void changeRes(int m, int n)
    {
	mesh_->resize(m, n);
	tesselateSurface();
    }

    /// Fetch mesh size
    void getRes(int& m, int& n)
    {
	n = mesh_->numStrips() + 1;
	m = mesh_->numVertices()/n;
    }

    /// Retesselate iso parametric curves
    void changeIsolines(bool isolines)
    {
	isolines_ = isolines;
	if (isolines_) {
	    tesselateIsolines();
	}
    }

    /// Fetch mesh of iso parametric curves
    void getIsolines(bool& isolines)
    {
	isolines = isolines_;
    }

    /// Change the number of iso parametric curves to compute
    void changeIsolineNumRes(int m, int n, int res)
    {
	uiso_ = m;
	viso_ = n;
	isores_ = res;
	tesselateIsolines();
    }
    /// Fetch the number of iso parametric curves to compute
    void getIsolineNumRes(int& m, int& n, int& res)
    {
	n = uiso_;
	m = viso_;
	res = isores_;
    }

private:
    void tesselateSurface();
    void tesselateIsolines();

    const ParamSurface& surf_;
    shared_ptr<RegularMesh> mesh_;
    std::vector<LineStrip> isolinestrips_;
    bool isolines_;
    int uiso_;
    int viso_;
    int isores_;
};

} // namespace Go


#endif // RECTANGULARSURFACETESSELATOR_H

