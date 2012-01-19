//===========================================================================
//                                                                           
// File: CurveTesselator.h                                                 
//                                                                           
// Created: Wed Nov 28 16:38:43 2001                                         
//                                                                           
// Author: Atgeirr F Rasmussen <atgeirr@sintef.no>
//                                                                           
// Revision: $Id: CurveTesselator.h,v 1.1 2008-03-27 16:13:09 vsk Exp $
//                                                                           
// Description:
//                                                                           
//===========================================================================

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

