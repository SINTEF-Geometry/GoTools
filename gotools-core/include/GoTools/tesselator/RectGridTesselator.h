//===========================================================================
//                                                                           
// File: RectGridTesselator.h                                              
//                                                                           
// Created: Wed Jan 19 13:14:14 2005                                         
//                                                                           
// Author: Atgeirr F Rasmussen <atgeirr@sintef.no>
//                                                                           
// Revision: 
//                                                                           
//===========================================================================

#ifndef _RECTGRIDTESSELATOR_H
#define _RECTGRIDTESSELATOR_H


#include "GoTools/tesselator/Tesselator.h"
#include "GoTools/geometry/RectGrid.h"
#include "GoTools/tesselator/QuadMesh.h"
#include <memory>

namespace Go
{

    /** Brief description. 
     *  Detailed description.
     */

class GO_API RectGridTesselator : public Tesselator
{
public:
  RectGridTesselator(const RectGrid& rg);
    
    virtual ~RectGridTesselator();
  
    virtual void tesselate();

    std::shared_ptr<QuadMesh> getMesh()
    {
	return quadmesh_;
    }

private:
    const RectGrid& rectgrid_;
    std::shared_ptr<QuadMesh> quadmesh_;
};

} // namespace Go


#endif // _RECTGRIDTESSELATOR_H

