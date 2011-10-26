//===========================================================================
//                                                                           
// File: PointOnEdge
//                                                                           
// Created: Dec. 08
//                                                                           
// Author: Vibeke Skytt
//                                                                           
// Revision: 
//                                                                           
// Description: 
//                                                                           
//===========================================================================

#ifndef _POINTONEDGE
#define _POINTONEDGE


#include "GoTools/utils/Point.h"

namespace Go
{

class ftEdge;

//===========================================================================
/** PointOnEdge - represents a point lying on an edge
 * 
 */
//===========================================================================
class PointOnEdge
    {

public:
	PointOnEdge(ftEdge* edge, double par);

	~PointOnEdge();

    const Point& position() const { return pt_; }
    ftEdge* edge() const { return edge_; }
    double par() const { return par_; }
    
    private:
    ftEdge* edge_;
    double par_;
    Point pt_;
};

} // namespace Go


#endif // _POINTONEDGE_H

