//===========================================================================
//                                                                           
// File: Path.h                                                       
//                                                                           
// Created: April 2010
//                                                                           
// Author: Vibeke Skytt
//                                                                           
// Revision: 
//                                                                           
// Description:
//                                                                           
//===========================================================================

#ifndef _PATH_H
#define _PATH_H

#include "GoTools/compositemodel/Vertex.h"
#include "GoTools/compositemodel/ftEdge.h"
#include <vector>


namespace Go
{
  /// Functions related to sequence of edges
  namespace Path
  {
    /// Estimate mid point, normal and radius defined by an edge
    /// sequence
    bool estimateHoleInfo(std::vector<ftEdge*> edges, Point& centre, 
			  Point& axis, double& radius);

    /// Identify a loops starting and ending in a given vertex in an ordered
    /// sequence of edges
    std::vector<ftEdge*> identifyLoop(std::vector<ftEdge*> edges, 
				      std::shared_ptr<Vertex> vx);

  }  // namespace Patch

}  // namespace Go


#endif // _PATH_H
