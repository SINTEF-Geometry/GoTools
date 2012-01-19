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
/** PointOnEdge - represents a point lying on an edge. Storage.
 * 
 */
//===========================================================================
class PointOnEdge
{

 public:
  /// Constructor
  /// \param the edge on which the point lies
  /// \param par associated edge parameter
  PointOnEdge(ftEdge* edge, double par);

  /// Destructor
  ~PointOnEdge();
	
  /// Point position
  const Point& position() const { return pt_; }
  /// Associated edge
  ftEdge* edge() const { return edge_; }
  /// Associated edge parameter
  double par() const { return par_; }
    
 private:
  ftEdge* edge_;
  double par_;
  Point pt_;
};

} // namespace Go


#endif // _POINTONEDGE_H

