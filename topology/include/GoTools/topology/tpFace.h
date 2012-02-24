//===========================================================================
//                                                                           
// File: tpFace.h                                                            
//                                                                           
// Created: Tue Mar 21 15:23:59 2000                                         
//                                                                           
// Author: Atgeirr F Rasmussen <atgeirr@sintef.no>
//                                                                           
// Revision: $Id: tpFace.h,v 1.27 2005-06-09 07:27:03 oan Exp $
//                                                                           
// Description:
//                                                                           
//===========================================================================

#ifndef _TOPFACE_H
#define _TOPFACE_H

#include "GoTools/utils/Values.h"
#include "GoTools/utils/BoundingBox.h"
#include <memory>
#include <vector>

namespace Go
{

class tpEdge;


//===========================================================================
/** tpFace -  
 * Minimal structure of template faceType when using FaceAdjacency.
 * The interface is limited to the functionality required to compute
 * adjacency information in a set of faces and to remove this face
 * from a face set, i.e. remove all topology pointers related to this
 * face.
 *
 * \author Atgeirr F Rasmussen <atgeirr@sintef.no>
 * \see Class1
 */
//===========================================================================

class tpFace
{
public:
    // These first members functions are needed by template class faceType when
    // using the tpTopologyTable.

    /// Destructor
    virtual ~tpFace();
    /// Compute the edges associated to this face or fetch already existing
    /// edges
    virtual std::vector<shared_ptr<tpEdge> > 
      createInitialEdges(double degenerate_epsilon=DEFAULT_SPACE_EPSILON) = 0;
    /// Return pointers to first part of all bd cvs.
    virtual std::vector<shared_ptr<tpEdge> > startEdges() = 0;
    /// Evaluate point on face
    virtual Point point(double u, double v) const = 0;
    /// Evaluate surface normal
    virtual Point normal(double u, double v) const = 0;
    /// The bounding box corresponding to this face
    virtual BoundingBox boundingBox() = 0;
    /// Return id, default id is -1
    virtual int getId() = 0;
    /// Remove all adjacency information related to this face
    virtual void isolateFace();
    //virtual std::vector<shared_ptr<tpEdge> > 
    //setOrientation(double degenerate_epsilon=DEFAULT_SPACE_EPSILON) = 0;
    //void turnFace(std::vector<tpFace*>& turned);

    // The following member functions are not needed for the tpTopologyTable
    // to work, but are natural in this setting.
    //virtual void turnOrientation() = 0; // Called from our version of turnFace.

};

} // namespace Go

#endif // _TOPFACE_H

