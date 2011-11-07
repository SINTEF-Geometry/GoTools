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

using namespace Go;

class tpEdge;


//===========================================================================
/** tpFace -  Short description.
 * Detailed description.
 * Minimal structure of template faceType when using tpTopologyTable.
 *
 * \author Atgeirr F Rasmussen <atgeirr@sintef.no>
 * \bug It isn't documented yet!
 * \see Class1
 */
//===========================================================================

class tpFace
{
public:
    // These first members functions are needed by template class faceType when
    // using the tpTopologyTable.

    virtual ~tpFace();
    virtual std::vector<std::shared_ptr<tpEdge> > 
      createInitialEdges(double degenerate_epsilon=DEFAULT_SPACE_EPSILON) = 0;
    virtual std::vector<std::shared_ptr<tpEdge> > startEdges() = 0;
    virtual Point point(double u, double v) const = 0;
    virtual Point normal(double u, double v) const = 0;
    virtual BoundingBox boundingBox() = 0;
    virtual int getId() = 0;
    virtual void isolateFace();
    //virtual std::vector<std::shared_ptr<tpEdge> > 
    //setOrientation(double degenerate_epsilon=DEFAULT_SPACE_EPSILON) = 0;
    //void turnFace(std::vector<tpFace*>& turned);

    // The following member functions are not needed for the tpTopologyTable
    // to work, but are natural in this setting.
    //virtual void turnOrientation() = 0; // Called from our version of turnFace.

};


#endif // _TOPFACE_H

