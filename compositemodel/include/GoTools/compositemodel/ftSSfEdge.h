//===========================================================================
//                                                                           
// File: ftSSfEdge.h                                                            
//                                                                           
// Created: 010802.
//                                                                           
// Author: Vibeke Skytt, SINTEF.
//                                                                           
// Revision: 
//                                                                           
// Description:
//                                                                           
//===========================================================================

#ifndef _FTSSFEDGE_H
#define _FTSSFEDGE_H


#include "GoTools/compositemodel/ftEdge.h"
#include "GoTools/compositemodel/ftFaceBase.h"
//#include "GoTools/model_toolbox/ftSuperSurface.h"

namespace Go
{


  class ftEdge;

  //===========================================================================
  /** ftSSfEdge - topological edge for Fantastic corresponding to ftSuperSurface
   * Detailed description.
   *
   * The ftSSfEdge is a half-edge implementation of a topological data structure.
   * It implements the ftEdgeBase interface by using the Go geometry library.
   *
   * \author Atgeirr F Rasmussen <atgeirr@sintef.no>
   * \bug Not tested
   * \see ftEdgeBase
   */
  //===========================================================================
  class ftSSfEdge : public ftEdgeBase
  {
  public:

    /** Constructor.
     * Detailed description.
     */
    ftSSfEdge(ftFaceBase* face, ftEdge *edge, int entry_id = -1);
    /// Empty destructor
    ~ftSSfEdge();

    virtual double tMin() const
    {
      return edg_->tMin();
    }
    virtual double tMax() const
    {
      return edg_->tMax();
    }
/*     virtual void turnOrientation(); */
/*     virtual void setOrientation(); */
/*     virtual bool isTurned(); */

    virtual void setReversed(bool is_reversed)
    {
      THROW("Not implemented!");
    }

    virtual bool isReversed()
    {
      THROW("Not implemented!");
    }

    virtual void setFace(ftFaceBase* face)
    {
      THROW("Not implemented!");
    }

    virtual ftFaceBase* face();

    // The bounding box is not exact, it is much too large...
    // (in some cases). It is implemented as the bounding box
    // of the WHOLE edgecurve instead of only the piece covered
    // by this halfedge.
    virtual BoundingBox boundingBox();

#if ((_MSC_VER > 0) && (_MSC_VER < 1300))
    virtual ftEdgeBase* split(double t);
#else
    virtual ftSSfEdge* split(double t);
#endif

    virtual int entryId() { return entry_id_; }
    virtual void setEntryId(int id) { entry_id_ = id; }

    virtual Point point(double t) const;
    virtual Point tangent(double t) const;
    virtual Point normal(double t) const;
    virtual Point normal(double t, Point& face_par_pt, double* face_seed) const;
    virtual void closestPoint(const Point& pt, double& clo_t,
			      Point& clo_pt, double& clo_dist,
			      double const *seed = 0) const;

    /// Return edge pointer
    virtual ftEdge* geomEdge(); 


  private:
    ftFaceBase* face_;
    ftEdge* edg_;

    int entry_id_;
  };

} // namspace Go


#endif // _FTSSFEDGE_H

