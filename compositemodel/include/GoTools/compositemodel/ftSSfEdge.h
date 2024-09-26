/*
 * Copyright (C) 1998, 2000-2007, 2010, 2011, 2012, 2013 SINTEF ICT,
 * Applied Mathematics, Norway.
 *
 * Contact information: E-mail: tor.dokken@sintef.no                      
 * SINTEF ICT, Department of Applied Mathematics,                         
 * P.O. Box 124 Blindern,                                                 
 * 0314 Oslo, Norway.                                                     
 *
 * This file is part of GoTools.
 *
 * GoTools is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Affero General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version. 
 *
 * GoTools is distributed in the hope that it will be useful,        
 * but WITHOUT ANY WARRANTY; without even the implied warranty of         
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the          
 * GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public
 * License along with GoTools. If not, see
 * <http://www.gnu.org/licenses/>.
 *
 * In accordance with Section 7(b) of the GNU Affero General Public
 * License, a covered work must retain the producer line in every data
 * file that is created or manipulated using GoTools.
 *
 * Other Usage
 * You can be released from the requirements of the license by purchasing
 * a commercial license. Buying such a license is mandatory as soon as you
 * develop commercial activities involving the GoTools library without
 * disclosing the source code of your own applications.
 *
 * This file may be used in accordance with the terms contained in a
 * written agreement between you and SINTEF ICT. 
 */

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

