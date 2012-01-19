//===========================================================================
//                                                                           
// File: ftSurfaceSetPoint.h                                                 
//                                                                           
// Created: Fri Feb  1 13:38:59 2002                                         
//                                                                           
// Author: Sverre Briseid <Sverre.Briseid@math.sintef.no>
//                                                                           
// Revision: $Id: ftSurfaceSetPoint.h,v 1.1 2009-01-23 13:34:32 vsk Exp $
//                                                                           
// Description: Subclass of ftSamplePoint; intended for use when a point is
//              member of multiple ftFaceBase's.
//                                                                           
//===========================================================================

#ifndef _FTSURFACESETPOINT_H
#define _FTSURFACESETPOINT_H

#include "GoTools/compositemodel/ftPointSet.h"
#include "GoTools/compositemodel/ftFaceBase.h"


namespace Go
{


  /** Subclass of ftSamplePoint; intended for use when a point is
   *              member of multiple ftFaceBase's.
   */
  class ftSurfaceSetPoint : public ftSamplePoint
  {

  public:
    /// Constructor
    ftSurfaceSetPoint(Vector3D xyz, int bnd);
    /// Constructor
    ftSurfaceSetPoint(Vector3D xyz, int bnd,
		      shared_ptr<ftFaceBase>& face, Vector2D par_pt);
    /// Destructor
    virtual ~ftSurfaceSetPoint();

    /// Return pointer to this point as a surface set point
    virtual ftSurfaceSetPoint* asSurfaceSetPoint()
	{
	    return this;
	}

    virtual bool containsFace(ftFaceBase* face) const;

    /// For a given point, add information about param values on another ftFaceBase.
    void addPair( shared_ptr<ftFaceBase>& face, const Vector2D& par_pt);

    /// For a given point, add another face. Parameter values are computed
    void addFace(shared_ptr<ftFaceBase>& face);

    /// Return the number of faces in which point is member of.
    int nmbFaces();

    /// Return pointer to face number i (error if not enough faces).
    shared_ptr<ftFaceBase> face(int i);

    /// Fetch the parameter vale in face number i of this point
    Vector2D parValue(int i);

    /// Fetch the parameter value in face face of this point
    Vector2D getPar(ftFaceBase* face);

    //     ftSurfaceSetPoint* clone() const
    //     { return new ftSurfaceSetpoint(*this); }

    /// Add face, parameter and connectivity info from the point other to this point
    /// NB! It is assumed that the position information is consistent
    void addInfo(ftSurfaceSetPoint* other);

    /// Change position of point and update parameter values accordingly
    void resetPosition(Vector3D pos, int bnd);

    /// Write parameter information to stream
    virtual
    void write2Dval(std::ostream& os) const;


  private:

    // Each ftSamplePoint may be a member of several patches (if on boundary).
    std::vector<std::pair<shared_ptr<ftFaceBase>, Vector2D> > par_pts_;

  };


} // namespace Go

#endif // _FTSURFACESETPOINT_H

