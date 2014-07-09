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
    void resetPosition(Vector3D pos);

    /// Write parameter information to stream
    virtual
    void write2Dval(std::ostream& os) const;


  private:

    // Each ftSamplePoint may be a member of several patches (if on boundary).
    std::vector<std::pair<shared_ptr<ftFaceBase>, Vector2D> > par_pts_;

  };


} // namespace Go

#endif // _FTSURFACESETPOINT_H

