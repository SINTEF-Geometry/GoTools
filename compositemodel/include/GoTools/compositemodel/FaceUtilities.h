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

#ifndef __FACEUTILITIES_H
#define __FACEUTILITIES_H

#include "GoTools/utils/Point.h"
#include "GoTools/compositemodel/ftSurface.h"
#include "GoTools/compositemodel/ftEdge.h"

namespace Go
{
  /// \brief Sample data related to one face. 
  /// The struct contains one point in a point set with information 
  /// about position, associated surface normal and curvature.
  /// Points lying on the face boundary know about its assoicated edge. If the
  /// surface normal and curvature information regarding boundary points is not
  /// unique, the assosiated data is set to be equal to MAX_DOUBLE
  struct SamplePointData
  {
    Point pos_;
    Point norm_;
    double mean_curvature_;
    ftSurface *face_;
    double face_par_[2];
    ftEdge *edge_;
    double edge_par_;

    /// Constructor for points at the face boundary
    SamplePointData(Point pos, Point norm, double curvature, 
		    ftSurface* face, double face_par_u, double face_par_v, 
		    ftEdge *edge, double edge_par)
    {
      pos_ = pos;
      norm_ = norm;
      mean_curvature_ = curvature;
      face_ = face;
      face_par_[0] = face_par_u;
      face_par_[1] = face_par_v;
      edge_ = edge;
      edge_par_ = edge_par;
    } 

    /// Constructor for points in the inner of the face
    SamplePointData(Point pos, Point norm, double curvature, 
		    ftSurface* face, double face_par_u, double face_par_v)
    {
      pos_ = pos;
      norm_ = norm;
      mean_curvature_ = curvature;
      face_ = face;
      face_par_[0] = face_par_u;
      face_par_[1] = face_par_v;
      edge_ = NULL;
      edge_par_ = -1.0;
    }

  };

  /// Utility functionality for faces in order to release the need for
  /// private functions in ftSurface
  namespace FaceUtilities
  {
    void getBoundaryData(ftSurface* face, int nmb_sample, 
			 std::vector<SamplePointData>& sample_points);

    void getInnerData(ftSurface* face, int nmb_sample_u, int nmb_sample_v, 
		      std::vector<SamplePointData>& sample_points);

    /// Enforce colinearity of coefficients of spline surfaces
    /// related to the given faces
    /// Return value = false : No modification
    bool enforceCoLinearity(ftSurface *face1, ftEdge *edge1,
			    ftSurface *face2,
			    double tol, double ang_tol);

    bool enforceVxCoLinearity(shared_ptr<Vertex> vx, 
			      double tol, double ang_tol);
  };
}

#endif    // #ifndef __FACEUTILITIES_H
