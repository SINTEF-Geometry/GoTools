//===========================================================================
//
// File : FaceUtilites.h
//
// Created: December 2011
//
// Author: Vibeke Skytt
//
// Revision: 
//
// Description:
//
//===========================================================================


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
