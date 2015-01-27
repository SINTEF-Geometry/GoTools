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

#ifndef _REGULARIZEFACE_H
#define _REGULARIZEFACE_H

#include "GoTools/compositemodel/ftSurface.h"
#include "GoTools/geometry/BoundedSurface.h"

namespace Go {

/// \brief Split one face into a number of 4-sided domains without inner trimming.
/// This class is intended for use in block structuring. One face with possible 
/// inner and outer trimming is split according to certain rules to result in a
/// face set with 4 sided faces although faces with less than 4 sides can occur.
/// A side is defined as a piece of the face boundary between two corners or between
/// vertices where there are more than one adjacent face.
/// The trimmed surfaces being output from this class can later be approximated by
/// spline surfaces.
/// The splitting procedure is recursive.

class RegularizeFace
{
 public:
  /// Constructor
  RegularizeFace(shared_ptr<ftSurface> face, 
		 double epsge, double angtol, double tol2,
		 bool split_in_cand=false);

  /// Constructor
  RegularizeFace(shared_ptr<ftSurface> face, 
		 double epsge, double angtol, double tol2, double bend,
		 bool split_in_cand=false);

  /// Constructor
  RegularizeFace(shared_ptr<ftSurface> face, 
		 shared_ptr<SurfaceModel> model,
		 bool split_in_cand=false);

    /// Destructor
  ~RegularizeFace();

  /// Set information about the centre of the current face to the regularization
  /// of sub faces
  void setAxis(Point& centre, Point& axis);

  void unsetAxis();

  /// Set info about splitting performed in opposite faces in a body.
  /// Used from RegularizeFaceSet.
  void setCandSplit(std::vector<std::pair<std::pair<Point,int>,
		    std::pair<Point,int> > >  cand_split)
  {
    cand_split_ = cand_split;
  }

  /// Set information about faces not to be the cause of T-joint splitting
  void setNonTjointFaces(std::vector<shared_ptr<ftSurface> >& faces)
  {
    nonTjoint_faces_ = faces;
  }

  void setSplitMode(int split_mode)
  {
    split_mode_ = split_mode;
  }

  /// Classify vertices according to significance. Mark vertices that should
  /// not trigger splitting
  void classifyVertices();

  /// Fetch result
  std::vector<shared_ptr<ftSurface> > getRegularFaces();

  /// Fetch info about removed seams
  std::vector<Point> getSeamJointInfo() const
    {
      return seam_joints_;
    }

  /// Decides whether T-joint splitting should be performed at face level
  void setDivideInT(bool divideInT)
  {
    divideInT_ = divideInT;
  }

  /// Fetch info about point corrspondance
  std::vector<std::pair<Point, Point> > fetchVxPntCorr()
    {
      return corr_vx_pts_;
    }

  private:

  /// Struct to store face hole information. Related to face regularization.
  struct hole_info
  {
    /// Estimated hole centre
    Point hole_centre_;
    /// Estimated hole axis
    Point hole_axis_;
    /// Estimated hole radius
    double hole_radius_;

    void setInfo(Point& centre, Point& axis, double radius)
    {
      hole_centre_ = centre;
      hole_axis_ = axis;
      hole_radius_ = radius;
    }
  };

  double epsge_;  // Geometry tolerance
  double angtol_; // Angular tolerance
  double tol2_;   // Neighbourhood tolerance
  double bend_;  // Corner tolerance

  shared_ptr<ftSurface> face_;  // Surface corresponding to face
  std::vector<shared_ptr<ftSurface> > sub_faces_;  // Current
  // division of face
  shared_ptr<SurfaceModel> model_;

  Point centre_;  // Weightpoint of hole centra
  Point axis_;    // Normal axis corresponding to weight point
  double radius_;

  int split_mode_;
  int divideInT_;
  bool top_level_;
  double isolate_fac_;

  std::vector<shared_ptr<Vertex> > vx_;
  std::vector<shared_ptr<Vertex> > corners_;

  bool split_in_cand_;
  std::vector<std::pair<std::pair<Point,int>, std::pair<Point,int> > >  cand_split_;

  std::vector<shared_ptr<Vertex> > non_sign_vx_;
  std::vector<shared_ptr<Vertex> > seam_vx_;

  std::vector<Point> seam_joints_;
  
  std::vector<std::pair<Point,Point> > corr_vx_pts_;

  std::vector<shared_ptr<ftSurface> > nonTjoint_faces_;

    // Perform division
  void divide();

  void splitInTJoints();

  std::vector<shared_ptr<ftSurface> > 
    divideInTjoint(shared_ptr<ftSurface> face,
		   std::vector<shared_ptr<Vertex> >& Tvx,
		   std::vector<shared_ptr<Vertex> >& corner);

void faceWithHoles(std::vector<std::vector<ftEdge*> >& half_holes);

  void faceOneHole(std::vector<std::vector<ftEdge*> >& half_holes);

  void faceOneHole2();

  shared_ptr<CurveOnSurface> 
    computeCornerSplit(shared_ptr<Vertex> corner,
		       std::vector<shared_ptr<Vertex> >& hole_vx,
		       std::vector<shared_ptr<Vertex> >& hole_vx2,
		       shared_ptr<BoundedSurface>& bd_sf,
		       const Point& close, bool outer_vx=true);

  std::vector<shared_ptr<ftSurface> >
    faceOuterBdFaces(std::vector<std::vector<ftEdge*> >& half_holes);

  void faceOuterBd(std::vector<std::vector<ftEdge*> >& half_holes);


  std::vector<shared_ptr<ftSurface> > 
    divideVertex(ftEdge* edge, double par,
		 std::vector<shared_ptr<Vertex> > cand_vx,
		 ftEdge* cand_edge);


  double getSegmentAngle(shared_ptr<Vertex> vx1,
			 shared_ptr<Vertex> vx2,
			 Point& pnt, Point& normal);

  bool getVertexProperties(shared_ptr<Vertex> vx, Point& parpnt,
			   double& ang, bool& T_joint);

 shared_ptr<Vertex> 
    getSignificantVertex(std::vector<shared_ptr<Vertex> > cand_vx);

 std::vector<shared_ptr<Vertex> >
   prioritizeCornerVx(std::vector<shared_ptr<Vertex> > cand_vx);

  std::vector<std::vector<ftEdge*> > getHalfHoles(int idx=0);

  std::vector<shared_ptr<ftSurface> > 
    isolateHole(const Point& seg_pnt, shared_ptr<Vertex> vx,
		std::vector<std::vector<ftEdge*> >& half_holes);

  std::vector<shared_ptr<ftSurface> > 
    initIsolateHole(std::vector<std::vector<ftEdge*> >& half_holes);
  
  void 
    selectCandidateSplit(shared_ptr<Vertex> select_vx,
			 std::vector<shared_ptr<Vertex> >& vx,
			 std::vector<shared_ptr<Vertex> >& cand_vx,
			 ftEdge*& cand_edge, bool keep_T_joints=true);

  void 
    selectCandidateSplit(ftEdge* edge,
			 std::vector<shared_ptr<Vertex> >& vx,
			 std::vector<shared_ptr<Vertex> >& cand_vx,
			 ftEdge*& cand_edge);

  bool
    sortAlongLine(std::vector<hole_info>& holes, Point& pnt,
		  Point& dir, std::vector<double>& parvals,
		  std::vector<int>& perm);

  void
    identifyParLine(std::vector<hole_info>& holes, Point& dir);

  bool
    sortRadially(std::vector<hole_info>& holes, const Point& wgt_pnt,
		 const Point& norm, std::vector<double>& angles, 
		 std::vector<int>& perm);

  std::vector<shared_ptr<ftSurface> >
    divideAcrossLine(std::vector<std::vector<ftEdge*> >& half_holes,
		     std::vector<hole_info>& holes, Point& pnt,
		     Point& dir, std::vector<int>& perm);

  std::vector<shared_ptr<ftSurface> >
    isolateHolesParametrically(std::vector<std::vector<ftEdge*> >& half_holes,
			       Point& dir);

  std::vector<shared_ptr<ftSurface> >
    isolateHolesRadially(std::vector<std::vector<ftEdge*> >& half_holes,
			 const Point& mid, const Point& axis,
			 int loop_idx,
			 std::vector<hole_info>& holes, 
			 std::vector<int>& perm);
  std::vector<shared_ptr<ftSurface> >
    isolateHolesRadially2(std::vector<std::vector<ftEdge*> >& half_holes,
			  const Point& mid, const Point& axis,
			  int loop_idx,
			  std::vector<hole_info>& holes, 
			  std::vector<int>& perm);

  std::vector<shared_ptr<ftSurface> >
    isolateOneHoleRadially(const Point& mid, const Point& axis,
			   hole_info& hole);

  std::vector<shared_ptr<ftSurface> >
    divideByPlanes(std::vector<Point>& pnts, std::vector<Point>& normals,
		   std::vector<std::vector<ftEdge*> >& half_holes,
		   const vector<double>& level_dist);

  std::vector<shared_ptr<ftSurface> >
    divideByPlanes(std::vector<Point>& pnts, 
		   const Point& mid, const Point& axis,
		   int loop_idx,
		   std::vector<std::vector<ftEdge*> >& half_holes,
		   std::vector<hole_info>& holes, 
		   double level_dist);

  std::vector<shared_ptr<ftSurface> >
    holeToHoleSplit(std::vector<vector<ftEdge*> >& half_holes,
		    std::vector<hole_info>& holes, 
		    std::vector<std::pair<int,int> >& hole_idx,
		    std::vector<double>& seg_lengts,
		    std::vector<std::pair<Point,Point> >& seg_endpt,
		    int loop_idx);

  void extractCandPt(Point mid, int hole_ix,
		     std::vector<shared_ptr<CurveOnSurface> >& seg,
		     Point cand_pt[], Point cand_par[]);

  bool 
    adjustTrimSeg(shared_ptr<CurveOnSurface>& trim_seg,
		  shared_ptr<Vertex> vx1,
		  const std::vector<ftEdge*>& edges1,
		  shared_ptr<Vertex> vx2,
		  const std::vector<ftEdge*>& edges2,
		  double len);

  int
    positionWeigthPoint(const Point& wgt_par);

  void removeHalfHoleVx(std::vector<shared_ptr<Vertex> >& vx,
			std::vector<std::vector<ftEdge*> >& half_holes);

  ftSurface* 
    identifySeamFaces(shared_ptr<ftSurface> face1,
		      std::vector<shared_ptr<ftSurface> > faces,
		      int& pardir);

  shared_ptr<ftSurface>
    mergeSeamFaces(ftSurface* face1, ftSurface* face2, int pardir);

  bool
    fetchPatternSplit(const Point& corner,
		      Point& parval1, Point& parval2,
		      bool use_input_point = true);

  bool
    fetchPatternSplit2(const Point& cant_pt,
		       Point& parval1, Point& parval2,
		       double level_dist, double level_ang,
		       const Point& normal, bool only_outer_bd);

  int nmbSplitPattern(const Point& p1, const Point& p2);  

  int fetchSplitPattern(const Point& p1, const Point& p2,
			std::vector<shared_ptr<Vertex> >& vx,
			std::vector<int>& vx_idx, int nmb_idx,
			std::vector<std::pair<Point,Point> >& pattern,
			std::vector<std::pair<int,int> >& pattern_vx);

  bool splitWithPatternLoop();

  void removeOuterCands(std::vector<shared_ptr<CurveOnSurface> >& cand_cvs);

  void
    removeInsignificantVertices(std::vector<shared_ptr<Vertex> >& vx,
				bool keep_T_joints = false,
				ftSurface *face=NULL);

  void mergeSeams(std::vector<shared_ptr<ftSurface> >& faces, int& nmb_faces,
		  std::vector<shared_ptr<ftSurface> >& faces2);
  void 
    splitTrimSegments(std::vector<shared_ptr<CurveOnSurface> >& segments);

  void getConcaveCorners(std::vector<shared_ptr<Vertex> >& corners, 
			 std::vector<shared_ptr<Vertex> >& concave_corners);

  std::vector<shared_ptr<ftSurface> >
    chopOffRegBlocks(std::vector<shared_ptr<Vertex> >& concave_corners);

  std::vector<shared_ptr<ftSurface> >
    connectToVertex(std::vector<shared_ptr<Vertex> >& concave_corners);

  bool checkRegularity(std::vector<shared_ptr<Vertex> >& cand_vx);

  /// Unset top level info
  void unsetTopLevel()
  {
    top_level_ = false;
  }

  bool checkTrimSegments(std::vector<shared_ptr<CurveOnSurface> >& trim_segments,
			 shared_ptr<Vertex> corner, Point pnt,
			 shared_ptr<BoundedSurface>& bd_sf,
			 bool outer_vx);

  void getSegPntRadialSplit(hole_info& hole1, hole_info& hole2,
			    double len_fac,
			    std::vector<Point>& seg_pnt,
			    std::vector<Point>& seg_norm,
			    std::vector<double>& min_dist);

  void snapToVertex(Point& pos, Point& par, 
		    vector<shared_ptr<Vertex> >& vxs,
		    double lim);

  void updateTrimSeg(std::vector<shared_ptr<CurveOnSurface> >& trim_segments,
		     const Point& pnt1, const Point& param1, bool at_vx1, 
		     int loop_idx1, const Point& pnt2, const Point& param2, 
		     bool at_vx2, int loop_idx2);

};

}  // namespace Go

#endif // _REGULARIZEFACE_H
