//==========================================================================
//                                                                          
// File: RegularizeFace.h
//                                                                          
// Created: April 2010
//                                                                          
// Author: Vibeke Skytt
//                                                                          
// Revision: 
//                                                                          
// Description: Split one face into a number of 4-sided domains without
//              inner trimming
//                                                                          
//==========================================================================

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

   /// Destructor
  ~RegularizeFace();

  /// Set information about the centre of the current face to the regularization
  /// of sub faces
  void setAxis(Point& centre, Point& axis);

  void unsetAxis();

  /// Set info about splitting performed in opposite faces in a body.
  /// Used from RegularizeFaceSet.
  void setCandSplit(std::vector<std::pair<Point,Point> >  cand_split)
  {
    cand_split_ = cand_split;
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

  Point centre_;  // Weightpoint of hole centra
  Point axis_;    // Normal axis corresponding to weight point
  double radius_;

  bool divideInT_;

  std::vector<shared_ptr<Vertex> > vx_;
  std::vector<shared_ptr<Vertex> > corners_;

  bool split_in_cand_;
  std::vector<std::pair<Point,Point> >  cand_split_;

  std::vector<shared_ptr<Vertex> > non_sign_vx_;
  std::vector<shared_ptr<Vertex> > seam_vx_;

  std::vector<Point> seam_joints_;
  
  std::vector<std::pair<Point,Point> > corr_vx_pts_;

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
		       bool outer_vx=true);

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
			 ftEdge*& cand_edge);

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
		   double level_dist);

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
    fetchPatternSplit(Point& corner,
		      Point& parval1, Point& parval2,
		      bool use_input_point = true);

  int nmbSplitPattern(const Point& p1, const Point& p2);  

  int fetchSplitPattern(const Point& p1, const Point& p2,
			std::vector<shared_ptr<Vertex> >& vx,
			std::vector<int>& vx_idx, int nmb_idx,
			std::vector<std::pair<Point,Point> >& pattern,
			std::vector<std::pair<int,int> >& pattern_vx);

  void splitWithPatternLoop();

  void removeOuterCands(std::vector<shared_ptr<CurveOnSurface> >& cand_cvs);

  void
    removeInsignificantVertices(std::vector<shared_ptr<Vertex> >& vx);

  void mergeSeams(std::vector<shared_ptr<ftSurface> >& faces, int& nmb_faces,
		  std::vector<shared_ptr<ftSurface> >& faces2);
  void 
    splitTrimSegments(std::vector<shared_ptr<CurveOnSurface> >& segments);
};

}  // namespace Go

#endif // _REGULARIZEFACE_H
