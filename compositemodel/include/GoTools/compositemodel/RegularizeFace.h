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

/// Split one face into a number of 4-sided domains without inner trimming
class RegularizeFace
{
 public:
  /// Constructor
  RegularizeFace(shared_ptr<ftSurface> face, 
		 double epsge, double angtol, double tol2);

  /// Constructor
  RegularizeFace(shared_ptr<ftSurface> face, 
		 double epsge, double angtol, double tol2, double bend);

   /// Destructor
  ~RegularizeFace();

  void setAxis(Point& centre, Point& axis);

  void setCandParams(std::vector<std::pair<Point,Point> >  cand_params)
  {
    cand_params_ = cand_params;
  }

  void classifyVertices();

  /// Fetch result
  std::vector<shared_ptr<ftSurface> > getRegularFaces();

  /// Fetch info about removed seams
  std::vector<Point> getSeamJointInfo() const
    {
      return seam_joints_;
    }

  void setDivideInT(bool divideInT)
  {
    divideInT_ = divideInT;
  }

  private:

  /// Struct to store face hole information.
  struct hole_info
  {
    Point hole_centre_;
    Point hole_axis_;
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

  std::vector<std::pair<Point,Point> >  cand_params_;

  std::vector<shared_ptr<Vertex> > non_sign_vx_;
  std::vector<shared_ptr<Vertex> > seam_vx_;

  std::vector<Point> seam_joints_;
  

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
		       shared_ptr<BoundedSurface>& bd_sf);

  void faceOuterBd(std::vector<std::vector<ftEdge*> >& half_holes);


  std::vector<shared_ptr<ftSurface> > 
    divideVertex(ftEdge* edge, double par,
		 std::vector<shared_ptr<Vertex> > cand_vx,
		 ftEdge* cand_edge);


  double getSegmentAngle(shared_ptr<Vertex> vx1,
			 shared_ptr<Vertex> vx2,
			 Point& pnt, Point& normal);

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
    fetchPatternSplit(shared_ptr<Vertex> corner,
		      Point& parval1, Point& parval2);
  void
    removeInsignificantVertices(std::vector<shared_ptr<Vertex> >& vx);

  void mergeSeams(std::vector<shared_ptr<ftSurface> >& faces, int& nmb_faces,
		  std::vector<shared_ptr<ftSurface> >& faces2);
  void 
    splitTrimSegments(std::vector<shared_ptr<CurveOnSurface> >& segments);
};

}  // namespace Go

#endif // _REGULARIZEFACE_H
