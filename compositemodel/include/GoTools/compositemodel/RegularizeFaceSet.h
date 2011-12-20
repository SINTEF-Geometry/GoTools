//==========================================================================
//                                                                          
// File: RegularizeFaceSet.h
//                                                                          
// Created: August 2010
//                                                                          
// Author: Vibeke Skytt
//                                                                          
// Revision: 
//                                                                          
// Description: Split a set of faces into a number of 4-sided domains without
//              inner trimming
//                                                                          
//==========================================================================

#ifndef _REGULARIZEFACESET_H
#define _REGULARIZEFACESET_H

#include "GoTools/compositemodel/ftSurface.h"
#include "GoTools/geometry/BoundedSurface.h"
#include "GoTools/geometry/BoundedUtils.h"

namespace Go {

  class ftEdge;
  class Vertex;
  class EdgeVertex;

/// Split a set of faces into a number of 4-sided domains without inner trimming.
class RegularizeFaceSet
{
 public:
  /// Constructor
  RegularizeFaceSet(std::vector<shared_ptr<ftSurface> > faces, 
		 double epsge, double angtol);
  /// Constructor
  RegularizeFaceSet(shared_ptr<SurfaceModel> model);
  /// Destructor
  ~RegularizeFaceSet();

  /// Set information
  void setFaceCorrespondance(int idx1, int idx2);

  /// Fetch result
  std::vector<shared_ptr<ftSurface> > getRegularFaces();

  shared_ptr<SurfaceModel> getRegularModel();

  private:
  shared_ptr<SurfaceModel> model_;

  std::vector<std::pair<int,int> > corr_faces_;

  std::vector<std::vector<std::pair<Point,Point> > > cand_params_;

  std::vector<Point> seam_joints_;

    // Perform division
  void divide();

  void splitInTJoints();

  std::vector<shared_ptr<ftSurface> > 
    divideInTjoint(shared_ptr<ftSurface>& face,
		   std::vector<shared_ptr<Vertex> >& Tvx,
		   std::vector<shared_ptr<Vertex> >& corner);


  void 
    selectCandidateSplit(shared_ptr<ftSurface> face,
			 shared_ptr<Vertex> select_vx,
			 std::vector<shared_ptr<Vertex> >& vx,
			 std::vector<shared_ptr<Vertex> >& cand_vx,
			 ftEdge*& cand_edge);

  double getSegmentAngle(shared_ptr<ftSurface> face,
			 shared_ptr<Vertex> vx1,
			 shared_ptr<Vertex> vx2,
			 Point& pnt, Point& normal);


  std::vector<std::pair<Point,Point> > 
    getEndParameters(std::vector<shared_ptr<ftSurface> >& faces);

  ftSurface*
    identifySeamFaces(shared_ptr<ftSurface> face1, int& pardir,
		      int& status);

  int
    getParameterDirection(ftSurface* face1, ftSurface* face2,
			  shared_ptr<ftEdge> edge1,
			  shared_ptr<ftEdge> edge2,
			  double eps, int& pardir);
  void
    removeInsignificantVertices(std::vector<shared_ptr<Vertex> >& vx);

  void seamClassification();

  void
    getSeamRadialEdge(ftSurface* face, 
		      std::vector<shared_ptr<EdgeVertex> >& edgevx, 
		      std::vector<std::pair<Point,Point> >& endpts);

  void
    attachRadialEdge(ftSurface* face, 
		     std::vector<shared_ptr<EdgeVertex> >& edgevx, 
		     std::vector<std::pair<Point,Point> >& endpts,
		     double tol);

  bool mergeSituation(ftSurface* face, 
		      shared_ptr<Vertex> vx,
		      ftSurface*& merge1,
		      int& dir1, double& val1, bool& atstart1,
		      ftSurface*& merge2,
		      int& dir2, double& val2, bool& atstart2,
		      shared_ptr<Vertex>& other_vx,
		      std::pair<Point, Point>& co_par1,
		      std::pair<Point, Point>& co_par2);

void prioritizeFaces(std::vector<std::shared_ptr<ftSurface> >& faces,
		     std::vector<int>& perm);
};

}  // namespace Go

#endif // _REGULARIZEFACESET_H
