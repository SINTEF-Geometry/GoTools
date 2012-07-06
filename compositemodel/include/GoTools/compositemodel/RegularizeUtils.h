//==========================================================================
//                                                                          
// File: RegularizeUtils.h
//                                                                          
// Created: April 2011
//                                                                          
// Author: Vibeke Skytt
//                                                                          
// Revision: 
//                                                                          
// Description: Utility functionality for RegularizeFace and 
//              RegularizeFaceSet
//                                                                          
//==========================================================================

#ifndef _REGULARIZEUTILS_H
#define _REGULARIZEUTILS_H

#include "GoTools/compositemodel/ftSurface.h"
#include "GoTools/geometry/BoundedSurface.h"

namespace Go {

  /// Utility functionality for RegularizeFace and RegularizeFaceSet   
  namespace RegularizeUtils
  {
    std::vector<shared_ptr<ftSurface> > 
      divideVertex(shared_ptr<ftSurface> face,
		   shared_ptr<Vertex> vx, 
		   std::vector<shared_ptr<Vertex> > cand_vx,
		   ftEdge* cand_edge,
		   double epsge, double tol2, double angtol,
		   double bend,
		   std::vector<shared_ptr<Vertex> > non_corner,
		   const Point& centre, const Point& axis);

    std::vector<shared_ptr<ftSurface> > 
      createFaces(std::vector<shared_ptr<BoundedSurface> >& sub_sfs,
		  shared_ptr<ftSurface>  face,
		  double epsge, double tol2, double angtol,
		  std::vector<shared_ptr<Vertex> > non_corner);

    void getDivisionPlane(shared_ptr<ftSurface> face,
			  shared_ptr<Vertex> vx,
			  double epsge,
			  Point& pnt,
			  Point& normal);

    void 
      getClosestBoundaryPar(shared_ptr<ftSurface> face,
			    shared_ptr<Vertex> vx,
			    std::vector<shared_ptr<ParamCurve> >& vx_cvs,
			    const Point& pnt,
			    double epsge,
			    int& close_idx, double& close_dist,
			    Point& close_par);
    bool
      cornerInShortestPath(shared_ptr<Vertex> vx1,
			   shared_ptr<Vertex> vx2,
			   shared_ptr<ftSurface> face,
			   double angtol);

    bool
      getPath(ftEdge* edg, shared_ptr<Vertex> vx,
	      shared_ptr<Vertex> last, shared_ptr<ftSurface> face,
	      std::vector<ftEdge*>& path);

    bool
      noExtension(shared_ptr<Vertex> vx, ftSurface* face,
		  shared_ptr<Vertex>& vx2, std::pair<Point, Point>& co_par1, 
		  std::pair<Point, Point>& co_par2, int& dir1, int& dir2,
		  double& val1, double& val2, double angtol, bool check_constant_curve);

    double getMaxParFrac(shared_ptr<ftSurface> face);

    int selectCandVx(shared_ptr<ftSurface> face,
		     shared_ptr<Vertex> vx, const Point& in_vec, 
		     vector<shared_ptr<Vertex> > cand_vx,
		     RectDomain& dom,
		     double epsge, double angtol,
		     const Point& centre, const Point& normal,
		     std::vector<shared_ptr<ParamCurve> >& vx_cvs,
		     double close_dist, const Point& close_pt,
		     double& cyl_rad);

    void checkTrimSeg(std::vector<shared_ptr<CurveOnSurface> >& trim_segments,
		      std::vector<shared_ptr<Vertex> >& next_vxs,
		      const Point& vx_point, double epsge);

    ftEdge* getOppositeBoundaryPar(shared_ptr<ftSurface> face,
				   shared_ptr<Vertex> vx, 
				   std::vector<shared_ptr<Vertex> >& corners,
				   double epsge, Point& point, 
				   double& par, double& dist);

    Point getInVec(shared_ptr<Vertex> vx, shared_ptr<ftSurface> face);
  }

}  // namespace Go

#endif // _REGULARIZEUTILS_H
