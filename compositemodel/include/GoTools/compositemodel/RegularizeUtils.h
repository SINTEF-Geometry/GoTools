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
		   std::vector<shared_ptr<Vertex> > non_corner);

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
  }

}  // namespace Go

#endif // _REGULARIZEUTILS_H
