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
    std::vector<std::shared_ptr<ftSurface> > 
      divideVertex(std::shared_ptr<ftSurface> face,
		   std::shared_ptr<Vertex> vx, 
		   std::vector<std::shared_ptr<Vertex> > cand_vx,
		   ftEdge* cand_edge,
		   double epsge, double tol2, double angtol,
		   std::vector<std::shared_ptr<Vertex> > non_corner);

    std::vector<std::shared_ptr<ftSurface> > 
      createFaces(std::vector<std::shared_ptr<BoundedSurface> >& sub_sfs,
		  std::shared_ptr<ftSurface>  face,
		  double epsge, double tol2, double angtol,
		  std::vector<std::shared_ptr<Vertex> > non_corner);

    void getDivisionPlane(std::shared_ptr<ftSurface> face,
			  std::shared_ptr<Vertex> vx,
			  double epsge,
			  Point& pnt,
			  Point& normal);

    void 
      getClosestBoundaryPar(std::shared_ptr<ftSurface> face,
			    std::shared_ptr<Vertex> vx,
			    std::vector<std::shared_ptr<ParamCurve> >& vx_cvs,
			    const Point& pnt,
			    double epsge,
			    int& close_idx, double& close_dist,
			    Point& close_par);
    bool
      cornerInShortestPath(std::shared_ptr<Vertex> vx1,
			   std::shared_ptr<Vertex> vx2,
			   std::shared_ptr<ftSurface> face,
			   double angtol);

    bool
      getPath(ftEdge* edg, std::shared_ptr<Vertex> vx,
	      std::shared_ptr<Vertex> last, std::shared_ptr<ftSurface> face,
	      std::vector<ftEdge*>& path);
  }

}  // namespace Go

#endif // _REGULARIZEUTILS_H
