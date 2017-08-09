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

#ifndef _REGULARIZEUTILS_H
#define _REGULARIZEUTILS_H

#include "GoTools/compositemodel/ftSurface.h"
#include "GoTools/compositemodel/ftEdge.h"
#include "GoTools/geometry/BoundedSurface.h"

namespace Go {

  /// Utility functionality for RegularizeFace and RegularizeFaceSet   
  namespace RegularizeUtils
  {
    std::vector<shared_ptr<ftSurface> > 
      divideVertex(shared_ptr<ftSurface> face,
		   shared_ptr<Vertex> vx, 
		   std::vector<shared_ptr<Vertex> >& cand_vx,
		   ftEdge* cand_edge,
		   std::vector<shared_ptr<Vertex> >& prio_vx,
		   double epsge, double tol2, double angtol,
		   double bend,
		   std::vector<shared_ptr<Vertex> >& non_corner,
		   const Point& centre, const Point& axis,
		   bool strong = false);

    std::vector<shared_ptr<CurveOnSurface> > 
      findVertexSplit(shared_ptr<ftSurface> face,
		      shared_ptr<Vertex> vx, 
		      std::vector<shared_ptr<Vertex> >& cand_vx,
		      ftEdge* cand_edge,
		      std::vector<shared_ptr<Vertex> >& prio_vx,
		      double epsge, double tol2, double angtol,
		      double bend,
		      std::vector<shared_ptr<Vertex> >& non_corner,
		      const Point& centre, const Point& axis,
		      shared_ptr<BoundedSurface>& bd_sf,
		      bool strong = false);

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
			    Point& close_par, int loop_idx=-1);
    bool
      cornerInShortestPath(shared_ptr<Vertex> vx1,
			   shared_ptr<Vertex> vx2,
			   shared_ptr<ftSurface> face,
			   double angtol);

    bool
      getPath(ftEdge* edg, shared_ptr<Vertex> vx,
	      shared_ptr<Vertex> last, shared_ptr<ftSurface> face,
	      std::vector<ftEdge*>& path);

    int
      noExtension(shared_ptr<Vertex> vx, ftSurface* face,
		  shared_ptr<Vertex>& vx2, std::pair<Point, Point>& co_par1, 
		  std::pair<Point, Point>& co_par2, int& dir1, int& dir2,
		  double& val1, double& val2, double angtol, bool check_constant_curve);

    bool
      mergeSituationContinuation(ftSurface* init_face, shared_ptr<Vertex> vx,
				 ftEdge* edge, double angtol);

    double getMaxParFrac(shared_ptr<ftSurface> face);

    int selectCandVx(shared_ptr<ftSurface> face,
		     shared_ptr<Vertex> vx, const Point& in_vec, 
		     vector<shared_ptr<Vertex> > cand_vx,
		     RectDomain& dom,
		     double epsge, double angtol,
		     const Point& centre, const Point& normal,
		     std::vector<shared_ptr<ParamCurve> >& vx_cvs,
		     double close_dist, const Point& close_pt,
		     double& cyl_rad, bool strong=false);

    void 
      adjustTrimSeg(std::vector<shared_ptr<CurveOnSurface> >& trim_segments,
		    Point *parval1, Point *parval2,
		    shared_ptr<ftSurface> face,
		    shared_ptr<BoundedSurface>& bd_sf,
		    std::vector<shared_ptr<Vertex> >& non_corner,
		    double tol, double epsge);

    void checkTrimSeg(std::vector<shared_ptr<CurveOnSurface> >& trim_segments,
		      std::vector<shared_ptr<Vertex> >& next_vxs,
		      const Point& vx_point, const Point& other_pt,
		      double epsge);

    void 
      checkTrimSeg2(std::vector<shared_ptr<CurveOnSurface> >& trim_segments,
		    const Point& vx_par1, const Point& vx_par2, 
		    double epsge);

    void 
      checkTrimSeg3(std::vector<shared_ptr<CurveOnSurface> >& trim_segments,
		    const Point& vx_par1, const Point& vx_par2, 
		    double epsge);

    void 
      checkTrimConfig(shared_ptr<ftSurface> face,
		      std::vector<shared_ptr<CurveOnSurface> >& trim_segments,
		      shared_ptr<Vertex> vx,
		      std::vector<shared_ptr<Vertex> >& corners,
		      double epsge);

    ftEdge* getOppositeBoundaryPar(shared_ptr<ftSurface> face,
				   shared_ptr<Vertex> vx, 
				   std::vector<shared_ptr<Vertex> >& corners,
				   double epsge, Point& point, 
				   double& par, double& dist);

    Point getInVec(shared_ptr<Vertex> vx, shared_ptr<ftSurface> face);

    shared_ptr<ParamCurve> checkStrightParCv(shared_ptr<ftSurface> face,
					     const Point& pos1,
					     const Point& pos2,
					     double epsge);
    shared_ptr<ParamCurve> checkStrightParCv(shared_ptr<ftSurface> face,
					     shared_ptr<Vertex> vx1, 
					     shared_ptr<Vertex> vx2,
					     double epsge);
    shared_ptr<ParamCurve> checkStrightParCv(shared_ptr<ftSurface> face,
					     shared_ptr<Vertex> vx1, 
					     const Point& mid,
					     double epsge);
    bool
      checkPath(shared_ptr<Vertex> vx1, shared_ptr<Vertex> vx2,
		shared_ptr<Vertex> vx, shared_ptr<ftSurface> face,
		double angtol);

    bool checkRegularity(std::vector<shared_ptr<Vertex> >& cand_vx,
			 shared_ptr<ftSurface> face,
			 bool checkConvex = true);

    std::vector<shared_ptr<Vertex> > endVxInChain(shared_ptr<ftSurface> face,
						  ftSurface* face1,
						  ftSurface* face2,
						  shared_ptr<Vertex> vx,
						  shared_ptr<Vertex> prev,
						  shared_ptr<Vertex> vx0,
						  std::vector<shared_ptr<Vertex> >& met_already);

    int traverseUntilTJoint(std::vector<ftSurface*> vx_faces,
			    shared_ptr<Vertex> vx,
			    shared_ptr<Vertex>& vx2,
			    std::vector<ftSurface*>& vx_faces2);

    void angleInEndpoints(shared_ptr<CurveOnSurface> seg,
			  shared_ptr<Vertex> vx1, 
			  shared_ptr<Vertex> vx2,
			  shared_ptr<ftSurface> face,
			  double& min_ang1, double& min_ang2);

    void getSourceCvs(std::vector<shared_ptr<ftEdge> >& all_edg,
		      std::vector<shared_ptr<ParamCurve> >& all_cvs);
  }

}  // namespace Go

#endif // _REGULARIZEUTILS_H
