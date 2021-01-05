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

#ifndef LRVOLSTITCH_H
#define LRVOLSTITCH_H

#include "GoTools/geometry/ParamSurface.h"
#include "GoTools/lrsplines3D/LRSplineVolume.h"
#include "GoTools/lrsplines3D/Direction3D.h"
#include "GoTools/lrsplines3D/Mesh3D.h"
#include "GoTools/lrsplines3D/LRBSpline3D.h"

namespace Go
{
  /// This class stitches a collection of LR B-spline surfaces to
  /// ensure C0 or C1 continuity between adjacent surfaces. The surfaces
  /// must be organized in a regular grid. Not all surfaces in the grid
  /// must exist.
  /// The parameter domain of the surfaces must join such that adjacent
  /// surfaces have adjacent, none-overlapping parameter domains without
  /// any gaps.
  /// The surfaces must be oriented consistently
  /// C1 continuity is implemented only for 1D surfaces
  /// Due to the restrictions on the surface configurations only the
  /// organizing function is made public. Otherwise, general functionality
  /// is made private as they do not check if the configuration 
  /// constraints are satisfied.

  class LRVolStitch
  {
  public:
    /// Perform stitching. The volumes are organized from front to back, bottom to
    /// top and from left to right
    void stitchRegVols(std::vector<shared_ptr<LRSplineVolume> >& vols,
		       int nmb_u, int nmb_v, int nmb_w, double eps,
		      int cont);
    /// We calculate the max distance (cont = 0) / angle (cont = 1) between corresponding edges.
    std::vector<double> analyzeContinuity(std::vector<shared_ptr<LRSplineVolume> >& vols,
					  int nmb_u, int nmb_v, int nmb_w, int cont,
					  int num_edge_samples = 20);

  private:

    struct VolCorner
    {
      VolCorner(shared_ptr<LRSplineVolume>& vol, int uix, int vix, int wix)
      {
	vol_ = vol;
	uix_ = uix;
	vix_ = vix;
	wix_ = wix;
      }

      shared_ptr<LRSplineVolume> vol_;
      int uix_, vix_, wix_;
    };

    
    struct VolEdge
    {
      VolEdge(shared_ptr<LRSplineVolume>& vol, Direction3D dir, int ix1, int ix2)
      {
	vol_ = vol;
	dir_ = dir;
	ix1_ = ix1;
	ix2_ = ix2;
      }

      shared_ptr<LRSplineVolume> vol_;
      Direction3D dir_;
      int ix1_, ix2_;
    };

    struct VolBdSf
    {
      VolBdSf(shared_ptr<LRSplineVolume>& vol, Direction3D dir, int ix)
      {
	vol_ = vol;
	dir_ = dir;
	ix_ = ix;
      }

      shared_ptr<LRSplineVolume> vol_;
      Direction3D dir_;
      int ix_;
    };

    
    void consistentSplineSpaces(std::vector<shared_ptr<LRSplineVolume> >& vols,
				int nmb_u, int nmb_v, int nmb_w, double eps,
				int cont);

    // Make sure that the refinements along all edges have a full tensor product structure for
    // the first 'element_width' elements.
    // For corners additional knots are inserted.
    void tensorStructure(shared_ptr<LRSplineVolume> vol, int element_width,
			 bool bdsfs[6]);

    bool matchSplineSpace(shared_ptr<LRSplineVolume> vol1, int bdsf1,
			  shared_ptr<LRSplineVolume> vol2, int bdsf2, 
			  int element_width, double tol);

    void checkCornerMatch(std::vector<VolCorner>& volc, 
			 double tol, std::vector<int>& nmb_match,
			 std::vector<Point>& corner_val, 
			 std::vector<Point>& param);
    int averageCorner(std::vector<VolCorner>& volc, double tol);

   int makeCornerC1(std::vector<VolCorner>& volc,
		    double tol);

   bool averageEdge(std::vector<VolEdge>& voledg,
		     int cont, double tol);

   bool averageBdSf(std::vector<VolBdSf>& volsf,
		     int cont, double tol);

   void fetchEdgeCorners(const VolEdge& voledg, Point& c1, Point& c2);

   void fetchSfCorners(shared_ptr<LRSplineVolume> vol, int bdsf, Point corner[]);

   void extractBoundaryBsplines(const VolEdge& voledg, int cont,
				std::vector<std::vector<LRBSpline3D*> >& bsplines);
   void extractBdSfBsplines(const VolBdSf& volsf, int cont,
			    std::vector<std::vector<LRBSpline3D*> >& bsplines);

   void defineRefinements(const Mesh3D& mesh1, const Mesh3D& mesh2,
			  Direction3D dir,
			  std::vector<std::pair<std::vector<GPos2D>,int> >& mrects,
			  bool at_start, int width, int ix,
			  std::vector<LRSplineVolume::Refinement3D>& refs);

   void adaptMeshRectangles(std::vector<std::pair<std::vector<GPos2D>,int> >& rects,
			    int ix, int start, int end, int t1, int t2, int del,
			    std::vector<std::pair<std::vector<GPos2D>,int> >& rects2);
   
   void writeElems(std::ostream& out, shared_ptr<LRSplineVolume>& vol);
   void writeBSplines(std::ostream& out, std::vector<LRBSpline3D*>& bsplines);
   void writeMrects(std::ostream& out, const Mesh3D& mesh,
		    std::vector<std::pair<std::vector<GPos2D>,int > >& mrects,
		    Direction3D dir, int ix);
   void makeCrossC1(std::vector<LRBSpline3D*>& bsp, Direction3D dir);
    
   void makeLineC1(LRBSpline3D* bsp[4], Direction3D dir);

    // Utility function for debugging.
    // Compute distance, tangent and normal differences in input parameter points.
    void volumeDifference(const LRSplineVolume* sf1, const Point& param1,
			  const LRSplineVolume* sf2, const Point& param2,
			  double& vol_dist, double& tang1_ang,
			  double& tang2_ang, double& tang3_ang);

};
};

#endif
