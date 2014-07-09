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
/// This class is intended for use in block structuring. A set of faces  with possible 
/// inner and outer trimming is split according to certain rules to result in a
/// face set with 4 sided faces although faces with less than 4 sides can occur.
/// A side is defined as a piece of the face boundary between two corners or between
/// vertices where there are more than one adjacent face.
/// The trimmed surfaces being output from this class can later be approximated by
/// spline surfaces.
class RegularizeFaceSet
{
 public:
  /// Constructor
  RegularizeFaceSet(std::vector<shared_ptr<ftSurface> > faces, 
		    double epsge, double angtol, 
		    bool split_in_cand = false);

  RegularizeFaceSet(std::vector<shared_ptr<ftSurface> > faces, 
		    double gap, double neighbour, 
		    double kink, double bend, 
		    bool split_in_cand = false);
   /// Constructor
  RegularizeFaceSet(shared_ptr<SurfaceModel> model, bool split_in_cand = false);
  /// Destructor
  ~RegularizeFaceSet();

  void setSplitMode(int split_mode)
  {
    split_mode_ = split_mode;
  }

  /// Set information
  void setFaceCorrespondance(int idx1, int idx2);

  /// Fetch result
  std::vector<shared_ptr<ftSurface> > getRegularFaces(bool reverse_sequence=false);

  /// Return the resulting face set as a surface model
  shared_ptr<SurfaceModel> getRegularModel(bool reverse_sequence=false);

  /// Fetch info about point corrspondance
  std::vector<std::pair<Point, Point> > fetchVxPntCorr()
    {
      return corr_vx_pts_;
    }

  /// Fetch info about adjacent surface models that are modified
  /// due to changes in the current model
  std::vector<SurfaceModel*> getModifiedAdjacentModels()
    {
      return modified_models_;
    }

  private:
  shared_ptr<SurfaceModel> model_;

  std::vector<std::pair<int,int> > corr_faces_;

  int split_mode_;
  bool split_in_cand_;
  std::vector<std::vector<std::pair<std::pair<Point, int>,
    std::pair<Point,int> > > > cand_split_;

  std::vector<Point> seam_joints_;

  std::vector<std::pair<Point,Point> > corr_vx_pts_;

  // Information about adjacent surface models that are
  // updated during the regularization of the current model
  std::vector<SurfaceModel*> modified_models_;

    // Perform division
  void divide(bool reverse_sequence);

  void splitInTJoints();

  std::vector<shared_ptr<ftSurface> > 
    divideInTjoint(shared_ptr<ftSurface>& face,
		   std::vector<shared_ptr<Vertex> >& Tvx,
		   std::vector<shared_ptr<Vertex> >& corner,
		   bool& changed);


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


  std::vector<std::pair<std::pair<Point,int>, std::pair<Point,int> > > 
    getEndSplit(shared_ptr<ftSurface> prev_face,
		std::vector<shared_ptr<ftSurface> >& faces);

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

  int mergeSituation(ftSurface* face, 
		      shared_ptr<Vertex> vx,
		      ftSurface*& merge1,
		      int& dir1, double& val1, bool& atstart1,
		      ftSurface*& merge2,
		      int& dir2, double& val2, bool& atstart2,
		      shared_ptr<Vertex>& other_vx,
		      std::pair<Point, Point>& co_par1,
		      std::pair<Point, Point>& co_par2);

void prioritizeFaces(std::vector<shared_ptr<ftSurface> >& faces,
		     std::vector<int>& perm);
void
  computeFaceCorrespondance(std::vector<shared_ptr<ftSurface> >& faces);

};

}  // namespace Go

#endif // _REGULARIZEFACESET_H
