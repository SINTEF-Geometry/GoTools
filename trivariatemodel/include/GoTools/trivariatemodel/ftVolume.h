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

#ifndef _FTVOLUME_H
#define _FTVOLUME_H

#include "GoTools/trivariate/ParamVolume.h"
#include "GoTools/compositemodel/SurfaceModel.h"
#include "GoTools/compositemodel/Body.h"

namespace Go
{

  class SurfaceOnVolume;
  class ParamCurve;

  /// Struct to store information about adjacency relations between two bodies
  struct VolumeAdjacencyInfo
  {
  public:
    /// True if adjacency was found, false if not
    bool adjacency_found_;
    /// True if two volumes meet in along an edge
    bool corner_adjacency_;
    /// Boundary index of first volume. Value 0-5 in order umin, umax, vmin, vmax, wmin, wmax
    int bd_idx_1_;
    /// Boundary index of second volume. Same interpretation as bd_idx_1_
    int bd_idx_2_;
    /// Set if corner_adjacency_ is true. 0-3 in the order umin, umax, vmin, vmax
    int edg_idx_1_;
    /// Set if corner_adjacency_ is true. 0-3 in the order umin, umax, vmin, vmax
    int edg_idx_2_;
    /// True if common boundaries are equally oriented in first parameter direction of first surface
    bool same_orient_u_;
    /// True if common boundaries are equally oriented in second parameter direction of first surface
    bool same_orient_v_;
    // same_orient_w: (bd_idx_1_ + bd_idx_2_)%2 == 1
    /// True if corner_adjacency_ is found and the boundary curves are equally oriented 
    bool same_orient_edge_;
    /// True if u-directions concide for both boundary surfaces (regardless of orientation)
    bool same_dir_order_;
    /// True if adjacency is found, user wanted to test if the surfaces meet
    /// in corner-to-corner configuration, and the test failed. Otherwise false.
    bool corner_failed_;

    /// Constructor
    VolumeAdjacencyInfo() :
    adjacency_found_(false), corner_adjacency_(false), corner_failed_(false)
      { }
    };  // End struct AdjacencyInfo


  /// \brief A topological solid with a trivariate geometry description

  class ftVolume : public Body
  {
  public:
    /// Constructor given a trivariate volume description
    ftVolume(shared_ptr<ParamVolume> vol, int id=-1);

    /// Constructor given a trivariate volume description and tolerances
    /// for topology analysis of the boundary shell
    ftVolume(shared_ptr<ParamVolume> vol, double gap_eps,
	     double kink_eps, int id=-1);

    ftVolume(shared_ptr<ParamVolume> vol, double gap_eps, double neighbour,
	     double kink_eps, double bend, int id=-1);

    /// Given a volume and boundary surfaces, create a possibly trimmed
    /// ftVolume
    ftVolume(shared_ptr<ParamVolume> vol, 
	     shared_ptr<SurfaceModel> shell,
	     int id=-1);

    /// Given a volume and a number of boundary shells, create a possibly trimmed
    /// ftVolume
    ftVolume(shared_ptr<ParamVolume> vol, 
	     std::vector<shared_ptr<SurfaceModel> > shells,
	     int id=-1);

    /// Create a ftVolume when no geometry description is known yet
    ftVolume(shared_ptr<SurfaceModel> shell,
	     int id=-1);

     /// Destructor
    ~ftVolume();

    /// Fetch geometric description
    shared_ptr<ParamVolume> getVolume()
      {
	return vol_;
      }

    /// Fetch Id, not necessarily uniquely set
    int getId()
    {
      return id_;
    }

    /// The bounding box corresponding to this solid
    virtual BoundingBox boundingBox() const;

    /// Check if the volume is represented as a spline
    bool isSpline() const;

    /// Fetch all radial edges belonging to this model
    std::vector<shared_ptr<EdgeVertex> > radialEdges() const;

    /// Fetch all radial edges common to this model and another one
    std::vector<shared_ptr<EdgeVertex> > 
      getCommonEdges(ftVolume *other) const;

    /** Fetch edges where no radial edge exist, i.e. there is no 
	ajacent volume along any boundary surface meeting in this edge.
	Only one occurance in an edge, twin-edge pair is returned */
    std::vector<shared_ptr<ftEdge> > uniqueNonRadialEdges() const;

    /// Get neighbouring bodies
    void getAdjacentBodies(std::vector<ftVolume*>& neighbours);

    /// Check if two bodies are neighbours, are splines and have a common
    /// spline space at the interface
    bool commonSplineSpace(ftVolume *other, double tol);

    /// Ensure that two neighbouring spline volumes have a common
    /// spline space at the interface
    bool makeCommonSplineSpace(ftVolume *other);

    /// Check if two volumes are splines and meet in a corner-to-corner
    /// configuration
    bool isCornerToCorner(shared_ptr<ftVolume> other,
			  double tol);

    /// Ensure that two spline volumes meet in a corner to corner
    /// configuration
    void splitAtInternalCorner(ftVolume* other,
			    std::vector<shared_ptr<ftVolume> >& new_vol1,
			    std::vector<shared_ptr<ftVolume> >& new_vol2,
			    double tol=DEFAULT_SPACE_EPSILON);

    /// DEPRECATED METHOD. USE THE OTHER getAdjacencyInfo() instead
    /// Fetch info on adjacency between neighbouring bodies
    /// return value true if adjacency is found and both bodies are non-trimmed or
    /// boundary trimmed
    /// bd1: boundary index of first volume
    ///  0 = umin, 1 = umax, 2 = vmin,  3 = vmax, 4 = wmin, 5 = wmax
    /// bd2: boundary index of second volume
    ///  0 = umin, 1 = umax, 2 = vmin,  3 = vmax, 4 = wmin, 5 = wmax
    /// orientation: Indicates whether the surfaces are equally oriented along
    /// the common boundary
    /// 0 = same orientation
    /// 1 = first parameter of the common boundary surface is opposite
    /// 2 = second parameter of the common boundary surface is opposite
    /// 3 = both parameters of the common boundary surface is opposite
    /// same_seq: Indicates if the two parameter directions of the common
    ///           boundary surface are given in the same sequence
    bool getAdjacencyInfo(ftVolume *other, double tol,
			  int& bd1, int& bd2, 
			  int& orientation, bool& same_seq);

    /// Fetch info on adjacency between neighbouring bodies
    VolumeAdjacencyInfo getAdjacencyInfo(ftVolume *other, double tol,
					 int adj_idx = 0, bool test_corner = false);

    VolumeAdjacencyInfo getCornerAdjacencyInfo(ftVolume *other, 
					       EdgeVertex* evx,
					       double tol, int adj_idx=0);

    /// Fetch info on adjacence between bodies meeting an a common
    /// edge (no common face)
    VolumeAdjacencyInfo getCornerAdjacencyInfo(ftVolume *other, double tol,
					       int adj_idx = 0);

    /// Look for adjacencey betweed degenerate volumes in a degenerate
    /// boundary surface
    bool checkDegAdjacency(ftVolume *other, 
			   shared_ptr<EdgeVertex> evx,
			   double tol,
			   shared_ptr<ParamSurface>& bdsf1, 
			   shared_ptr<ParamSurface>& bdsf2);

    /// Given two adjacent spline volumes represented in the same spline space,
    /// return local enumeration of the coefficients. This surface is given 
    /// first and the other volume second
    bool getCorrCoefEnumeration(ftVolume *other, double tol,
				std::vector<std::pair<int,int> >& enumeration);

    /// Given a vertex, find if this vertex is associated this body and compute
    /// the parameter value of the body associated with the vertex
    bool getVertexPosition(shared_ptr<Vertex> vx, Point& param) const;

    /// Given a vertex, find if this vertex is associated this body, compute
    /// the parameter value of the body associated with the vertex, and if this
    /// body is a spline volume, return the corner number and the number 
    /// of the associated spline coefficient
    /// corner: 0=(umin,vmin,wmin), 1=(umax,vmin,wmin), 2=(umin,vmax,wmin),
    ///         3=(umax,vmax,wmin), 4=(umin,vmin,wmax), 5=(umax,vmin,wmax),
    ///         6=(umin,vmax,wmax), 7=(umax,vmax,wmax)
    bool getVertexEnumeration(shared_ptr<Vertex> vx, 
			      Point& param, int& corner,
			      int& coef_nmb) const;

    /// Return the boundary number of outer boundary surfaces that follow 
    /// volume boundaries
    /// return value = true: All free boundaries are associated with volume
    ///                      boundaries
    ///              = false:Some free boundaries are not associated with 
    ///                      volume boundaries
    bool getFreeBoundaryInfo(double tol, std::vector<int>& free_boundaries);

    /// For a spline volume, get the local enumeration of coefficients
    /// along a specified boundary
    /// return value: true if spline, false if not
    bool getBoundaryCoefEnumeration(int bd, std::vector<int>& enumeration);

    /// Information about whether or not the volume is trimmed and how it
    /// is trimmed
    /// Check if the volume is boundary trimmed (not trimmed). The boundary
    /// surfaces themselves may be trimmed
    bool isBoundaryTrimmed() const;

    /// Check if all boundary surfaces correspond to an iso-parameter in
    /// the volume
    bool isIsoTrimmed() const;

    /// Check if all boundary surfaces spproximately correspond to an 
    /// iso-parameter in the volume
/*     bool isIsoTrimmed(double eps) const; */

    /// Update the volume by regularizing all boundary shells,
    /// i.e. all faces in all boundary shells of all connected volumes
    /// should be 4-sided, and no T-joints are allowed
    bool regularizeBdShells(std::vector<std::pair<Point,Point> >& corr_vx_pts,
			    std::vector<SurfaceModel*>& modified_adjacent,
			    int split_mode = 1, bool pattern_split = false);

    /// Check if this volume has 6 boundary surfaces that may act
    /// as the boundary surfaces of a non-trimmed spline volume
    bool isRegularized() const;

    /// Modify a regularized, possibly trimmed volume to become non-trimmed
    bool untrimRegular(int degree);

/*     /// Split this and the corresponding volume with regard to the */
/*     /// intersections between the boundary surfaces corresponding to */
/*     /// these two volumes */
/*     std::vector<shared_ptr<ftVolume> > */
/*       splitVolumes(shared_ptr<ftVolume> other, double eps); */

    /// Divide a trimmed volume into a set of regular volumes
    /// NB! The boundary shells of the corresponding volume model (or this
    /// volume) must be regular.
    /// Not all configurations are handled. If an unknown configuration
    /// occur, nothing is returned
    /// Should be called from VolumeModel.
    /// Ruins the current ftVolume.
    std::vector<shared_ptr<ftVolume> > 
      replaceWithRegVolumes(int degree,
			    std::vector<SurfaceModel*>& modified_adjacent,
			    bool performe_step2=true,
			    int split_mode=1,
			    bool pattern_split=false);
    
    /// Update boundary shells to reflect changes in the geometric volume
    /// while maintaining topology information
    void 
      updateBoundaryInfo();

    /// Debug
    bool checkBodyTopology();


  private:
     /// Geometric description of volume
    shared_ptr<ParamVolume> vol_;
    int id_;

    std::vector<shared_ptr<ftEdge> > missing_edges_;  // Private storage

    /// Private method to create the boundary shell
    shared_ptr<SurfaceModel> 
      createBoundaryShell(double eps, double tang_eps);

    std::vector<shared_ptr<ftSurface> >
      getBoundaryFaces(shared_ptr<ParamVolume> vol,
		       double eps, double tang_eps);

    /// Sort boundary faces in a regular ftVolume
    bool 
      sortRegularSurfaces(std::vector<shared_ptr<ParamSurface> >& sorted_sfs,
			  std::vector<std::pair<int,double> >& classification);

    shared_ptr<SurfaceOnVolume> 
      getVolSf(shared_ptr<ParamSurface>& surf) const;
    
    std::vector<std::pair<int, double> >
      getMidCurveIntersections(shared_ptr<ParamCurve> curve,
			       std::vector<shared_ptr<ParamSurface> >& sfs,
			       double tol) const;

    shared_ptr<ParamVolume> 
      createByLoft(shared_ptr<ParamSurface> sf1,
		   shared_ptr<ParamSurface> sf2, 
		   double tol,  int pardir);

    shared_ptr<ParamVolume> 
      createByCoons(std::vector<shared_ptr<ParamSurface> >& sfs,
		    std::vector<std::pair<int,double> >& classification,
		    double tol, int degree);

    bool
      getCoonsCurvePairs(std::vector<shared_ptr<ParamSurface> >& sfs, 
			 std::vector<std::vector<std::pair<shared_ptr<ParamCurve>,shared_ptr<ParamCurve> > > >& curves,
			 std::vector<std::vector<int> >& indices);

    void getCoonsBdCurves(std::vector<std::pair<shared_ptr<ParamCurve>,shared_ptr<ParamCurve> > >& cvs,
			  std::vector<int>& indices,
			  std::vector<std::pair<int,double> >& classification,
			  double tol, int degree,
			  std::vector<shared_ptr<SplineCurve> >& coons_cvs);
    
    std::vector<shared_ptr<ftSurface> >  
      generateMissingBdSurf(int degree,
			    std::vector<std::pair<Point,Point> >& corr_vx_pts,
			    bool perform_step2, bool smooth_connections);

    void makeSurfacePair(std::vector<ftEdge*>& loop,
			 int degree,
			 shared_ptr<ftSurface>& face1,
			 shared_ptr<ftSurface>& face2,
			 std::vector<std::pair<ftEdge*,ftEdge*> >& replaced_wires);

    void getEdgeCurves(std::vector<ftEdge*>& loop, 
		       std::vector<shared_ptr<ParamCurve> >& space_cvs,
		       std::vector<Point>& joint_points);

    ftEdge*  getLeftLoopEdge(ftSurface* face, Body *bd,
			     shared_ptr<EdgeVertex> radial);

    bool  doSwapEdges(ftSurface* face, ftEdge* edge1, ftEdge *edge2);

    std::vector<std::vector<ftEdge*> > 
      getMissingSfLoops(std::vector<std::pair<Point,Point> >& corr_vx_pts,
			bool perform_step2, bool smooth_connections);

    bool loopExisting(std::vector<ftEdge*>& loop, 
		      std::vector<std::vector<ftEdge*> >& curr_loops);

    std::vector<shared_ptr<ftEdge> > getStartEdges();

    std::vector<std::vector<ftEdge*> > getLoop(shared_ptr<ftEdge> start_edge);
    
    bool getLoopEdges(std::vector<ftEdge*>& loop, 
		      shared_ptr<Vertex> start_vx,
		      shared_ptr<Vertex> vx,
		      int max_nmb=4);

    bool sameFace(std::vector<ftEdge*>& loop);

    bool checkPlaneLoop(std::vector<ftEdge*>& loop);

    std::vector<shared_ptr<ftVolume> > 
      createRegularVolumes(std::vector<shared_ptr<ftSurface> > bd_faces);

    void sortCoonsPatchBdCvs(std::vector<shared_ptr<ParamCurve> >& cvs,
			     std::vector<shared_ptr<ParamCurve> >& space_cvs,
			     double tol);

    void
      moveVolParCv(shared_ptr<ParamCurve>& pcv,
		   shared_ptr<ParamCurve>& spacecv,
		   const Point& dir, double tol);

    void 
      getCurrConnectedModel(std::vector<shared_ptr<ftSurface> >& face,
			    size_t idx,
			    std::vector<shared_ptr<ftSurface> >& curr_set,
			    std::vector<shared_ptr<ftSurface> >& all_sets) const;

    void 
      replaceParamVolume(shared_ptr<ParamVolume> vol, 
			 std::vector<shared_ptr<ParamSurface> >& sorted_sfs,
			 bool loft_sequence);
    int 
      findFaceMatch(shared_ptr<ftSurface> face,
		    std::vector<shared_ptr<ftSurface> >& cand_matches);

    std::vector<std::pair<int, int> >  
      oppositeSfs(shared_ptr<SurfaceModel> model);

    void  removeSeamFaces();

    void eraseMissingEdges();

    shared_ptr<ParamCurve> makeMissingEdgeCv(shared_ptr<Vertex> vx1,
					     shared_ptr<Vertex> vx2);

    void simplifyOuterBdShell(int degree);

    int mergeSituation(ftSurface* face1, ftSurface* face2,
		       shared_ptr<Vertex> vx1, shared_ptr<Vertex> vx2,
		       int& dir1, double& val1, bool& atstart1, 
		       int& dir2, double& val2, bool& atstart2, 
		       std::pair<Point, Point>& co_par1, 
		       std::pair<Point, Point>& co_par2);

   std::vector<ftSurface*> getMergeCandFaces(shared_ptr<ftSurface> curr,
					     std::vector<std::pair<shared_ptr<Vertex>,
					      shared_ptr<Vertex> > >& common_vxs);

   void estMergedSfSize(ftSurface* face1, ftSurface* face2,
			shared_ptr<Vertex> vx1,shared_ptr<Vertex> vx2,
			double& len_frac, double& other_frac, double& sf_reg);
  };



} // namespace Go


#endif // _FTVOLUME_H
