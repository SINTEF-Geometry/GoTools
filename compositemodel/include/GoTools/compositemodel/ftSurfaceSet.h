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

#ifndef _FTSURFACESET_COMP_H
#define _FTSURFACESET_COMP_H

//===========================================================================
//===========================================================================

#include <vector>             // Standard library STL vector
#include <string>             // Standard library string
#include "GoTools/compositemodel/ftMessage.h"        // Output message (warning or error)
#include "GoTools/igeslib/ftTangPriority.h"
#include "GoTools/compositemodel/ftPointSet.h"
#include "GoTools/topology/tpTolerances.h"
#include "GoTools/compositemodel/ftFaceBase.h"
#include "GoTools/compositemodel/ftEdgeBase.h"
#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/geometry/SplineCurve.h"
#include "GoTools/compositemodel/ttlPoint.h"
#include "GoTools/compositemodel/ttlTriang.h"
#include "GoTools/compositemodel/ftSurfaceSetPoint.h"
#include "GoTools/compositemodel/ftPlane.h"
#include "GoTools/compositemodel/ftSSfEdge.h"


namespace Go
{


  class ftSurfaceSet;
  // from_t, to_t, surf*, index (ordering: (vmin, vmax, umin, umax)).
  typedef std::pair<std::pair<double, double>,
    std::pair<ftFaceBase*, int> > BoundaryPiece;

  /** A set of surfaces is viewed as one surface.
   * Member functions assume object is handled by a tpTopologyTable.
   */
  class GO_API ftSurfaceSet : public ftFaceBase
  {
  public:
  
    /// Constructor
    ftSurfaceSet(const std::vector<shared_ptr<ParamSurface> >& surfaces,
		 tpTolerances& topeps, double approxeps);


    /// Destructor
    ~ftSurfaceSet();

    // Evaluation and interrogation. Overridden from ftFaceBase.
    virtual std::vector<shared_ptr<ftEdgeBase> > 
      createInitialEdges(double degenerate_epsilon = DEFAULT_SPACE_EPSILON,
			 double kink = 0.00015,
			 bool no_split = false) = 0;

/*     virtual std::vector<shared_ptr<ftEdgeBase> >  */
/*       setOrientation(double degenerate_epsilon = DEFAULT_SPACE_EPSILON); */

    virtual std::vector<shared_ptr<ftEdgeBase> > startEdges();

    virtual Point point(double u, double v) const = 0;
    virtual Point normal(double u, double v) const = 0;

    virtual BoundingBox boundingBox();

    virtual void setId(int id);
    virtual int getId();

/*     virtual void turnOrientation(); */
/*     virtual bool getOrientation(); */

    virtual shared_ptr<ParamSurface> surface()
      { return shared_ptr<ParamSurface>(surf_); }

    virtual ftMessage createSurf(double& max_error, double& mean_error) = 0;
    virtual void getError(double& max_error, double& mean_error) = 0;

    /// Priority type
    virtual ftTangPriority getPrioType() const = 0;

    /// Faces or parts of faces in faces_ are trimmed if they lie above plane.
    /// Remember to update top_table_ and first_edges_!
    ftMessage trimWithPlane(const ftPlane& plane);

    /// Return the surfaces belonging to the surface set
    void getSurfaces(vector<shared_ptr<ParamSurface> >& surfaces);

    // Modify the current surface. The surface may be representing the
    // spline space without containing coefficient data. Parameter iteration,
    // iteration on the weight on approximation compared to smoothness, and
    // refinement of the spline space is included.
    ftMessage modifySurface(int maxiter, int max_update_iter, 
			    std::vector<int>& edge_derivs,
			    ftPointSet& points,
			    double& max_error, double& mean_error);

    double getApproxweight()
    { return approxweight_; }

    void setApproxweight(double approxweight)
    { approxweight_ = approxweight; }

    void setAdditionalCornerPts(std::vector<Point>& add_corner_pts)
    { additional_corner_pts_ = add_corner_pts; }

    /// Closest point between this face and a point
    virtual void closestPoint(const Point& pt,
		      double&  clo_u,
		      double&  clo_v, 
		      Point& clo_pt,
		      double&  clo_dist,
		      double   epsilon) const
    {
      THROW("Not implemented");
    }
  
  protected:

    //------------------ Data members ----------------------
    tpTolerances toptol_;
    //tpTopologyTable<ftEdgeBase, ftFaceBase>* top_table_;
    std::vector<shared_ptr<ftFaceBase> > faces_;
    shared_ptr<SplineSurface> surf_;
    //std::vector<shared_ptr<ftEdgeBase> > first_edges_;
    /// The boundary loops limiting this face, the first loop represents
    /// the initial boundary
    std::vector<shared_ptr<Loop> > boundary_loops_;

    double approxtol_;
    double approxweight_;  // Weight used on point approximation in smoothing
    // We store approximated boundary curves. 
    std::vector<shared_ptr<SplineCurve> > approx_bd_curves_;
    //     bool is_turned_;
    //     int id_;
    std::vector<Point> additional_corner_pts_; // In space. User defined (assuming # corners in set < 4).

    // Protected member functions to be used by sub classes.

    // Approximate only the edges in the group of surfaces
    virtual bool sampleOnlyEdges() const = 0;

    // We expect the set of surfaces to define one boundary (outer) only.
    // Function returns :vector defining bd, with corners as given by corners.
    // kink_tol defines whether two consecutive edges meet in a kink (corner).
    // Assumes first_edges_ has been set.
    ftMessage getOuterLoop(std::vector<ftEdgeBase*>& outer_loop, std::vector<int>& corners,
			   double bend_tol, double kink_tol);

    // Sample data points representing the shape of the current
    // ftSurfaceSet
    ftMessage
      fetchSamplePoints(const std::vector<ftEdgeBase*>& edgeloop, std::vector<int>& corner,
			std::vector<int>& cn, ftPointSet& points);

    // Get initial boundary data points from the input surface set. Elements in
    // points->points_ are really of type ftSurfaceSetPoint. Must cast when using.
    // Number of samples either based on length or on curvature. Default = length.
    void getInitBndData(std::vector<ftEdgeBase*>& edgc, ftPointSet& points, int cn[]);

    // Get initial inner data points from the input surface set.
    // Routine only samples points, topology set by separate routine.
    // Version 2 samples along iso-lines in one direction, possibly differing #pts.
    ftMessage getInitInnerData2(ftPointSet& points);
    ftMessage getInitInnerData(ftPointSet& points, int max_sample);

    // Given collection of points where only bd points have some topological
    // information, we construct a topology for all points by using TTL.
    // Warning! Casting inside function (assumes points contain ftSurfaceSetPoint's).
    ftMessage updatePointTopology(ftPointSet& points);

    // Merge the input surfaces into one large surface.
    // Rutine calls subclasses different makeSurface routines.
    ftMessage merge(double& max_error, double& mean_error,
		    double bend_tol, double kink_tol);

    // Return approximative lengths of edges (oriented vmin, vmax, umin, umax).
    std::vector<double> estimateSurfSides(SplineSurface* surf);
    //     void estimateSurfSize(SplineSurface* surf,
    // 			  double& length1, double& length2);

    // We give surf_ parameter domain based on the geometry (lengths) given by
    // outer loop.
    ftMessage reparametrizeSurf(const std::vector<ftEdgeBase*>& first_edges,
				const std::vector<int>& corners,
				ftPointSet& points, double umin=0, double umax=0,
				double vmin=0, double vmax=0);

    // Fetch the curves and cross tangents from the boundary edges.
    // If require_cord_length_param == true, we use appr of cvs even if not within approxtol_.
    // Currently the gridding relies on the surface being cord length parametrized.
    ftMessage
      getBoundaryConditions(const std::vector<ftEdgeBase*>& edgeloop, 
			    std::vector<int>& corner, 
			    std::vector< shared_ptr<SplineCurve> >& bd_curves, 
			    std::vector< shared_ptr<SplineCurve> >& cross_curves,
			    bool compute_cross_curves,
			    std::vector<BoundaryPiece>& bdpiece,
			    bool require_cord_length_param = false);

    // As long as the boundary conditions are uniform (no neighbouring 
    // surface or the same neighbouring surface), compute the size 
    // of the current type of condition, and return the condition.
    // cv is newed and handled by smart ptr.
    ftMessage getNextEdgePieceInfo(std::vector<ftEdgeBase*>::iterator& first_edge,
				   std::vector<ftEdgeBase*>::iterator& last_edge,
				   ftFaceBase*& adjsurf, SplineCurve*& cv,
				   Point& ptmin, Point& ptmax,
				   // 				 double& tmin, double& tmax,
				   double& length);

    // Join cross boundary curves along an edge into one curve, and
    // create missing curves. If no curve exist return a dummy shared pointer.
    shared_ptr<SplineCurve> 
      joinCrossCurves(std::vector< shared_ptr<SplineCurve> >& cross_curves,
		      double startpar, double endpar);

    // Return a part of the boundary of the current surface if the
    // boundary is stored separately
    shared_ptr<SplineCurve> 
      getBoundaryPiece(Point& ptmin, Point& ptmax, double eps);

  private:

    // Construct the topology information regarding the input geometry,
    // based on the information computed in top_table_.
    void buildTopology();

    //  Fetch all boundary loops of the current ftSuperSurface. Start
    //  from first_edges_.
    void BoundaryLoops(std::vector< std::vector<ftEdgeBase*> >& loopvec);

    // Create a 4-sided surface approximating the input surfaces, and interpolate
    // the boundary (if it is sufficiently smooth). Called by subclasses.
    // loop describes outer boundary of surface set, corner indicates corner edges.
    virtual ftMessage makeSurface(const std::vector<ftEdgeBase*>& loop,
				  std::vector<int>& corner,
				  double& max_error, double& mean_error) = 0;

    virtual int nmbToEval(ftEdge* edge, double tmin, double tmax) = 0;

    // Add point to points, update structure, called by getInitBndData().
    void addBndPoint(shared_ptr<ftSurfaceSetPoint>& ftpnt,
		     PointIter& prevpt, ftPointSet& points);

    // For given point, we add reference to face in point (we also add par value).
    // Called by getInitBndData() & getInitInnerData().
    void addFaceRef(ftSurfaceSetPoint& point, shared_ptr<ftFaceBase>& face,
		    double* seed);

    /// Sample interior points of an edge, called by getInitBndData().
    /// nmb_eval refers to total sample points along edge. Replace?
    void getEdgeInnerData(ftEdgeBase* curr_edge, PointIter& prevpt,
			  ftPointSet& points, std::vector<ftEdgeBase*>& edgc,
			  int cn[], bool set_second, int nmb_eval);

    // We require that distance from interior points (sampled) to straight line
    // segment (defined by points in tmin and tmax) is within max_dist.
    // By recursion we decide # times we have to split the edge in order to satisfy
    // requirement (sampling is to be performed uniformly in parameter domain).
    int nmbSampleSegments(const ftEdgeBase& edge,
			  double tmin, double tmax, double max_dist);

    // Collecting all pts in pt_set w/index matching idx. Assuming that all pts are sampled
    // from bd of sf given by idx; we may then traverse that bd given a sample_pt on bd.
    // We may assume that members in pt_set are connected to to other pts with the same index.
    std::vector<Vector2D> getBdTrain(ftPointSet& pt_set, shared_ptr<ftFaceBase>& face);

    // Returned par_pt is NULL unless edge is a CurveOnSurface.
    void getParPoint(ftEdgeBase* edge, double tpar, shared_ptr<Point>& par_pt);

    // Given input of a general parametric surface, return a containing surface.
    // If input is a SplineSurface, we return the copy. For a trimmed surface we find
    // boundingbox of trim cv, then return the subsurface of the underlying SplineSurface.
    SplineSurface* createContainingSurface(ParamSurface* sf);

  };

} // namespace Go


#endif // _FTSURFACESET_COMP_H
 
