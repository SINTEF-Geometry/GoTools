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

#ifndef _FTCHARTSURFACE_COMP_H
#define _FTCHARTSURFACE_COMP_H

#include <vector>             // Standard library STL vector
#include <string>             // Standard library string
#include "GoTools/compositemodel/ftMessage.h"        // Output message (warning or error)
#include "GoTools/igeslib/ftTangPriority.h"
#include "GoTools/compositemodel/ftPointSet.h"
#include "GoTools/compositemodel/ftFaceBase.h"
#include "GoTools/compositemodel/ftEdgeBase.h"
#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/compositemodel/ttlPoint.h"
#include "GoTools/compositemodel/ttlTriang.h"
#include "GoTools/compositemodel/ftSurfaceSet.h"
#include "GoTools/compositemodel/ftCurve.h"
#include "GoTools/compositemodel/cmUtils.h"

#include "GoTools/compositemodel/ftPlanarGraph.h"



namespace Go
{


  //===========================================================================
  /** ftChartSurface -  Short description.
   * Detailed description.
   *
   */
  //===========================================================================

  class ftChartSurface : public ftSurfaceSet
  {
  public:
  
    /// Constructor
    ftChartSurface(const std::vector<shared_ptr<ParamSurface> >& surfaces,
		   tpTolerances& topeps, double approxeps,
		   double curvature_tol, int m = 0, int n = 0);

    /// Destructor
    ~ftChartSurface();

    /// Warning! Assumes the return edges are taken care of, i.e. stored for
    /// further use. Otherwise first_edges_ will be invalid!
    /// To be called when topology is set up; essentially performed from constructor.
    virtual std::vector<shared_ptr<ftEdgeBase> > 
      createInitialEdges(double degenerate_epsilon = DEFAULT_SPACE_EPSILON,
			 double kink = 0.00015,
			 bool no_split = false);

    /// Evaluation and interrogation. Overridden from ftFaceBase.
    virtual Point point(double u, double v) const;
    /// Given par values in surf_, return corresponding par values in faces_[i].
    /// If use_input_face == true, we do not search in graph but go directly to face.
    Point point(double& u, double& v, shared_ptr<ftFaceBase>& face,
		double* seed = NULL, bool use_input_face = false) const;
    virtual Point normal(double u, double v) const;

    ///  Merge the current surface patches to produce one SplineSurface
    virtual ftMessage createSurf(double& max_error, double& mean_error);
    virtual void getError(double& max_error, double& mean_error);

    /// Priority type
    virtual ftTangPriority getPrioType() const;

    void setPrioType(ftTangPriority type);

    void setGridResolution(int m, int n);

    /// We may also decide to add grid pts along an edge.
    bool gridCreated(int& m, int& n) const;

    void setEdgeScales(std::vector<std::pair<Point, double> >& edge_scales);

    void setSymmDistrFunctions(bool symm_distr_functions);

    /// Flicking on grid_distr_functions_ so that they both agree along and towards common bd.
    /// For all members of topology it is expected that both createSurf() has been called.
    ftMessage prepareGrid(RotationInfo* rot_info = NULL, ftCurve* outer_bd = NULL);

    /// User may choose to let bd_edges of all sfs in the topology be in correspondence
    /// through a rotation of the entire set.
    /// For all members of topology it is expected that both createSurf() & prepareGrid()
    /// have been called.
    ftMessage createGrid(RotationInfo* rot_info = NULL, ftCurve* outer_bd = NULL);

    /// The created grid is returned (copied).
    /// The pts are returned in ccw order.
    /// Ordering is: b, t, l, r
    std::vector<shared_ptr<ftSamplePoint> > getEdgeGrid(int edge_ind);

    void writeGrid(std::ostream& os) const;

    /// Write to file LineCloud describing grid. Still awaiting exchange format.
    void writeDebugGrid(std::ostream& os) const;

    vector<shared_ptr<FaceConnectivity<ftEdgeBase> > > getInnerEdgeCont() const;
      
  protected:

    ftTangPriority prio_type_;
    //------------------ Data members ----------------------
    //    double top_tol_; // tpTolerances exist as class inherits from ftSurfaceSet.
    //    double approxtol_;  // The approximation surface should be relatively similar to
    // the subsurfaces. Large, applies only to corner points?
    double curvature_tol_; // If edge is too curved, we sample inner points of edge.
    ftPlanarGraph graph_; // As it refers to surf; created when surface is created.
    //?    bool sample_only_edges_
    double maxerror_;
    double meanerror_;

    // Surface is sampled in a m_*n_ grid (giving a (m_-1)*(n_-1) block structure).
    // m_ should be wrt the longest direction (i.e. avg of corr edges), n_ the shortest.
    // If m_ and n_ differ; as soon as surf_ is created values are swapped if needed.
    int m_;
    int n_;
    // For each edge we define a distribution function which remaps the parameter values along edge.
    std::vector<shared_ptr<SplineCurve> > grid_distr_functions_; // b, t, l, r
    // The pair consists of a space pt corresponding to an edge of the rectangular surface and a
    // scaling scalar used to override default distr near an edge.
    // As value was set prior to creation of sf_ we may not refer to a bd idx.
    std::vector<std::pair<Point, double> > edge_scales_;
    // If symm_distr_functions_ == true, initial distr f's in u(v) direction will coincide.
    bool symm_distr_functions_;
    std::vector<bool> secn_distr_;
    // grid_pts_ are stored in a m_*n_ matrix.
    std::vector<std::vector<shared_ptr<ftSurfaceSetPoint> > > grid_pts_; // m_*n_ pts, initialized in createGrid().

    /// Approximate only the edges in the group of surfaces
    virtual bool sampleOnlyEdges() const;

  private:

    // Create a 4-sided surface approximating the input surfaces, and
    // interpolating the boundary (if it is sufficiently smooth).
    // loop describes outer boundary of surface set, corner indicates corner edges.
    virtual ftMessage makeSurface(const std::vector<ftEdgeBase*>& loop,
				  std::vector<int>& corner,
				  double& max_error, double& mean_error);

    // For each edge, we approximate it's sampled points with a spline curve.
    // Basis is then returned, to be used for initial surface.
    void constructBasises(ftPointSet& pointset, const std::vector<int>& corners,
			  double epsge,
			  BsplineBasis& basis_u, BsplineBasis& v_basis);

    // We use only those points which are on a boundary (or subsurface-boundary).
    // Warning! Elements of points are really of type ftSurfaceSetPoint! We cast...
    ftMessage createGraph(std::vector<ftSamplePoint*>& graph_nodes); //ftPointSet& points);

    // Helper function returning number of points to evaluate in the inner of edge.
    virtual int nmbToEval(ftEdge* edge, double tmin, double tmax);

    // Following the boundary as given by first_point, we sample points between
    // of samples on boundary. We do not mark them as bnd points, or update any
    // neighbourhood structure. We sample nmb_per_int points between two bnd points.
    // Used when approximating the sf.
    // @@ Points are sampled on straight segments between boundary points. However,
    // bnd points were sampled so that each segment is toptol_.neighbour_-straight.
    void addOuterBoundaryPoints(ftPointSet& points);

    // Not to be called until surf_ is created.
    void createGridDistrFunctions();

    //     // Test routine to check performace of ftPlanarGraph::locateInGraph().
    //     // Writes output to screen.
    //     void testLocateInGraph(ftPointSet& points);

    // Assuming that underlying surface of face1 & face2 is a SplineSurface.
    ftMessage getMatchingEdges(Point from_space_pt, Point to_space_pt,
			       ftFaceBase* face1, ftFaceBase* face2,
			       int& edge1, double& par1, double& par2, 
			       int& edge2, double& par3, double& par4,
			       double epsgeo);

    shared_ptr<SplineCurve> 
      createOneDistrFunction(double par1[], double par2[], bool isudir);

    void modifyOneDistrFunction(shared_ptr<SplineCurve> func,
				double scale, bool atStart, int nbints);   

    // Modify the grid distribution functions to get smooth
    // transitions across block boundaries and correspondence
    // along block boundaries
    ftMessage modifyGridDistrFunctions(std::vector<BoundaryPiece>& bdpiece);

    ftMessage modifyGridDistrFunctions(RotationInfo* rot_info = NULL, ftCurve* outer_bd = NULL);

    // Make sure that the grid_distr_functions_ match along neighbouring faces.
    ftMessage updateGridDistribution(RotationInfo* rot_info = NULL, ftCurve* outer_bd = NULL);

    // Assuming that opposite indicates opposite directions for the two cvs.
    void adaptGridDistrFunc(SplineCurve *from_func,
			    SplineCurve *to_func, bool opposite);
    void averageAdjGridDistrFunc(SplineCurve *from_func,
				 SplineCurve *to_func, bool opposite);
    void averageGridDistrFunc(SplineCurve *func1, bool start1, double l1,
			      bool secn1, int nbints1, SplineCurve *func2, 
			      bool start2, double l2, bool secn2, int nbints2);
    // The general version of the prevoius function.
    void averageGridDistrFunc(std::vector<SplineCurve*>& funcs, std::vector<bool> start,
			      std::vector<double> lengths, std::vector<bool> secn,
			      const std::vector<int>& nbints);
    // Ordering of bdidx: (b, t, l, r)
    SplineCurve* getGridDistrFunc(int bdidx, double par1, double par2, double par_tol);

    SplineCurve* getPriorDistrFunc(int bdidx, bool opposite, bool& start,
				   double& length, bool& secn, int& nbints);

    void sampleGridPts();

    ftMessage getMatchingEdges(std::vector<ftEdgeBase*> local_outer_loop, std::vector<int> corners,
			       RotationInfo* rot_info, ftCurve* total_outer_loop,
			       std::vector<std::pair<int, int> >& matching_edges,
			       std::vector<ftChartSurface*>& matching_faces,
			       std::vector<double>& rot_angles);


    ftMessage getMatchingEdges(std::vector<ftEdgeBase*> outer_loop, std::vector<int> corners,
			       std::vector<std::pair<int, int> >& matching_edges,
			       std::vector<std::pair<int, int> >& matching_grid_res,
			       std::vector<ftChartSurface*>& matching_faces,
			       std::vector<double>& rot_angles);

    // To be called after sf_ has been created. Uses values in edge_scales_. Returns (b, t, l, r).
    std::vector<double> getEdgeScales();

  };


} // namespace Go


#endif // _FTCHARTSURFACE_COMP_H
 
