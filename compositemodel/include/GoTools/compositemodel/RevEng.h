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

#ifndef _REVENG_H
#define _REVENG_H

#include "GoTools/compositemodel/ftPointSet.h"
#include "GoTools/compositemodel/SurfaceModel.h"
#include "GoTools/compositemodel/RevEngRegion.h"
#include "GoTools/utils/Point.h"
#include "GoTools/utils/BoundingBox.h"
#include <string.h>

namespace Go
{
  class RevEngPoint;
  class HedgeSurface;
  class RevEngEdge;

  /// Elementary surface types to be recognized
  enum
  {
   PLANE, CYLINDER, SPHERE, CONE, TORUS
  };

  /// Characterization of model surface. Note that the default choice is ROUGH,
  /// and that other choices are currently not supported
  enum
  {
   SMOOTH=0, MEDIUM_ROUGH, ROUGH
  };

  /// Collection of surface information used in computation of global parameters
  struct SurfaceProperties
  {
    int sfix_;
    ClassType type_;
    ClassType prev_type_;
    Point dir_, loc_;
    double rad1_, rad2_;
    int num_points_;
    int surfflag_;

    SurfaceProperties(int ix, ClassType type, int num, int surfflag, Point& dir, 
		      Point& loc, ClassType prev_type=Class_Unknown, double rad1=-1,
		      double rad2=-1)
    {
      sfix_ = ix;
      type_ = type;
      num_points_ = num;
      surfflag_ = surfflag;
      dir_ = dir;
      loc_ = loc;
      prev_type_ = prev_type;
      rad1_ = rad1;
      rad2_ = rad2;
    }
  };

  /// Information computed in adaptToMainAxis. Specifies model axes including
  /// axis location and number of points supporting the information
  struct AxisInfo
  {
    Point axis_;
    std::vector<std::pair<Point,int> > plane_loc_;
    std::vector<std::pair<Point,int> > rotational_loc_;

    AxisInfo(Point& axis)
    {
      axis_ = axis;
    }

    void addPlaneLocation(Point& loc, int num)
    {
      plane_loc_.push_back(std::make_pair(loc,num));
    }
    
    void addRotationalLocation(Point& loc, int num)
    {
      rotational_loc_.push_back(std::make_pair(loc,num));
    }
  };
    
  /** RevEng -  Reverse engineering engine. Workflow funcationality and 
      storage of data.
   * 
   */
  
  class RevEng
  {
  public:
    /// Default constructor
    RevEng();
    
    /// Constructor
    /// \param tri_sf Initial triangulated surface
    RevEng(shared_ptr<ftPointSet> tri_sf);

    /// Destructor
    ~RevEng();

    /// Enhance points (triangle vertices) with estimated surface normal and
    /// curvature information
    void enhancePoints();

    /// Classify points with respect to estimated curvature information
    /// Possible classifications: edge point (high curvature), peak (negative mean
    /// curvature and positive Gauss), ridge (negative mean, zero Gauss),
    /// sridge (negative mean and Gauss), non (zero mean, positive Gauss, not expected),
    /// flat (zero mean and Gauss), minsurf (zero mean, negative Gauss),
    /// pit (positive mean and Gauss), valley (positive mean, zero Gauss)
    /// and svalley (positive mean, negative Gauss)
    void classifyPoints();

    /// Group connected points with the same classification into regions
    void segmentIntoRegions();

    /// Create first surfaces from large regions, simultanously dismissing
    /// points that do not belong to the region. Possible surfaces: planes,
    /// cylinders and cones
    void initialSurfaces();

    /// Add adjacent regions to regions with surfaces if appropriate
    void growSurfaces();

    /// Identify coordinate axes corresponding to the current model and
    /// update surfaces as appropriate (if their surface normal/axis almost
    /// corresponds to a coordinate axis and the approximation error stays
    /// limited)
    void updateAxesAndSurfaces();

    /// Compute intersection edges between defined neighbouring and almost
    /// neighbouring surfaces provided that these surfaces
    void firstEdges();

    /// Second attempt to create surfaces. Some more regions might have grown to a
    /// size that allows surface creation. Context information from adjacent surfaces
    /// are used to guide the surface recognition. In addition to planes, cylinders and
    /// cones are spheres and torii created. Context information is used to split
    /// large composite regions. Additional edges may be created.
    void surfaceCreation(int pass=1);

    /// Identify common surface axes/normals and define a local coordinate system.
    /// Update surfaces accordingly
    void updateRegionsAndSurfaces(int& ix, std::vector<RevEngRegion*>& grown_regions,
				  std::vector<HedgeSurface*>& adj_surfs);

    /// Second round in identifying main axes in the model. In addition to axis
    /// direction is location registered. Surfaces are updated accordingly.
    void adaptToMainAxis();

    /// Third attempt to create surfaces. Small regions with similar characteristica
    /// are merged and subject to surface recognition. The required region size is
    /// reduced. Possible surfaces: planes, cylinders and cones. Additional edges
    /// may be created.
    void smallRegionSurfaces();

    /// Create edge blends of type cylinder and torus associated with identified
    /// edges. Blend surfaces with almost similar radius are made consistent.
    void manageBlends1();

    /// Bound blend surface along the edge. Identify corner blends of type
    /// torus and 4-sided free form surface. Identify and extract trimming edges
    /// from blend surface and define associated trimming edges for adjacent surfaces
    void manageBlends2();

    /// Compute missing trimming edges, between surfaces and towards regions without
    /// surfaces. Extract relevant part of trimming edges, add missing edges as appropriate
    /// and combine edges into trimming loops. Define bounded surfaces and associated faces.
    void trimSurfaces();

    /// Perform topology analysis and create surface model that can be exported in a g22 file
    shared_ptr<SurfaceModel> createModel();

    /// Compute approximation tolerance based on model properties
    void setApproxTolerance();

    /// Set approximation tolerance from application
    void setApproxTol(double eps)
    {
      approx_tol_= eps;
    }
    
    /// Enquire defined approximation tolerance
    double getApproxTol()
    {
      return approx_tol_;
    }

    /// Enquire edge classification type. Current: always curvature
    int getEdgeClassificationType()
    {
      return edge_class_type_;
    }

    /// Enquire parameter for edge classification. Relates to estimated curvature
    /// radius in point and average length of triangle edges
    double getCfac()
    {
      return cfac_;
    }

    /// Set parameter for edge classification
    void setCfac(double cfac)
    {
      cfac_ = cfac;
    }

    /// Enquire point classification type. Current: always curvature
    int getClassificationType()
    {
      return classification_type_;
    }

    /// Enquire limit for when mean curvature is considered to be zero
    double getMeanCurvatureZero()
    {
      return zero_H_;
    }

    /// Set limit for when mean curvature is considered to be zero
    void setMeanCurvatureZero(double zero_H)
    {
      zero_H_ = zero_H;
    }

    
    /// Enquire limit for when Gauss curvature is considered to be zero
    double getGaussCurvatureZero()
    {
      return zero_K_;
    }

    /// Set limit for when Gauss curvature is considered to be zero
    void setGaussCurvatureZero(double zero_K)
    {
      zero_K_ = zero_K;
    }

    /// Enquire preference for primary surfaces relative to free form surfaces. Always
    /// primary, free form is currently disabled exept for some corner blends
    int getElementaryPreferLevel()
    {
      return prefer_elementary_;
    }

    /// Enquire model characterization involved in the definition of the approximation tolerance.
    /// Currently: always ROUGH
    int getModelCharacterization()
    {
      return model_character_;
    }

    /// The main axes of the model is computed by updateAxesAndSurfaces(). Prelimenary axes
    /// my be set prior to this function, but they should not be altered afterwards
   void setMainAxis(Point mainaxis[3])
    {
      mainaxis_[0] = mainaxis[0];
      mainaxis_[1] = mainaxis[1];
      mainaxis_[2] = mainaxis[2];
    }

    /// Enquire local coordinate axes
    void getMainAxis(Point mainaxis[3])
    {
      mainaxis[0] = mainaxis_[0];
      mainaxis[1] = mainaxis_[1];
      mainaxis[2] = mainaxis_[2];
    }

    /// Enquire current number of regions
    int numRegions()
    {
      return (int)regions_.size();
    }

    /// Access specified region
    shared_ptr<RevEngRegion> getRegion(int ix)
    {
      return regions_[ix];
    }

    /// Enquire number of identified surfaces
    int numSurfaces()
    {
      return (int)surfaces_.size();
    }

    /// Access specified surface
    shared_ptr<HedgeSurface> getSurface(int ix)
    {
      return surfaces_[ix];
    }

    /// Enquire minimum number of points in a region for being considered for surface
    /// recognition
    double getMinPointRegion()
    {
      return min_point_region_;
    }

    /// Store current stage for later to restart the computation from this stage (not
    /// accessible from manageBlends2 and later)
    void storeGrownRegions(std::ostream& os);

    /// Read stored stage
    void readGrownRegions(std::istream& is);

    /// Debug functionality
    /// Write point groups and associated surfaces to files depending on the
    /// the number of points in the groups (uses min_point_region_ to distinguish with
    /// respect to size)
    void writeRegionStage(std::ostream& of, std::ostream& ofm, std::ostream& ofs) const;

    /// Write point regions with an assoicated surface and the surface to file
    void writeRegionWithSurf(std::ostream& of) const;

    /// Write edges and point regions appropriated to an associated blend surface to file
    void writeEdgeStage(std::ostream& of) const;
    
  private:
    /// Characterizes the surface of the current model. Currently always set to ROUGH
    int model_character_;

    /// Triangulated surface
    shared_ptr<ftPointSet> tri_sf_;

    /// Average length of triangle edges
    double mean_edge_len_;

    /// Group points assumed to be associated one surface
    std::vector<shared_ptr<RevEngRegion> > regions_;

    /// Points with different properties than their neighbouring points.
    /// Regions with one point can also occur
    std::vector<RevEngPoint*> single_points_;

    /// Recognized surfaces
    std::vector<shared_ptr<HedgeSurface> > surfaces_;

    /// Recognized edges. Intersection curves between surfaces with additional information
    std::vector<shared_ptr<RevEngEdge> > edges_;

    /// Final surface model including topology
    shared_ptr<SurfaceModel> sfmodel_;

    /// Bounds the point cloud (triangulated surface)
    BoundingBox bbox_;
    
    /// Minimum number of neighbouring connected points for estimating point properties
    /// (surface normal and curvature)
    int min_next_;
    
    // Stop searching for neighbouring points to estimate point properties when this number
    // is reached
    int max_next_;
    
    /// Factor for radius in which to search for neighbouring points
    double rfac_;

    int edge_class_type_ = CURVATURE_EDGE;
    int classification_type_ = CLASSIFICATION_CURVATURE;
    
    /// A vertex is considered an edge points if the estimated curvature radius is 
    /// less than cfac_ times the average length of triangulation edges in the vertex
    double cfac_;
    
    /// Limit for when the cone angle corresponding to triangle normals indicate an edge
    /// Not active
    double norm_ang_lim_;

    /// Limit for when the cone angle correspondingto triangle normals indicate a plane
    /// Not active
    double norm_plane_lim_;  

    /// When mean curvature is considered zero
    double zero_H_;  

    /// When Gauss curvature is considered zero
    double zero_K_;

    /// Minimum number of points in a region to be considered for surface recognition
    int min_point_region_;

    /// Approximation tolerance
    double approx_tol_; 

    /// Intersection tolerance
    double int_tol_;

    /// Angular tolerance used in topology build. 5.0*anglim_ is used to check whether
    /// the estimated normal in a point corresponds to the normal of a recognized surface
    /// in the point
    double anglim_;

    /// Used to decide if a point is an outlier
    int max_nmb_outlier_;

    /// Preferance in selecting the surface associated to a point cloud, primary surface or
    /// best fit surface. Default is 0. 0 = always, 1 = preferred, 2 = best accuracy
    int prefer_elementary_; 

    /// Coordinated axes associated to the model
    Point mainaxis_[3];

    /// Identified model axes including position
    std::vector<AxisInfo> model_axis_;

    /// Collection of information related to smallRegionSurfaces()
    struct SmallSurface
    {
      int axis_ix_, pos_ix_, lev_ix_;
      vector<vector<RevEngPoint*> > assos_points_;
      vector<shared_ptr<ElementarySurface> > surfs_;
      vector<BoundingBox> bbox_;
      int type_;   // 1 = plane1, 2=plane2, 3=rotational, 4=from remaining

      SmallSurface(int ix1, int ix2, int ix3, int type,
		   vector<shared_ptr<ElementarySurface> >& surfs)
      {
	axis_ix_ = ix1;
	pos_ix_ = ix2;
	lev_ix_ = ix3;
	type_ = type;
	surfs_ = surfs;
      }

      void addPoints(vector<RevEngPoint*>& points, BoundingBox bb)
      {
	assos_points_.push_back(points);
	bbox_.push_back(bb);
      }
    };
    
    /// Set edge classification type. 
    void setEdgeClassificationType(int edge_class_type)
    {
      edge_class_type_ = edge_class_type;
    }
    
    /// Set point classification type
    void setClassificationType(int classification_type)
    {
      classification_type_ = classification_type;
    }

    void setElementaryPreferLevel(int preferlevel)
    {
      prefer_elementary_ = preferlevel;
    }
    
    void setModelCharacterization(int character)
    {
      //model_character_ = std::min(ROUGH, std::max(SMOOTH, character));
      model_character_ = std::min(2, std::max(0, character));
    }

     double getInitApproxTol();
    void edgeClassification();
    void curvatureFilter();

    void initParameters();
    void updateParameters();
    bool recognizeOneSurface(int& ix, int min_point_in, double angtol,
			     int pass);
    void recognizeSurfaces(int min_point_in, int pass);
    void recognizeEdges(bool only_curve=false);
    bool createBlendSurface(int ix);
    void adjustPointRegions(int min_point_in);

    void computeAxisFromCylinder(Point initaxis[3], int min_num, double max_ang,
				 Point axis[3], int num_points[3]);
    
    void computeAxisFromPlane(Point initaxis[3], int min_num, double max_ang,
			      Point axis[3], int num_points[3]);
    
    void surfaceExtractOutput(int idx,
			      std::vector<std::vector<RevEngPoint*> > out_groups,
			      std::vector<HedgeSurface*> prev_surfs);

    bool segmentComposite(int& ix, int min_point_in, double angtol);
    bool segmentByPlaneGrow(int ix, int min_point_in, double angtol);
    bool segmentByAxis(int ix, int min_point_in);
    bool segmentByContext(int ix, int min_point_in, double angtol, bool first);
    void growSurface(int& ix, int pass = 1);
    void growBlendSurface(int& ix);
    void growMasterSurface(int& ix);

    void defineSmallRegionSurfaces();
    
    void growSmallRegionSurface(int& ix);

    bool identifySmallRotational(std::vector<RevEngPoint*>& points,
				 Point midp, Point loc, Point axis, Point Cx,
				 double ppar1, double ppar2,
				 std::vector<shared_ptr<ElementarySurface> >& sfs);

    bool identifySmallPlanar(std::vector<RevEngRegion*>& groups,
			     Point loc, Point axis, Point Cx,
			     double ppar1, double ppar2, double delta,
			     std::vector<shared_ptr<ElementarySurface> >& sfs);

    bool identifySmallPlanar(std::vector<RevEngPoint*>& groups,
			     Point loc, Point axis, Point Cx,
			     double ppar1, double ppar2, double delta,
			     std::vector<shared_ptr<ElementarySurface> >& sfs);
    
    shared_ptr<ElementarySurface>
    defineElemSurf(std::vector<RevEngPoint*>& points, std::vector<RevEngPoint*>& in_points,
		   BoundingBox& bbox, std::vector<RevEngPoint*>& remain);

    void planarAtPlane(shared_ptr<Plane> axis_plane,
		       std::vector<RevEngPoint*>& points,
		       std::vector<HedgeSurface*>& sfs,
		       std::vector<shared_ptr<RevEngRegion> >& plane_sf_reg,
		       std::vector<shared_ptr<HedgeSurface> >& plane_sf_hedge);
    
    void integrateInSmallSurfs(std::vector<shared_ptr<RevEngRegion> >& small_sf_reg,
			       std::vector<RevEngRegion*>& nosf_reg,
			       std::vector<RevEngRegion*>& include_reg);
    
    void extractSmallSurfs(SmallSurface& small_surf,
			   std::vector<shared_ptr<RevEngRegion> >& small_sf_reg,
			   std::vector<shared_ptr<HedgeSurface> >& small_sf_hedge,
			   std::vector<RevEngRegion*>& nosf_reg,
			   std::vector<RevEngPoint*>& non_assigned_pts);

    void doAdaptToAxis();

    bool axisUpdate(int ix, double max_ang, double angtol);
    
    void adjustWithMainAxis(std::vector<Point>& axes, std::vector<int>& num_pts);
    Point planarFit(std::vector<int>& sf_ix, Point axis);

    Point rotationalFit(std::vector<int>& sf_ix, Point axis, Point Cx,
			std::vector<RevEngEdge*>& nopar_edgs);

    void collectAxis(std::vector<SurfaceProperties>& sfprop);

    void computeLocFunc(RevEngPoint* pt, std::vector<RevEngPoint*>& points,
			Point& vec1, Point& vec2, double radius);

    int setSmallRegionNumber();
    
    void storeParams(std::ostream& os) const;
    void readParams(std::istream& is);
    void setBoundingBox();
    
    void checkConsistence(std::string text) const;

    std::vector<shared_ptr<RevEngEdge> >
    defineEdgesBetween(size_t ix1,shared_ptr<ElementarySurface>& surf1,
		       Point& dir1, size_t ix2, shared_ptr<ElementarySurface>& surf2,
		       Point& dir2, bool only_curve=false,
		       double lenlim=-1.0, bool check_common = true);
    
    shared_ptr<RevEngEdge> 
    defineOneEdge(size_t ix1, shared_ptr<ElementarySurface>& surf1,
		  Point& dir1, size_t ix2,
		  shared_ptr<ElementarySurface>& surf2, Point& dir2,
		  shared_ptr<CurveOnSurface>& int_cv1,
		  shared_ptr<CurveOnSurface>& int_cv2,
		  double width, std::vector<RevEngRegion*>& common_reg,
		  bool only_curve, bool check_common);
    
    RevEngPoint* getDistantPoint(shared_ptr<CurveOnSurface>& cv,
				 double tmin, double tmax, double dd,
				 double width,
				 std::vector<RevEngPoint*>& points);

    void extendBlendAssociation(size_t ix);
    
    bool setBlendEdge(size_t ix);
    
    double
    computeCylinderRadius(std::vector<std::vector<RevEngPoint*> > blend_pts,
			  double width, const Point& pos, const Point& axis,
			  const Point& dir1, const Point& dir2);
   shared_ptr<Cylinder>
    createCylinderBlend(std::vector<std::vector<RevEngPoint*> > blend_pts,
			double rad1, const Point& pos, const Point& axis,
			const Point& dir1, const Point& dir2, int sgn);
    
    double
    computeTorusRadius(std::vector<std::vector<RevEngPoint*> >& blend_pts,
		       shared_ptr<CurveOnSurface>& cv,
		       const Point& locp, const Point& normal,
		       shared_ptr<ElementarySurface> rotational,
		       double width, bool plane_out, bool rot_out);
    
    void getTorusParameters(shared_ptr<ElementarySurface> planar,
			    shared_ptr<ElementarySurface> rotational,
			    double radius, int sgn1, int sgn2, double& Rrad, 
			    Point& centre, Point& normal, Point& Cx);
    
    double
    computeTorusRadius(std::vector<std::vector<RevEngPoint*> >& blend_pts,
		       shared_ptr<CurveOnSurface>& cv,
		       shared_ptr<ElementarySurface> elem1,
		       shared_ptr<ElementarySurface> elem2,
		       double width, bool out1, bool out2, int sgn,
		       double& d2);
    
  bool
    getTorusParameters(shared_ptr<ElementarySurface> elem1,
		       shared_ptr<ElementarySurface> elem2, Point pos,
		       double radius, double d2, bool out1, bool out2, int sgn,
		       double& Rrad, Point& centre, Point& normal, Point& Cx,
		       bool check_common = true);

  shared_ptr<Torus>
    torusBlend(std::vector<std::vector<RevEngPoint*> >& blend_pts,
	       std::vector<shared_ptr<CurveOnSurface> >& cvs,
	       const Point& locp, const Point& normal,
	       shared_ptr<ElementarySurface> rotational,
	       double width, bool plane_out, bool rot_out);
    
    shared_ptr<Torus>
    torusBlend(std::vector<std::vector<RevEngPoint*> >& blend_pts,
	       shared_ptr<CurveOnSurface>& cv,
	       shared_ptr<ElementarySurface> elem1,
	       shared_ptr<ElementarySurface> elem2,
	       double width, bool out1, bool out2, int sgn);
    
    void setBlendBoundaries(RevEngRegion *reg);

    void equalizeBlendRadii();

    void equalizeAdjacent(size_t ix1, size_t ix2);
    
    void updateBlendRadius(size_t ik, double radius);
    
    bool defineTorusCorner(size_t ix);

    void defineMissingCorner(std::vector<RevEngRegion*>& cand_adj);

    bool createTorusBlend(size_t ix);

    bool suitcaseCorner(std::vector<RevEngRegion*>& adj_blends,
			RevEngEdge *rev_edge);

    void extractOutPoints(std::vector<RevEngPoint*>& points, shared_ptr<ParamSurface> surf,
			  std::vector<int>& cv_ix,
			  double tol, double angtol,
			  std::vector<std::vector<RevEngPoint*> >& move2adj,
			  std::vector<RevEngPoint*>& remain);

    void storeClassified(std::ostream& os) const;
    void readClassified(std::istream& is);

  };

} // namespace Go

#endif // _REVENG_H
