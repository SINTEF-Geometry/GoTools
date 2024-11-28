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

#ifndef _REVENGREGION_H
#define _REVENGREGION_H

#include "GoTools/compositemodel/RevEngPoint.h"
#include "GoTools/compositemodel/ImplicitApprox.h"
#include "GoTools/utils/BoundingBox.h"
#include "GoTools/geometry/ClassType.h"
//#include "GoTools/compositemodel/HedgeSurface.h"
#include <set>


namespace Go
{
  class HedgeSurface;
  class ftEdge;
  //class RevEngPoint;
  class Circle;
  class SplineCurve;
  class CurveOnSurface;
  class Plane;
  class Cylinder;
  class Sphere;
  class Cone;
  class Torus;
  class SplneSurface;
  class RevEngEdge;
  
  /// Method to classify points (vertices). Default and currently the only choice: CLASSIFICATION_CURVATURE
  enum
    {
     CLASSIFICATION_UNDEF, CLASSIFICATION_CURVATURE, CLASSIFICATION_SHAPEINDEX, CLASSIFICATION_POINTASSOCIATION
    };

  /// Method used in edge classification. Default and currently the only choice: CURVATURE_EDGE
  enum
    {
     TRIANGULATION_EDGE, PCATYPE_EDGE, CURVATURE_EDGE, CNESS_EDGE, RPFAC_EDGE
    };

  /// Preference for elementary. Default and recommended choice: ALWAYS_ELEM
  enum
    {
     ALWAYS_ELEM, PREFER_ELEM, BEST_ACCURACY
    };

  /// Surface flag. Quality of approximation.
  /// ACCURACY_OK: Requirements met for distance between points and surface and angle between normals
  /// ANGULAR_DEVIATION: Requirements met for distance between points and surface 
  /// PROBABLE_HELIX: As for ANGULAR_DEVIATION in the cylinder and cone case
  /// ACCURACY_POOR: Approximation accuracy not met, but the surface is relatively close to the points
  /// NOT_SET: No surface or large distance
  enum
    {
     ACCURACY_OK, ANGULAR_DEVIATION, PROBABLE_HELIX, ACCURACY_POOR, NOT_SET
    };

  /// Surface history. Whether or not the surface is modified after first computation
  enum
    {
     INITIAL, AXIS_ADAPTED, ADJACENT_ADAPTED
    };

  /// Nessecary information to create a swept spline surface. Currently disabled
  struct SweepData
  {
    int type_;  // Linear = 1, rotational = 2, cylinderlike = 3, conelike = 4
    shared_ptr<SplineCurve> profile_;
    Point location_;
    Point added_info_;
    double radius_;
    double angle_;
    double maxdist_;
    double avdist_;
    int num_in_;

    SweepData(int type, shared_ptr<SplineCurve> profile, Point location, 
	      Point info2, double maxdist, double avdist, int num_in,
	      double radius = -1.0, double angle = 0.0)
    {
      type_ = type;
      profile_ = profile;
      location_ = location;
      added_info_ = info2;
      maxdist_ = maxdist;
      avdist_ = avdist;
      num_in_ = num_in;
      radius_ = radius;
      angle_ = angle;
    }
  };

  struct SegmentData
  {
    int type_; // Around axis = 1
    Point loc_;
    Point axis_;
    double min_dist_;
    double max_dist_;

    SegmentData(int type, Point& loc, Point& axis, double mind, double maxd)
    {
      type_ = type;
      loc_ = loc;
      axis_ = axis;
      min_dist_ = mind;
      max_dist_ = maxd;
    }
  };
  
  /** RevEngRegion - Point groups after segmenting a triangulated surface in regions believed
      to be appropriate for representaton by one CAD type surface. The regions are dynamic, points
      are removed and added as appropriate. Class entities are primarily accessed by RevEng, but
      also from RevEngEdge, RevEngPoint and HedgeSurface
   * 
   */

  class RevEngRegion
  {
  public:
    /// Constructor
    RevEngRegion(int edge_class_type);

    RevEngRegion(int classification_type, int edge_class_type);

    /// Constructor given a group of RevEngPoints
    RevEngRegion(int classification_type, int edge_class_type,
		 std::vector<RevEngPoint*>& points);

    /// Destructor
    ~RevEngRegion();

    /// Set entity identity. Used in storing of stage (see RevEng)
    void setId(int Id)
    {
      Id_ = Id;
    }

    /// Enquire entity identity
    int getId()
    {
      return Id_;
    }

    /// Enquire type
     int getClassificationType()
    {
      return classification_type_;
    }

     int getEdgeClassificationType()
    {
      return edge_class_type_;
    }

    /// Check if the group is compatible with a surface type. NB! Not stable information
    bool isCompatible(ClassType classtype, int sfcode);
    
    /// Extend group with one point
    void addPoint(RevEngPoint* point);

    /// Extend group with points.
    /// NB! The group must be associated a surface
    /// NB! No testing on whether the points actually belongs to this group.
    /// Accuracy statistics is computed
    bool addPointsToGroup(std::vector<RevEngPoint*>& points,
			  double tol, double angtol, bool compute_accuracy = true);
    
    /// Remove one point. NB! Leaves overview information invalid.
    void removePoint(RevEngPoint* point);

    /// Extend current region with the points of another region. Update class
    /// parameters
    /// \param dist_ang If given: update new points with distance to surface and angular
    /// deviation between normals
    void addRegion(RevEngRegion* reg,
		   std::vector<std::pair<double, double> >& dist_ang,
		   double maxd=0.0, double avd=0.0, int num_inside=-1,
		   int num_inside2=-1);
    
    /// Update class parameters dependent on point properties
    void updateInfo(double tol=-1.0, double angtol=-1.0);

    /// Extend region with adjacent points having the same classification
    void collect(RevEngPoint *pt, RevEngRegion* prev=0);

    /// Enquire number of points in region
    int numPoints()
    {
      return (int)group_points_.size();
    }

    /// Get specified point
    RevEngPoint* getPoint(int ix)
    {
      if (ix < 0 || ix >= (int)group_points_.size())
	return 0;
      else
	return group_points_[ix];
    }

    /// Iterator to start of points in region
    std::vector<RevEngPoint*>::iterator pointsBegin()
    {
      return group_points_.begin();
    }
    
   /// Iterator past points in region
    std::vector<RevEngPoint*>::iterator pointsEnd()
    {
      return group_points_.end();
    }

    /// Return all points
    std::vector<RevEngPoint*> getPoints()
    {
      return group_points_;
    }

    /// Return bounding box surrounding all points
    const BoundingBox& getBbox()
    {
      return bbox_;
    }

    /// Return bounding box surrounding the parameter points corresponding to the
    /// points. Only if the region is associated a surface
    BoundingBox getParameterBox();

    /// Enquire the fraction of LocFunc normals that is closer to the average LocFunc normal
    /// in the region than the angular tolerance (used to identify possible planar regions)
    double getFracNorm()
    {
      return frac_norm_in_;
    }
    
   /// Extract a sub group of points that can be approximated by a plane and compute the plane
    /// \param mainaxis Local coordinate system for model
    /// \param tol Approximation tolerance
    /// \param plane_pts Extracted points
    /// \param plane_out Computed plane (if a sufficient number of points are found)
    void growLocalPlane(Point mainaxis[3], double tol,
			std::vector<RevEngPoint*>& plane_pts,
			shared_ptr<Plane>& plane_out);

    /// If the region is associated a plane, a cylinder or a cone, extend the point group
    /// with points from adjacent regions if feasible. It is ensured that all affected
    /// regions are connected
    /// \param mainaxis Local coordinate system for model
    /// \param tol Approximation tolerance
    /// \param grown_regions Regions completely absorbed in the current
    /// \param adj_surfs Surfaces associated with the grown regions
    /// \param added_groups Connected groups of points extracted from adjacent regions
    /// in order to maintain connectivity for these regions
    void growPlaneOrCyl(Point mainaxis[3], int min_pt_reg,
			double tol, double angtol,
			std::vector<RevEngRegion*>& grown_regions,
			std::vector<HedgeSurface*>& adj_surfs,
			std::vector<std::vector<RevEngPoint*> >& added_groups);

    /// Extract points with a significant distance to the surface and a large angular deviation
    /// \param added_groups Connected sub groups of removed points
    void removeLowAccuracyPoints(int min_pt_reg, double tol, double angtol,
				 std::vector<std::vector<RevEngPoint*> >& added_groups);

    /// Include adjacent region is current provided that current is associated with
    /// a surface and that the accuracy is sufficient. Accuracy towards surface is computed
    /// \param adj Region to include
    /// \param mainaxis Local coordinate system for model
    /// \param tol Approximation tolerance
    /// \param angtol Angular tolerance
    /// \param grown_regions Region absorbed in the current
    /// \param adj_surfs Surface associated with the grown region
    /// \param return Whether or not the region is included
    bool includeAdjacent(RevEngRegion* adj, Point mainaxis[3], 
			 double tol, double angtol,
			 std::vector<RevEngRegion*>& grown_regions,
			 std::vector<HedgeSurface*>& adj_surfs);

    /// Grow current region with adjacent regions. Current region must be associated
    /// with a surface and the accuracy must be sufficient
    /// \param mainaxis Local coordinate system for model
    /// \param min_pt_reg Associated with surface flag. Not used
    /// \param tol Approximation tolerance
    /// \param angtol Angular tolerance
    /// \param grown_regions Regions absorbed in the current. To be removed in RevEng
    /// \param adj_surfs Surfaces associated with the grown regions
    /// \param adj_edgs Surfaces associated with the grown regions
    void growWithSurf(Point mainaxis[3], int min_pt_reg,
		      double tol, double angtol,
		      std::vector<RevEngRegion*>& grown_regions,
		      std::vector<HedgeSurface*>& adj_surfs,
		      std::vector<RevEngEdge*>& adj_edgs,
		      bool use_base=false);

    /// Grow current region associated with a blend surface with adjacent regions not associated
    /// any surface. The accuracy must be sufficient
    void growBlendSurf(std::vector<RevEngRegion*>& next_blend, double tol,
		       double angtol, std::vector<RevEngRegion*>& grown_regions,
		       std::vector<std::vector<RevEngPoint*> >& added_regions);

    /// Check potential accuracy of integrating the region other into current.
    /// Current must be associated a surface
    int getGrowAccuracy(RevEngRegion *other, double tol,
			double angtol, double& maxdist,
			double& avdist, int& num_in, int& num2_in,
			double& maxdist2, double& avdist2, int& num_in2,
			int& num2_in2, std::vector<double>& parvals,
			std::vector<std::pair<double,double> >& distang);

    /// If feasible, merge adjacent regions without an associated surface into current.
    /// Current is not expected to have an associated surface. The combined regions must
    /// be possible to approximate with a plane
    bool mergePlanarReg(double zero_H, double zero_K, double tol,
			Point mainaxis[3],
			std::vector<RevEngRegion*>& grown_regions);

    /// Include adjacent regions with the same type of associated surface into current
    /// if the accuracy is sufficient. The current surface may be updated
    void mergeAdjacentSimilar(double tol, double angtol,
			      std::vector<RevEngRegion*>& grown_regions,
			      std::vector<HedgeSurface*>& adj_surfs,
			      std::vector<RevEngEdge*>& adj_edgs);

    /// Try to approximate the region points with a plane, a cylinder or a cone.
    /// If a surface is recognized will deviant points be extracted
    void
    initPlaneCyl(int min_point, int min_pt_reg,
		 double tol, double angtol, Point mainaxis[3],
		 double zero_H, double zero_K,
		 std::vector<shared_ptr<HedgeSurface> >& hedgesfs,
		 std::vector<vector<RevEngPoint*> >& out_groups,
		 std::vector<RevEngPoint*>& single_pts, bool& repeat);
    
    /// Not connected parts of the region is extracted and split into
    /// connected sub groups
    void splitRegion(std::vector<std::vector<RevEngPoint*> >& separate_groups);

    /// Split region with planar points according to the direction of the point normals.
    /// To divide the region into groups feasible for being represented by planes with
    /// varying normal
    void splitPlanar(double lim_cone, int min_point_reg, 
		     std::vector<std::vector<RevEngPoint*> >& other_groups,
		     std::vector<RevEngPoint*>& single);
    
    /// Include adjacent regions into current if the accuracy is sufficient. The current surface
    /// is not updated
    void joinToCurrent(double tol, double angtol, int small_lim,
		       std::vector<RevEngRegion*>& adapted_regions);

    /// Merge adjacent regions if consistent. The regions are not associated a surface.
    /// Time consuming and currently not in use
    void joinRegions(Point mainaxis[3], double approx_tol, double anglim,
		     std::vector<RevEngRegion*>& adapted_regions);

    /// Extract points with deviant accuracy. The removed points are divided into
    /// connected groups. The current region must be associated a surface
    void extractOutPoints(std::vector<std::pair<double, double> >& dist_ang,
			  double tol, double angtol,
			  std::vector<std::vector<RevEngPoint*> >& out_groups);

    /// Remove specified points from the current region. Information related to the
    /// region in the points is unset
    void removeAndUpdatePoints(vector<RevEngPoint*>& points);

    /// Identify points associated to a blend surface corresponding to an edge
    /// between this region and another region. 
    void extractOutOfEdge(shared_ptr<CurveOnSurface>& cv,
			  std::vector<shared_ptr<CurveOnSurface> >& intcv,
			  double radius, double tol, double angtol,
			  std::vector<RevEngPoint*>& out_points);

    /// Identify points on the "other side" of an edge between this region and another
    /// region. 
    void extractOutOfEdge2(std::vector<shared_ptr<CurveOnSurface> >& intcv,
			   double tol, double angtol,
			  std::vector<RevEngPoint*>& out_points);

    /// Identify points with deviant accuracy. The limit for being deviant is
    /// set dependent on tol and avd
    void identifyDistPoints(std::vector<std::pair<double, double> >& dist_ang,
			    double tol, double maxd, double avd,
			    std::vector<RevEngPoint*>& dist_points);

    /// Check if the region points are connected trough edges
    bool isConnected();

    /// Split into connected groups
    /// \param move A sub set of the region points to check for connectivity.
    /// \param out_groups Connected sub groups of move
    /// \param outer Whether or not only sub groups of move at the boundary of this
    /// regions is to be added to out_groups
    /// \param inner The remaining points in move if outer = true
   void connectedGroups(std::vector<RevEngPoint*>& move,
			std::vector<std::vector<RevEngPoint*> >& out_groups,
			bool outer, std::vector<RevEngPoint*>& inner);

    /// Extract specified points from the region. The specified points are
    /// divided into connected groups. If outer=true, points from move internal to
    /// the regions are kept
   void extractSpesPoints(std::vector<RevEngPoint*>& move,
			   std::vector<std::vector<RevEngPoint*> >& out_groups,
			   bool outer=false);

    /// Extract points with deviant accuracy. The removed points are divided into
    /// connected groups. The current region must be associated a surface
    void extractOutPoints(int dir, double tol, double angtol,
			  double angtol2,
			  std::vector<std::vector<RevEngPoint*> >& out_groups);
    
    /// Remove specified points from the current region. Points are not updated
    void removePoints(std::vector<RevEngPoint*>& remove);

    /// Enquire region classification. Unstable information
   int getClassification()
    {
      if (classification_type_ == CLASSIFICATION_CURVATURE)
	return group_points_[0]->C1_surf();
      else
	return -1;
    }

    /// Check if the regions is classified for being represented by a cylinder. Unstable information
   bool cylindertype()
    {
      if (classification_type_ == CLASSIFICATION_CURVATURE)
	return (group_points_[0]->C1_surf() == C1_RIDGE ||
		group_points_[0]->C1_surf() == C1_VALLEY);
      else
	return false;
    }
    
    /// Check if the regions is classified for being represented by a plane. Unstable information
   bool planartype()
    {
      if (classification_type_ == CLASSIFICATION_CURVATURE)
	return (group_points_[0]->C1_surf() == C1_FLAT);
      else
	return false;
    }

    /// Approximate region points with a plane if feasible, and register the plane
    /// in the region. Possibly extract devinant points in connected groups
    /// \param  mainaxis Local coordinate system for model
    /// \param tol Approximation tolerance
    /// \param min_pt Not used
    /// \param min_pt_reg Not used
    /// \param angtol Angular tolerance
    /// \param prefer_elementary Used in checking if an already existing surface should be replaced
    /// \param hedgesfs New surface
    /// \param prevsfs Possibly previous surface being replaced
    /// \param out_groups Deviant points being removed from region
    bool extractPlane(Point mainaxis[3],
		      double tol, int min_pt, int min_pt_reg, double angtol,
		      int prefer_elementary,
		      std::vector<shared_ptr<HedgeSurface> >& hedgesfs,
		      std::vector<HedgeSurface*>& prevsfs,
		      std::vector<std::vector<RevEngPoint*> >& out_groups);

    /// Approximate region points with a cylinder if feasible, and register the cylinder
    /// in the region. Possibly extract devinant points in connected groups
    /// \param tol Approximation tolerance
    /// \param min_pt Minimum number of points left after extracting deviant points to create
    /// a surface
    /// \param min_pt_reg Not used
    /// \param angtol Angular tolerance
    /// \param prefer_elementary Used in checking if an already existing surface should be replaced
    /// \param hedgesfs New surface
    /// \param prevsfs Possibly previous surface being replaced
    /// \param out_groups Deviant points being removed from region
    /// \param repeat Indicates that points are removed and another try to fit a cylinder is recommended
    bool extractCylinder(double tol, int min_pt, int min_pt_reg, 
			 double angtol, int prefer_elementary,
			 std::vector<shared_ptr<HedgeSurface> >& hedgesfs,
			 std::vector<HedgeSurface*>& prevsfs,
			 std::vector<std::vector<RevEngPoint*> >& out_groups,
			 bool& repeat);

    /// Check if the regions is suitable for being represented by a plane. Unstable information
    bool feasiblePlane(double zero_H, double zero_K) const;

    /// Check if the regions is suitable for being represented by a cylinder. Unstable information
    bool feasibleCylinder(double zero_H, double zero_K) const;

    /// Approximate region points with a sphere if feasible, and register the sphere
    /// in the region. Possibly extract devinant points in connected groups
    /// \param  mainaxis Local coordinate system for model
    /// \param tol Approximation tolerance
    /// \param min_pt Not used
    /// \param min_pt_reg Not used
    /// \param angtol Angular tolerance
    /// \param prefer_elementary Used in checking if an already existing surface should be replaced
    /// \param hedgesfs New surface
    /// \param prevsfs Possibly previous surface being replaced
    /// \param out_groups Deviant points being removed from region
     bool extractSphere(Point mainaxis[3],
		      double tol, int min_pt, int min_pt_reg, 
		       double angtol, int prefer_elementary,
		       std::vector<shared_ptr<HedgeSurface> >& hedgesfs,
		       std::vector<HedgeSurface*>& prevsfs,
		       std::vector<std::vector<RevEngPoint*> >& out_groups);

    /// Approximate region points with a linear swept spline surface if feasible, and register 
    /// the surface in the region. NB! Currently not in use
    /// \param  mainaxis Local coordinate system for model
    /// \param tol Approximation tolerance
    /// \param min_pt Not used
    /// \param min_pt_reg Not used
    /// \param angtol Angular tolerance
    /// \param prefer_elementary Used in checking if an already existing surface should be replaced
    /// \param hedgesfs New surface
    /// \param prevsfs Possibly previous surface being replaced
    bool extractLinearSweep(double tol, int min_pt, int min_pt_reg, 
			    double angtol, int prefer_elementary,
			    std::vector<shared_ptr<HedgeSurface> >& hedgesfs,
			    std::vector<HedgeSurface*>& prevsfs);

    /// Approximate region points with a cone if feasible, and register the plane
    /// in the region. Possibly extract devinant points in connected groups
    /// \param tol Approximation tolerance
    /// \param min_pt Not used
    /// \param min_pt_reg Not used
    /// \param angtol Angular tolerance
    /// \param prefer_elementary Used in checking if an already existing surface should be replaced
    /// \param hedgesfs New surface
    /// \param prevsfs Possibly previous surface being replaced
    /// \param out_groups Deviant points being removed from region
    bool extractCone(double tol, int min_pt, int min_pt_reg, 
		     double angtol, int prefer_elementary,
		     std::vector<shared_ptr<HedgeSurface> >& hedgesfs,
		     std::vector<HedgeSurface*>& prevsfs,
		     std::vector<std::vector<RevEngPoint*> >& out_groups);

    /// Approximate region points with a torus if feasible, and register the torus
    /// in the region. Possibly extract devinant points in connected groups
    /// \param  mainaxis Local coordinate system for model
    /// \param tol Approximation tolerance
    /// \param min_pt Not used
    /// \param min_pt_reg Not used
    /// \param angtol Angular tolerance
    /// \param prefer_elementary Used in checking if an already existing surface should be replaced
    /// \param hedgesfs New surface
    /// \param prevsfs Possibly previous surface being replaced
    /// \param out_groups Deviant points being removed from region
   bool extractTorus(Point mainaxis[3],
		      double tol, int min_pt, int min_pt_reg, 
		      double angtol, int prefer_elementary,
		      std::vector<shared_ptr<HedgeSurface> >& hedgesfs,
		      std::vector<HedgeSurface*>& prevsfs,
		      std::vector<std::vector<RevEngPoint*> >& out_groups);

    /// Approximate region points with a torus using context information if feasible, and 
    /// register the torus in the region. Possibly extract devinant points in connected groups
    /// \param  mainaxis Local coordinate system for model
    /// \param tol Approximation tolerance
    /// \param min_pt Not used
    /// \param min_pt_reg Minimum number of points left after extracting deviant points to create
    /// a surface
    /// \param angtol Angular tolerance
    /// \param prefer_elementary Used in checking if an already existing surface should be replaced
    /// \param hedgesfs New surface
    /// \param prevsfs Possibly previous surface being replaced
    /// \param out_groups Deviant points being removed from region
    bool contextTorus(Point mainaxis[3],
		      double tol, int min_pt, int min_pt_reg, 
		      double angtol, int prefer_elementary,
		      std::vector<shared_ptr<HedgeSurface> >& hedgesfs,
		      std::vector<HedgeSurface*>& prevsfs,
		      std::vector<std::vector<RevEngPoint*> >& out_groups);

    /// Approximate region points with a cylinder using context information if feasible, and 
    /// register the cylinder in the region. Possibly extract devinant points in connected groups
    /// \param  mainaxis Local coordinate system for model
    /// \param tol Approximation tolerance
    /// \param min_pt Not used
    /// \param min_pt_reg Minimum number of points left after extracting deviant points to create
    /// a surface
    /// \param angtol Angular tolerance
    /// \param prefer_elementary Used in checking if an already existing surface should be replaced
    /// \param hedgesfs New surface
    /// \param prevsfs Possibly previous surface being replaced
    /// \param out_groups Deviant points being removed from region
    bool contextCylinder(Point mainaxis[3],
			 double tol, int min_pt, int min_pt_reg,
			 double angtol, int prefer_elementary,
			 std::vector<shared_ptr<HedgeSurface> >& hedgesfs,
			 std::vector<HedgeSurface*>& prevsfs,
			 std::vector<std::vector<RevEngPoint*> >& out_groups);


    /// Approximate region points with a sphere or a torus using context information from an
    /// adjacent cylinder if feasible, and register the surface in the region. 
    /// Possibly extract devinant points in connected groups
    /// \param  mainaxis Local coordinate system for model
    /// \param tol Approximation tolerance
    /// \param min_pt Used to check if an approximation is accepted
    /// \param min_pt_reg Not used
    /// \param angtol Angular tolerance
    /// \param prefer_elementary Used in checking if an already existing surface should be replaced
    /// \param hedgesfs New surface
    /// \param prevsfs Possibly previous surface being replaced
    /// \param out_groups Deviant points being removed from region
    bool adjacentToCylinder(Point mainaxis[3],
			    double tol, int min_pt, int min_pt_reg,
			    double angtol, int prefer_elementary,
			    std::vector<shared_ptr<HedgeSurface> >& hedgesfs,
			    std::vector<HedgeSurface*>& prevsfs,
			    std::vector<std::vector<RevEngPoint*> >& out_groups);

    /// Return closest region point to a given position along with the distance between
    /// the points
    RevEngPoint* closestPoint(const Point& pos, double& dist);

    /// Return all region points that has neighbours belonging to a different region
    std::vector<RevEngPoint*> extractBdPoints();

    /// Return all region points that has neighbours belonging to a specified region
    std::vector<RevEngPoint*> extractBdPoints(std::vector<RevEngRegion*> regions);
    
    /// Return all region points that has neighbours belonging to more than one region
    /// different from current
    std::vector<RevEngPoint*> extractBranchPoints();

    /// Identify adjacent regions with and without an associated surface that are
    /// candidates for being included in the current region
    void getAdjCandMerge(std::vector<RevEngRegion*>& adj_surf,
			 std::vector<RevEngRegion*>& adj_nosurf);

    /// Estimate the length and width of a potential blend surface adjacent to this
    /// region along the curves in cvs
    void estimateBlendDimensions(std::vector<shared_ptr<CurveOnSurface> >& cvs,
				 std::vector<RevEngPoint*>& bd_points,
				 double tol, double distlim, 
				 std::vector<std::pair<double,double> >& t1_t2,
				 std::vector<double>& width, int& num_in_lim);

    /// Extract points lying within a specified distance to a given curve where
    /// the projection onto the curve falls within parameter interval <tmin, tmax>
    void getNearPoints(shared_ptr<CurveOnSurface>& cvs,
		       double& tmin, double& tmax, double width, double angtol,
		       std::vector<RevEngPoint*>& nearpoints);

    /// Extract points lying within a specified distance to a given set of curves
    void getNearPoints2(std::vector<RevEngPoint*>& points,
			shared_ptr<CurveOnSurface>& cv,
			double width, std::vector<RevEngPoint*>& nearpoints);

    /// Extract those points in points that satisfies given conditions with
    /// respect to the region surface
     std::vector<RevEngPoint*>
    removeOutOfSurf(std::vector<RevEngPoint*>& points,
		    double tol, double angtol, bool outer, double& min_dist);

    /// Grow current region, which has a surface that is realistically bounded by
    /// constant parameter curves (blend surface) with points from the region adjacent
    /// if they satisfy the accuracy constraints
    void growInDomain(RevEngRegion *adjacent, double tol, double angtol);

    /// Check if a current surface can possibly be replaced by a surface with
    /// better accuracy
    bool tryOtherSurf(int prefer_elementary, bool replace);
    
    /// Approximate region points with a spline surface if feasible, and register 
    /// the surface in the region. NB! Currently not in use
    /// \param tol Approximation tolerance
    /// \param min_pt Used in checking if the surface is feasible
    /// \param min_pt_reg Not used
    /// \param angtol Angular tolerance
    /// \param prefer_elementary Used in checking if an already existing surface should be replaced
    /// \param hedgesfs New surface
    /// \param prevsfs Possibly previous surface being replaced
    bool extractFreeform(double tol, int min_pt, int min_pt_reg,
			 double angtol, int prefer_elementary,
			 std::vector<shared_ptr<HedgeSurface> >& hedgesfs,
			 std::vector<HedgeSurface*>& prevsfs,
			 std::vector<std::vector<RevEngPoint*> >& out_groups);

    /// Divide composite region into more consistent pieces by extracting a planar piece
    /// and requiring the remaining parts to be connected
    void segmentByPlaneGrow(Point mainaxis[3], double tol,
			    double angtol, int min_pt, 
			    std::vector<shared_ptr<HedgeSurface> >& hedgesfs,
			    std::vector<HedgeSurface*>& prevsfs,
			    std::vector<std::vector<RevEngPoint*> >& out_groups);
    
    /// Divide composite region into more consistent pieces sorting points with respect from
    /// surface normals of adjacent planar surfaces. Regions with planar points are grown from
    /// this region. Cylindrical surfaces may be recognized.
    /// and requiring the remaining parts to be connected
    bool segmentByPlaneAxis(Point mainaxis[3], int min_point_in,
			    int min_pt_reg,
			    double tol, double angtol, int prefer_elementary,
			    std::vector<RevEngRegion*>& adj_planar,
			    std::vector<shared_ptr<HedgeSurface> >& hedgesfs,
			    std::vector<shared_ptr<RevEngRegion> >& added_reg,
			    std::vector<HedgeSurface*>& prevsfs,
			    std::vector<std::vector<RevEngPoint*> >& added_groups);

    /// Set pointer to associated surface. Remove possible previous surface
    void setHedge(HedgeSurface* surface)
    {
      associated_sf_.clear();
      associated_sf_.push_back(surface);
      computeDomain();
    }

    /// Add an additional associated surface. Obsolete!
    void addHedge(HedgeSurface* surface)
    {
      associated_sf_.push_back(surface);
      computeDomain();
    }

    /// Enquire if a region is associated a surface
    bool hasSurface()
    {
      return (associated_sf_.size() > 0);
    }

    /// Enquire the number of associated surfaces. Always one!
    int numSurface()
    {
      return (int)associated_sf_.size();
    }

    /// Return specified surface
    HedgeSurface* getSurface(int ix)
    {
      return (ix < 0 || ix>=(int)associated_sf_.size()) ? 0 : associated_sf_[ix];
    }

    /// Remove associated surface
    void clearSurface()
    {
      if (associated_sf_.size() > 0)
	associated_sf_.clear();
      maxdist_ = avdist_ = 0.0;
      num_inside_ = num_inside2_ = 0;
    }

    /// Compute accuracy and set associated surface
    void setAssociatedSurface(shared_ptr<ParamSurface>& surf,
			      double tol, double angtol, int min_pt_reg,
			      shared_ptr<HedgeSurface>& hedge);

    /// Modify region point collection. Used when the radius of a blend cylinder or
    /// torus is modified
    /// \param points_out Points to extract
    /// \param points_in Points to include
    void updateWithPointsInOut(std::vector<RevEngPoint*>& points_out,
			       std::vector<RevEngPoint*>& points_in,
			       double tol, double angtol);

    /// Move points from points into blend_points if appropriate depending on
    /// the distance to cvs and the flag in_blend
    void sortBlendPoints(std::vector<RevEngPoint*>& points,
			 std::vector<shared_ptr<CurveOnSurface> >& cvs,
			 double distance, bool in_blend,
			 std::vector<RevEngPoint*>& blend_points);
    
    /// Move points from points into this region or other region if appropriate
    /// Points with a distance to cvs larger than distance are moved
    void sortBlendPoints(std::vector<RevEngPoint*>& points,
			 std::vector<shared_ptr<CurveOnSurface> >& cvs,
			 double distance, RevEngRegion* other,
			 std::vector<RevEngPoint*>& blend_points1,
			 std::vector<RevEngPoint*>& blend_points2);

    /// Enquire flag on surface accuracy
    int getSurfaceFlag()
    {
      return surfflag_;
    }

    /// Set flag on surface accuracy
    void setSurfaceFlag(int surfflag)
    {
      surfflag_ = surfflag;
    }

    /// Enquire information about the history of the associated surface. Information not
    /// consistently carried out
    int getAdaptionHistory()
    {
      return surf_adaption_;
    }

    /// Return containing domain of the parameter values of the region points in the
    /// parameter domain of the associated surface.
    void getDomain(double dom[4])
    {
      for (int ka=0; ka<4; ++ka)
	dom[ka] = domain_[ka];
    }


    /// Set containing domain of parameter points
    void setDomain(double dom[4])
    {
      for (int ka=0; ka<4; ++ka)
	domain_[ka] = dom[ka];
    }

    /// Check if a fitted surface is accurate enough with respect to given tolerances,
    /// number of region points and computed accuracy. Old function. Use defineSfFlag
    bool accuracyOK(int min_pt, double tol, int num_inside, double avdist)
    {
      return (num_inside > min_pt && num_inside > (int)group_points_.size()/2 &&
	      avdist <= tol);
    }

    /// Check if a fitted surface is accurate enough with respect to given tolerances
    /// and computed accuracy. Old function. Use defineSfFlag
    bool accuracyOK(int num_points, int min_pt, double tol, int num_inside,
		    double avdist)
    {
      return (num_inside > min_pt && num_inside > num_points/2 &&
	      avdist <= tol);
    }

    /// Check if a fitted surface is accurate enough with respect to given tolerances,
    /// number of region points and computed accuracy.
    /// \param min_point Function returns not set if the number of OK points is less than min_point
    /// \param tol Approximation tolerance
    /// \param num_in Number of points that satisfy the distance and angular tolerances
    /// \param num_in2 Number of points that satisfy the distance tolerance
    /// \param avd Average absolute distance between the points and the surface
    /// \param type_cyl True if the surface is a cylinder or a cone
    int defineSfFlag(int min_point, double tol, int num_in,
		     int num_in2, double avd, bool type_cyl);
    
    /// Check if a fitted surface is accurate enough with respect to given tolerances,
    /// given number of points and computed accuracy.
    /// To be applied for a point group different from any region points
    /// \param num_points Number of points tested for accuracy
    /// \param min_point Function returns not set if the number of OK points is less than min_point
    /// \param tol Approximation tolerance
    /// \param num_in Number of points that satisfy the distance and angular tolerances
    /// \param num_in2 Number of points that satisfy the distance tolerance
    /// \param avd Average absolute distance between the points and the surface
    /// \param type_cyl True if the surface is a cylinder or a cone
    int defineSfFlag(int num_points, int min_point, double tol, int num_in,
		     int num_in2, double avd, bool type_cyl);

    /// Compute containing parameter domain of region points. Requires an associated surface
    void computeDomain();

    /// Bounding box for region points
    const BoundingBox& boundingBox()
    {
      return bbox_;
    }

    /// Enquire direction cone for LocFuncnormal_ (normal based on local approximating function)
    /// in the region points
    const DirectionCone& getNormalCone()
    {
      return normalcone_;
    }

    /// Enquire direction cone for triangulation normal in the region points
     const DirectionCone& getNormalConeTriang()
    {
      return normalcone2_;
    }

    /// Enquire range of principal curvatures in the region points
   void getPrincipalCurvatureInfo(double& mink1, double& maxk1, double& mink2, double& maxk2)
    {
      mink1 = mink1_;
      maxk1 = maxk1_;
      mink2 = mink2_;
      maxk2 = maxk2_;
    }

    /// Enquire average Gauss and mean curvatures in the region points
    /// \param avH Average mean curvature
    /// \param avK Average Gauss curvature
    /// \param MAH Average absolue mean curvature
    /// \param MAK Average absolute Gauss curvature
   void getAvCurvatureInfo(double& avH, double& avK, double& MAH, double& MAK)
    {
      avH = avH_;
      avK = avK_;
      MAH = MAH_;
      MAK = MAK_;
    }

    /// Enquire average LocFuncnormal_ (normal based on local approximating function)
    Point getMeanNormal()
    {
      return avnorm_;
    }

    /// Enquire average triangle normal
    Point getMeanNormalTriang()
    {
      return avnorm2_;
    }
    
    /// Set accuracy information of region points
    /// \param maxdist Maximum distance
    /// \param avdist Average absolute distance
    /// \param num_inside Number of points satisfying distance and angular tolerance,
    /// \param num_inside2 Mumber of points satisfying distance tolerance
    void setAccuracy(double maxdist, double avdist, int num_inside,
		     int num_inside2);

    /// Enquire accuracy information of region points
    /// \param maxdist Maximum distance
    /// \param avdist Average absolute distance
    /// \param num_inside Number of points satisfying distance and angular tolerance,
    /// \param num_inside2 Mumber of points satisfying distance tolerance
    void getAccuracy(double& maxdist, double& avdist, int& num_inside,
		     int& num_inside2)
    {
      maxdist = maxdist_;
      avdist = avdist_;
      num_inside = num_inside_;
      num_inside2 = num_inside2_;
    }

    /// Enquire average absolute distance betweeb region points and surface
    double getAverageDist()
    {
      return avdist_;
    }

    /// Enquire number of points satisfying distance and angular tolerance
    int getNumInside()
    {
      return num_inside_;
    }
    
    /// Enquire number of points satisfying distance tolerance
    int getNumInside2()
    {
      return num_inside2_;
    }

    /// Enquire maximum distance between region points and surface
    double getMaxSfDist()
    {
      return maxdist_;
    }

    /// Return distance to surface and angle between surface and point normal for all points
    void getDistAndAng(std::vector<std::pair<double,double> >& distang);

    /// In traversing points. Set as visited/not visited
    void setVisited(bool visited)
    {
      visited_ = visited;
    }

    /// In traversing points. Check if visited
    bool visited()
    {
      return visited_;
    }

    /// Check if the region is likely to correspond to a plane. Unstable information
    bool possiblePlane(double angtol, double inlim);
    
    /// Check if the region is likely to correspond to a cylinder. Unstable information
    bool possibleCylinder(double angtol, double inlim);
    
    /// Check if the region is likely to correspond to a cone. Unstable information
    bool possibleCone(double angtol, double inlim);
    
    /// Check if the region is likely to correspond to a torus. Unstable information
   bool possibleTorus(double angtol, double inlim);

    /// Check if the regions contains information for defining a swept surface
    bool hasSweepInfo()
    {
      return (sweep_.get() != 0);
    }

    /// Check if the regions contains information for defining a swept surface, and
    /// whether it is linear or rotational (not implemented) sweep
    int sweepType()
    {
      return (sweep_.get() ? sweep_->type_ : 0);
    }

    /// Enquire if the regions contains information on how to split this composite
    /// region into consistent pieces
    bool hasDivideInfo()
    {
      return (seg_info_.size() > 0);
    }

    /// Enquire number of information entities on how to split this composite
    /// region into consistent pieces
    int numDivideInfo()
    {
      return (int)seg_info_.size();
    }

    /// Define adjacency between this region and regions having points that are
    /// neighbours to this regions points
    void setRegionAdjacency();

    /// Update information about neighbouring regions
    void updateRegionAdjacency();

    /// Include this region into an adjacent region if feasible
    bool integrateInAdjacent(double mean_edge_len, int min_next,
			     int max_next, double tol, double angtol,
			     int max_nmb_outlier, RevEngRegion* taboo=0);

    /// Change parameterization of associated surface of type plane and update
    /// parameter values of region points
    void setPlaneParam(int min_pt_reg, Point mainaxis[3], double tol, double angtol);

    /// Given model coordinate system and axes from surfaces in adjacent region,
    /// modify the surface in this region to adapt to the most feasible of the candidate
    /// axes. Change surface if the accuracy of the new surface is sufficient
    bool updateSurfaceWithAxis(int min_pt_reg, Point adj_axis,
			       Point mainaxis[3], int ix, double tol,  
			       double angtol, Point pos);

    /// Add new adjacent region to the set of adjacent regions. NB! The adjacent region is not
    /// updated with information about this region
    void addAdjacentRegion(RevEngRegion* adj_reg)
    {
      //adj_reg->addAdjacentRegion(this);
      adjacent_regions_.insert(adj_reg);
    }

    /// Remove adjacent region from the set of adjacent regions. NB! The adjacent region is not
    /// updated 
    void removeAdjacentRegion(RevEngRegion* adj_reg)
    {
      //adj_reg->removeAdjacentRegion(this);
      if (std::find(adjacent_regions_.begin(), adjacent_regions_.end(), adj_reg) != adjacent_regions_.end())
	adjacent_regions_.erase(adj_reg);
      // else
      // 	std::cout <<"Something wrong in adjacent regions" << std::endl;
    }

    /// Remove all information about region adjacency. NB! adjacent regions are not updated
    void clearRegionAdjacency()
    {
      adjacent_regions_.clear();
    }

    /// Check if this region is adjacent to adj_reg
    bool isAdjacent(RevEngRegion* adj_reg)
    {
      return (adjacent_regions_.find(adj_reg) != adjacent_regions_.end());
    }

    /// Check if this region is adjacent to a regions that is adjacent to adj_reg
    bool isNextToAdjacent(RevEngRegion* adj_reg)
    {
      for (auto it=adjacent_regions_.begin(); it!=adjacent_regions_.end(); ++it)
	{
	  bool adjacent = (*it)->isAdjacent(adj_reg);
	  if (adjacent)
	    return true;
	}
      return false;
    }

    /// Fetch all regions adjacent to this region and adj_reg 
    std::vector<RevEngRegion*> commonAdjacent(RevEngRegion* adj_reg);

    /// Enquire number of adjacent regions
    int numAdjacentRegions()
    {
      return (int)adjacent_regions_.size();
    }

    /// Fetch all adjacent regions
    void getAdjacentRegions(std::vector<RevEngRegion*>& adjacent)
    {
      adjacent.insert(adjacent.end(), adjacent_regions_.begin(), adjacent_regions_.end());
    }

    /// Remove current regions from the adjacency information of all adjacent regions
     void removeFromAdjacent()
    {
      for (auto it=adjacent_regions_.begin(); it!=adjacent_regions_.end(); ++it)
	(*it)->removeAdjacentRegion(this);
    }

    /// Collect adjacent regions with and associated planar surface
    std::vector<RevEngRegion*> fetchAdjacentPlanar();

    /// Collect adjacent regions with and associated cylindrical surface
    std::vector<RevEngRegion*> fetchAdjacentCylindrical();

    /// Divide composite region into simpler pieces using information from
    /// adjacent regions
    bool segmentByAdjSfContext(Point mainaxis[3], int min_point_in, 
			       int min_pt_reg, double tol, double angtol,
			       std::vector<RevEngRegion*>& adj_planar,
			       std::vector<std::vector<RevEngPoint*> >& added_groups);

    /// Sort region points into groups with point normal being close to parallel or
    /// close to orthogonal to a set of axes (axis) and a group of remaining points
    bool sortByAxis(vector<Point>& axis, double tol, double axisang, double planeang,
		    std::vector<std::vector<RevEngPoint*> >& groups1,
		    std::vector<std::vector<RevEngPoint*> >& groups2,
		    std::vector<RevEngPoint*>& remaining);

    /// Identify cylindrical sub groups in a composite region and define cylindrical
    /// surfaces. The remaining points are split into connected groups
    bool extractCylByAxis(Point mainaxis[3], int min_point, int min_pt_reg,
			  double tol, double angtol, int prefer_elementary,
			  std::vector<shared_ptr<HedgeSurface> >& hedgesfs,
			  std::vector<shared_ptr<RevEngRegion> >& added_reg,
			  std::vector<std::vector<RevEngPoint*> >& out_groups,
			  std::vector<RevEngPoint*>& single_pts);

    /// Use segmentation information (SegmentData) to split composite region
    bool divideWithSegInfo(int seg_ix, int min_pt_reg,
			   std::vector<std::vector<RevEngPoint*> >& sep_groups,
			   std::vector<RevEngPoint*>& single_pts);

    /// Return surfaces axis/surface normal from the adjacent region with surface and
    /// most region points
    Point directionFromAdjacent(double angtol);

    /// Divide composite region using direction from adjacent and split region points
    /// according to the relation of their point normal to this direction
    bool segmentByDirectionContext(int min_point_in, double tol,
				   const Point& direction, double angtol,
				   std::vector<std::vector<RevEngPoint*> >& added_groups);

    /// Check if context information (surfaces in adjacent regions) indicate that this
    /// region could belong to a blend surface
    bool potentialBlend(double angtol);

    /// Fetch adjacent regions to this region that has a blend relation to cvs
    void neighbourBlends(std::vector<shared_ptr<CurveOnSurface> >& cvs,
			 double width, double tol,
			 std::vector<RevEngRegion*>& new_blends);

    /// Fetch adjacent regions with primary surfaces or primary base surfaces and surfaces
    void getAdjacentElemInfo(std::vector<std::pair<shared_ptr<ElementarySurface>, RevEngRegion*>  >& adj_elem,
			     std::vector<std::pair<shared_ptr<ElementarySurface>, RevEngRegion*>  >& adj_elem_base);

    /// Flag the region with information about the region from which it was extracted
    void setPreviousReg(RevEngRegion *prev)
    {
      prev_region_ = prev;
    }

    /// Extend the pool of trimming curves in this region with curves identifying
    /// boundaries between this region and adjacent regions with no surface or other
    /// regions where no boundaries between this region and the other region is defined
    void extendBoundaries(double mean_edge_len, int min_point_reg, double tol, 
			  double angtol, Point mainaxis[3]);
    

    /// Check if the region has a base surface defined. Partly displaced concept
    bool hasBaseSf()
    {
      return basesf_.get();
    }

    /// Fetch base surface. Used in the context of defining a model coordinate system
    shared_ptr<ParamSurface> getBase()
    {
      return basesf_;
    }
    
     /// Fetch base surface including accuracy information. 
    void getBase(shared_ptr<ParamSurface>& base, double& maxd, double& avd,
		 int& num_in, int& num_in2)
    {
      base = basesf_;
      maxd = maxdist_base_;
      avd = avdist_base_;
      num_in = num_in_base_;
      num_in2 = num_in_base2_;
    }
    
     /// Fetch accuracy information related to base surface. 
    void getBaseDist(double& maxd, double& avd, int& num_in, int& num_in2)
    {
      maxd = maxdist_base_;
      avd = avdist_base_;
      num_in = num_in_base_;
      num_in2 = num_in_base2_;
    }

    /// Update region with surface modified to adapt to a defined model axis
    /// and set related parameterization and accuracy information
    void updateSurfaceAndInfo(shared_ptr<ParamSurface> surf,
			      double tol, double angtol,
			      std::vector<double>& parvals,
			      std::vector<std::pair<double,double> >& dist_ang,
			      std::vector<RevEngEdge*>& nopar_edgs);

    /// Compute the approximation accuracy of this region with respect to a number of
    /// primary surfaces. Select the best and return accuracy information
    int checkSurfaceAccuracy(std::vector<shared_ptr<ElementarySurface> >& sfs, double tol,
			     double angtol, double& maxd, double& avd, int& num_in,
			     int& num2_in, int& sf_flag);

    /// Parameterize region points with respect to associated surface and update parameterization
    /// and accuracy information in points and region
    void parameterizePoints(double tol, double angtol);
    
    /// For each curve in cvs, fetch information on the part of the curve lying in the
    /// vicinity of the region points
    bool getCurveRestriction(std::vector<shared_ptr<CurveOnSurface> >& cvs,
			     double tol, double anglim,
			     std::vector<std::pair<double, double> >& endpars);

    /// Mark that this region should be associated to a future blend surface represented
    /// by blend_edge
    void setAssociatedBlend(RevEngEdge* blend_edge)
    {
      associated_blend_ = blend_edge;
    }

    /// Check if this region is associated with a blend
    bool hasAssociatedBlend()
    {
      return (associated_blend_ != 0);
    }

    /// Fetch the edge representing the blend with which this region is associated
    RevEngEdge* getAssociatedBlend()
    {
      return associated_blend_;
    }

    /// Remove the edge representing a blend from this region
    void removeAssociatedBlend()
    {
      associated_blend_ = 0;
    }

    /// Check if this region is adjacent to any (future) blend surfaces
    bool hasRevEdges()
    {
      return (rev_edges_.size() > 0);
    }

    /// Check the number of (future) blend surfaces associated to this region
    int numRevEdges()
    {
      return (int)rev_edges_.size();
    }

    /// Fetch a given edge representing a (future) blend surface
    RevEngEdge* getRevEdge(int ix)
    {
      if (ix < 0 || ix >= (int)rev_edges_.size())
	return 0;
      else
	return rev_edges_[ix];
    }

    /// Fetch all edges representing (future) blend surfaces
    std::vector<RevEngEdge*> getAllRevEdges()
    {
      return rev_edges_;
    }

    /// Extend the collection of edges representing (future) blend surfaces
    void addRevEdge(RevEngEdge* edge)
    {
      rev_edges_.push_back(edge);
    }

    /// Remove edge from the collection of edges representing (future) blend surfaces
    void removeRevEngEdge(RevEngEdge *edg)
    {
      auto it = std::find(rev_edges_.begin(), rev_edges_.end(), edg);
      if (it != rev_edges_.end())
	rev_edges_.erase(it);
    }

    /// Check if this region and the region other meets in a (future) blend surface
    bool commonRevEdge(RevEngRegion *other);
    
    /// Check if this region and the region other are separated by a trimming curve
    /// with twins
    bool commonTrimEdge(RevEngRegion *other);

    /// Mark that is region corresponds to a blend surface and set pointer to the edge
    /// representing this blend surface
    void setBlendEdge(RevEngEdge* edge)
    {
      blend_edge_ = edge;
    }

    /// Check if this region corresponds to a blend surface
    bool hasBlendEdge()
    {
      return (blend_edge_ != 0);
    }

    /// Fetch the edge containing information on e.g. the surfaces between which
    /// the surface associated to this region is a blend
    RevEngEdge* getBlendEdge()
    {
      return blend_edge_;
    }

    /// Remove blend information from this region
    void removeBlendEdge()
    {
      blend_edge_ = 0;
    }

    /// Increase the collection of trimming curves in this region with edge
    void addTrimEdge(shared_ptr<ftEdge> edge)
    {
      trim_edgs_.push_back(edge);
    }

    /// Check if this region has any trimming information
    bool hasTrimEdges()
    {
      return (trim_edgs_.size() > 0);
    }

    /// Fetch the number of trimming curve associated to this region
    int numTrimEdges()
    {
      return (int)trim_edgs_.size();
    }

    /// Fetch all trimming curves associated to this region
    std::vector<shared_ptr<ftEdge> > getTrimEdges()
    {
      return trim_edgs_;
    }

    /// Restrict the size of associated RevEngEdges with respect to the
    /// extent of the region points
    void adaptEdges();

    /// Mark this region to be removed from the region pool in Reveng
    void setRemove()
    {
      to_be_removed_ = true;
    }

    /// Check if this region is to be removed from the region pool in Reveng
    bool toBeRemoved()
    {
      return to_be_removed_;
    }

    /// Trim associated surface with respect to the trimming curves of this region,
    /// and replace the surface with a BoundedSurface if feasible
    bool trimSurface(double tol);


    /// Include the points of the region reg into this region, update parameterization
    /// according to parvals and accuracy information according to maxd, avd, num_inside
    /// num_inside2 and dist_ang, and update region adjacency information
    void includeAdjacentRegion(RevEngRegion* reg, double maxd, double avd,
			       int num_inside, int num_inside2,
			       std::vector<double>& parvals,
			       std::vector<std::pair<double, double> >& dist_ang,
			       std::vector<RevEngRegion*>& added_adjacent);

    /// Check if the current region surface should be replaced by a blend surface
    void checkEdgeAssociation(double tol, int min_point_reg,
			      std::vector<HedgeSurface*>& removed_sfs);

    /// Write debug information to file
    void writeRegionInfo(std::ostream& of);
    void writeRegionPoints(std::ostream& of);
    void writeAdjacentPoints(std::ostream& of);
    void writeUnitSphereInfo(std::ostream& of);
    void writeSubTriangulation(std::ostream& of);
    void writeSurface(std::ostream& of);

    /// Store current stage of region to file
    void store(std::ostream& os) const;

    /// Read region information from file
    void read(std::istream& is, shared_ptr<ftPointSet>& tri_sf,
	      std::vector<int>& associated_sf_id);
	      
  private:
    /// Unique id for region. Used in storing and reading data structure to and from file
    int Id_;

    /// Triangle vertices with additional information collected into this region
    /// Initially points belonging to classified segment. Updated during reverse engineering workflow
    std::vector<RevEngPoint*> group_points_;

    /// Selected method for classification in work flow
    int classification_type_;
    int edge_class_type_;

    /// Surface associated to this region. Only one surface is stored
    std::vector<HedgeSurface*> associated_sf_;

    /// Flag representing the accuracy with which the associated surface approximates the region
    /// points
    int surfflag_;

    /// History of surface approximation. Information not  consistently carried out
    int surf_adaption_;

    /// Rectangular domain surrounding the parameter pairs corresponding to the region points
    /// and the associated surface
    double domain_[4];

    /// Range of principal curvatures for region points
    double mink1_, maxk1_, mink2_, maxk2_;

    /// Average curvature in points: avH_ = mean, avK_ = Gauss, MAH_ = absolute valu of mean,
    /// MAK_ = absolute value of Gauss
    double avH_, avK_, MAH_, MAK_;

    /// Bounding box for region points
    BoundingBox bbox_;

    /// Direction cone for LocFunc normal in region points
    DirectionCone normalcone_;
    
    /// Direction cone for triangulation normal in region points
    DirectionCone normalcone2_;

    /// Fraction of LocFunc normal in points being closes to the average LocFunc normal
    /// than the given angular tolerance
    double frac_norm_in_;
    
    /// Fraction of triangulation normal in points being closes to the average triangulation normal
    /// than the given angular tolerance
    double frac_norm_in2_;

    /// Average LocFunc normal
    Point avnorm_;

    /// Average triangulation normal
    Point avnorm2_;

    /// Maximum and average absolute value of distance between region points and surface
    double maxdist_, avdist_;
    
    /// Number of points satisfying distance and angular tolerance and number of points
    /// satisfyng the distance tolerance
    int num_inside_, num_inside2_;

    /// Regions adjacent to this one
    std::set<RevEngRegion*> adjacent_regions_;

    /// Set if the region is extracted from a previous region to avoid it being integrated
    /// in the same region again
    RevEngRegion* prev_region_;

    /// Alternative approximating surface with accuracy information
    /// Can be set if the surface is modified due to adaption to model axes
    shared_ptr<ParamSurface> basesf_;
    double maxdist_base_, avdist_base_;
    int num_in_base_, num_in_base2_;

    /// Associated edges representing (future) blend surfaces between this region surface
    /// and another region surface
    std::vector<RevEngEdge*> rev_edges_;

    /// The region points is associated blend surface represented by this edge
    RevEngEdge* associated_blend_;

    /// The surface of this region is a blend surface associated to blend_edge
    RevEngEdge* blend_edge_;

    /// Collection of trimming edges intended for eventually bounding the region surface
    std::vector<shared_ptr<ftEdge> > trim_edgs_;

    /// Information that can be used to fit this regions points with a swept spline surface
    shared_ptr<SweepData> sweep_;

    /// Indicates that this region is visited in a search
    bool visited_;

    /// Indicates if this region is marked to be removed from the region pool in RevEng
    bool to_be_removed_;

    /// Information suitable for dividing this region that is found to be composite into
    /// simpler regions
    std::vector<shared_ptr<SegmentData> > seg_info_;;
    
    struct grow_cand
    {
      RevEngRegion *cand_;
      double maxd_, avd_, avang_;
      int num_in_, num2_in_, ang_in_;

      grow_cand(RevEngRegion* cand, double maxd, double avd, double avang,
		int num_in, int num2_in, int ang_in)
      {
	cand_ = cand;
	maxd_ = maxd;
	avd_ = avd;
	avang_ = avang;
	num_in_ = num_in;
	num2_in_ = num2_in;
	ang_in_ = ang_in;
      }
    };

    void integrateGrowCand(std::vector<grow_cand>& cand,
			   Point mainaxis[3], double tol,
			   double angtol, std::vector<RevEngRegion*>& grown_regions,
			   std::vector<HedgeSurface*>& adj_surfs);
    
    void analyseNormals(double tol, Point& normal, Point& centre, double& radius);
    void analysePlaneProperties(Point avnorm, double angtol,
				std::vector<RevEngPoint*>& in,
				std::vector<RevEngPoint*> out);
    void analyseCylinderProperties(Point avvec, double angtol,
				   std::vector<RevEngPoint*>& in,
				   std::vector<RevEngPoint*> out);
    void configSplit(std::vector<RevEngPoint*>& points,
		     std::vector<double>& param,
		     shared_ptr<Cylinder> cyl,
		     shared_ptr<SplineCurve> spl, double tol,
		     std::vector<std::vector<RevEngPoint*> >& configs);
    shared_ptr<Plane> computePlane(std::vector<RevEngPoint*>& points,
				   const Point& norm_dir, Point mainaxis[3]);
    shared_ptr<Cylinder>
    computeCylinder(std::vector<RevEngPoint*>& points, double tol);
    void analyseCylProject(shared_ptr<Cylinder> cyl, double tol,
			   std::vector<std::vector<RevEngPoint*> >& configs);
    void analyseCylRotate(shared_ptr<Cylinder> cyl, double tol,
			  double avdist, int num_in, double& avdist_lin,
			  int& num_in_lin, double& avdist_cub,
			  int& num_in_cub, shared_ptr<Cone>& cone);

    shared_ptr<Sphere> computeSphere(Point mainaxis[3], Point adj_axis,
				     std::vector<RevEngPoint*>& points);
    shared_ptr<Cone> computeCone(std::vector<RevEngPoint*>& points, Point& apex);
    shared_ptr<Torus> computeTorus(std::vector<RevEngPoint*>& points,
				   std::vector<Point>& adj_axis,
				   double tol, double angtol);
    shared_ptr<SplineSurface> computeLinearSwept(double tol, shared_ptr<SplineCurve>& profile,
						 Point& pt1, Point& pt2);
    shared_ptr<SplineSurface> computeFreeform(std::vector<RevEngPoint*>& points,
					      double tol);
    shared_ptr<SplineSurface> updateFreeform(std::vector<RevEngPoint*>& points,
					     double tol);
    void getPCA(double lambda[3], Point& eigen1, Point& eigen2, Point& eigen3);
    void getPCA(std::vector<RevEngPoint*>& points,
		double lambda[3], Point& eigen1, Point& eigen2, Point& eigen3);
    shared_ptr<SplineSurface> surfApprox(vector<RevEngPoint*>& points,
					 const BoundingBox& bbox);
    void splitCylinderRad(const Point& pos, const Point& axis,
			  const Point& Cx, const Point& Cy,
			  int nmb_split, std::vector<Point>& centr,
			  std::vector<double>& rad);
    void approximationAccuracy(std::vector<RevEngPoint*>& points,
			       shared_ptr<ParamSurface> surf,
			       double tol, double angtol,
			       double& maxd, double& avd,
			       std::vector<RevEngPoint*>& in,
			       std::vector<RevEngPoint*>& out);
    bool parameterizeOnSurf(std::vector<RevEngPoint*>& points,
			    shared_ptr<ParamSurface> surf,
			    std::vector<double>& data,
			    std::vector<double>& param,
			    int& inner1, int& inner2, bool& close1, bool& close2);
   bool parameterizeOnSurf(shared_ptr<ParamSurface> surf,
			    std::vector<double>& data,
			    std::vector<double>& param,
			    int& inner1, int& inner2, bool& close1, bool& close2);
    bool reparameterize(std::vector<RevEngPoint*>& points,
			std::vector<double>& param, std::vector<double>& param2,
			double& umin, double& umax, double& vmin, double& vmax);

    bool reparameterize(std::vector<double>& param, std::vector<double>& param2,
			double& umin, double& umax, double& vmin, double& vmax);

    void getParExtent(double curr[2], int pdir, std::vector<std::vector<int> >& raster,
		      int& i1, int& i2);

    void defineRaster(std::vector<double>& param, int nmb_div,
		      std::vector<std::vector<int> >& raster, double& umin,
		      double& umax, double& vmin, double& vmax);

    void extendInCorner(std::vector<double>& data, std::vector<double>& param,
			double umin, double umax, double vmin, double vmax);

    bool computeIntegrateInfo(std::vector<RevEngPoint*>& points, RevEngRegion *adj_reg,
			      double tol, double angtol, double radius, bool local_approx, 
			      int min_next, int max_next, int max_nmb_outlier, 
			      bool& outlier, int& nmb_pt_adj, double& maxdist, 
			      double& avdist, int& nmb_in, double& maxdist_adj, 
			      double& avdist_adj, int& nmb_in_adj,
			      int& nmb_in_adj2);

        bool
    analyseTorusContext(std::vector<std::pair<shared_ptr<ElementarySurface>, RevEngRegion*> >& adj,
			double tol, double angtol, std::vector<size_t>& adjacent_ix,
			int& plane_ix, int& cyl_ix, Point& pos, Point& axis,
			Point& Cx, double& R1, double& R2, double cyl_dom[4],
			bool& outer, bool& analyse_rotated);
    bool
    analyseCylinderContext(std::vector<std::pair<shared_ptr<ElementarySurface>, RevEngRegion*> >& adj,
			   double tol, double angtol, Point mainaxis[3], int mode,
			   std::vector<std::pair<shared_ptr<ElementarySurface>, RevEngRegion*> >& adj_planar,
			   Point& pos, Point& axis, Point& Cx, double& rad,
			   std::vector<RevEngPoint*>& cyl_pts);

    void computeFracNorm(double angtol, Point mainaxis[3], int nmb_axis[3],
			 double& in_frac1, double& in_frac2);
    
    shared_ptr<ElementarySurface>
    helicalComponent(shared_ptr<Cylinder> cyl,  double tol,
		     double angtol, int min_point,
		     int min_pt_reg, double avdist, int num_in2,
		     std::vector<std::pair<double,double> >& dist_ang,
		     std::vector<RevEngPoint*>& remaining,
		     std::vector<std::vector<RevEngPoint*> >& extracted,
		     double& maxd_out, double& avd_out,
		     int& num_in_out, int& num2_in_out,
		     std::vector<double>& parvals_out,
		     std::vector<std::pair<double,double> >& distang_out);

    bool defineHelicalInfo(shared_ptr<Cylinder> cyl,  double tol,
			   double angtol, int min_point, int min_pt_reg,
			   double avdist, int num_in1, int num_in2,
			   std::vector<std::pair<double,double> >& dist_ang,
			   std::vector<shared_ptr<HedgeSurface> >& hedgesfs,
			   std::vector<std::vector<RevEngPoint*> >& out_groups,
			   std::vector<RevEngPoint*>& single_pts);
    
    bool defineConeFromCyl(shared_ptr<Cylinder> cyl, double tol,
			   double angtol, int min_pt_reg,
			   double avdist, int num_in, int num2_in,
			   std::vector<shared_ptr<HedgeSurface> >& hedgesfs,
			   std::vector<std::vector<RevEngPoint*> >& out_groups,
			   std::vector<RevEngPoint*>& single_pts);
    
    bool integratePlanarPoints(std::vector<Point>& dir,
			       std::vector<std::vector<RevEngPoint*> >& groups,
			       std::vector<std::pair<shared_ptr<ElementarySurface>,RevEngRegion*> >& adj_elem,
			       double tol, double angtol,
			       std::vector<RevEngPoint*>& remaining);
    
    bool defineCylindricalRegs(Point mainaxis[3],
			       std::vector<std::vector<RevEngPoint*> >& groups,
			       int min_point, int min_pt_reg,
			       double tol, double angtol,
			       std::vector<shared_ptr<RevEngRegion> >& added_reg,
			       std::vector<shared_ptr<HedgeSurface> >& hedgesfs,
			       std::vector<RevEngPoint*>& remaining);
    
    void axisFromAdjacent(double angtol, std::vector<Point>& axis);
    
    void identifyOutPoints(std::vector<std::pair<double,double> >& distang,
			   double tol, double angtol, double angtol2,
			   std::vector<vector<RevEngPoint*> >& out_groups,
			   std::vector<RevEngPoint*>& remaining);

    void getRemainingPoints(std::vector<RevEngPoint*>& curr_pts,
			    std::vector<RevEngPoint*>& remaining);

    void adaptOneEdge(shared_ptr<ftEdge>& edge, double dom[4]);

    bool arrangeEdgeLoop(double tol, std::vector<int>& adjusted);
    
    shared_ptr<CurveOnSurface>
    constParSfCv(shared_ptr<ParamSurface> surf, int dir,
		 double par, int bd, double t1, double t2);

    void getDomainBoundaries(double tol, double angtol,
			     std::vector<std::pair<int, double> >& bd_par1,
			     std::vector<std::pair<int, double> >& bd_par2);
    
    void blendGrowFromAdjacent(RevEngRegion* adjacent,
			       std::vector<int>& pt_ix, double tol,
			       double angtol,
			       std::vector<RevEngPoint*>& grow_pt);
    
    // Select a point appropriate for starting the search for a sub group of points
    // that can be approximated by a plane
    RevEngPoint* seedPointPlane(int min_next, double rfac, double angtol);

    void updateRegion(double approx_tol, double anglim,
		      std::vector<RevEngRegion*>& adapted_regions,
		      std::vector<shared_ptr<RevEngRegion> >& outdiv_regions);

     void identifyAngPoints(std::vector<std::pair<double, double> >& dist_ang,
			    double tol, double disttol,
			    std::vector<RevEngPoint*>& ang_points);
    
     void identifyAngPoints(std::vector<std::pair<double, double> >& dist_ang,
			    double tol, double dtol, double dtol2,
			    std::vector<RevEngPoint*>& ang_points,
			    std::vector<RevEngPoint*>& remaining);
    
    void removeOtherPoints(std::vector<RevEngPoint*>& keep,
			   std::vector<HedgeSurface*>& prevsfs,
			   std::vector<std::vector<RevEngPoint*> >& out_groups);
    
    bool contextCylinder(Point mainaxis[3],
			 double tol, int min_pt, int min_pt_reg, 
			 double angtol, int prefer_elementary,
			 std::vector<std::pair<shared_ptr<ElementarySurface>, RevEngRegion*> >& adj_elem,
			 std::vector<shared_ptr<HedgeSurface> >& hedgesfs,
			 std::vector<HedgeSurface*>& prevsfs,
			 std::vector<std::vector<RevEngPoint*> >& out_groups,
			 int mode=1);

    RevEngPoint* closestParPoint(const Point& parpt, double& dist);

    std::vector<RevEngPoint*> extractNextToAdjacent(RevEngRegion* reg);

    void extractPointsAtSeam(std::vector<RevEngPoint*>& seam_pts1,
			     std::vector<RevEngPoint*>& seam_pts2, bool along_udir);
    
   void testBlendGap(std::vector<shared_ptr<CurveOnSurface> >& cvs,
		      double tmin, double tmax, double tdel, double width,
		      std::vector<std::pair<double,double> >& not_gap);

    void growFromNeighbour(Point mainaxis[3], int min_pt_reg,
			   std::vector<RevEngPoint*>& seed, double tol,
			   double angtol, RevEngRegion *neighbour,
			   bool do_update=true);

    shared_ptr<ParamSurface> surfaceWithAxis(std::vector<RevEngPoint*>& points,
					     Point axis, Point pos,
					     Point mainaxis[3]);

    bool adjustWithCylinder(Point mainaxis[3],
			    double tol, double angtol, int min_pt_reg,
			    std::vector<std::vector<RevEngPoint*> >& out_groups,
			    std::vector<RevEngRegion*>& grown_regions,
			    std::vector<HedgeSurface*>& adj_surfs);
    
    void getAdjInsideDist(shared_ptr<ParamSurface> surf, double dom[4],
			  double tol, RevEngRegion* reg,
			  double& avd, double& ava, int& nn,
			  std::vector<RevEngPoint*>& adjpts,
			  std::vector<double>& par_and_dist);

    bool identifySignificantAxis(std::vector<std::pair<shared_ptr<ElementarySurface>, RevEngRegion*> >& adj,
				 Point& pos, Point& axis, Point& axis2);

    void analyseRotated(Point& pos, Point& axis, Point& axis2);
    
    shared_ptr<Torus>
    analyseRotatedTorus(Point& pos, Point& Cx, Point& normal,
			double tol, double angtol);
    
    bool isInBlend(std::vector<shared_ptr<CurveOnSurface> >& cvs,
		   double width, double tol, int& in_blend);

    std::vector<RevEngPoint*> sortPtsSeq(double mean_edge_len,
					 std::vector<RevEngPoint*>& seq_pts,
					 std::vector<RevEngPoint*>& sub_pts);

    vector<RevEngPoint*>  extractBdOutPoints(shared_ptr<SplineCurve>& crv,
					     std::vector<RevEngPoint*>& seq_pts,
					     double tol);

    void adjustWithSurf(Point mainaxis[3], int min_pt_reg, double tol, double angtol);

    double getFracNormTriang()
    {
      return frac_norm_in2_;
    }
    
    void setBaseSf(shared_ptr<ParamSurface> base, double maxd, double avd,
		   int num_in, int num_in2)
    {
      basesf_ = base;
      maxdist_base_ = maxd;
      avdist_base_ = avd;
      num_in_base_ = num_in;
      num_in_base2_ = num_in2;
    }

    // Check if this region and the region adj both is connected to a triangulation
    // vertex (not in any region) that is defined as an edge point
    bool hasEdgeBetween(RevEngRegion* adj);

    void checkReplaceSurf(Point mainaxis[3], int min_pt_reg, double tol,
			  double angtol, bool always=false);

    void computeSurface(std::vector<RevEngPoint*>& points,
			Point mainaxis[3], double tol, double angtol,
			ClassType classtype, shared_ptr<ParamSurface>& updated,
			shared_ptr<ParamSurface>& updated2, bool& cyllike);
    
    void  curveApprox(std::vector<Point>& points, double tol,
		      shared_ptr<Circle> circle,
		      std::vector<double>& param,
		      shared_ptr<SplineCurve>& curve, Point& xpos);

    void getAdjacentBlends(std::vector<RevEngRegion*>& adj_blends);
    
    void rangeAlongAxis(const Point& pos, const Point& axis, double& tmin,
			double& tmax);
    
    shared_ptr<Plane>
    planarComponent(Point vec, int min_point, int min_pt_reg,
		    double tol, double angtol, Point mainaxis[3],
		    std::vector<RevEngPoint*>& remaining,
		    std::vector<std::vector<RevEngPoint*> >& extracted);
    
   bool planarComponent(Point vec, int min_point, int min_pt_reg, double tol,
			 double angtol, Point mainaxis[3],
			 std::vector<shared_ptr<HedgeSurface> >& hedgesfs,
			 std::vector<vector<RevEngPoint*> >& out_groups,
			 std::vector<RevEngPoint*>& single_pts,
			 bool create_surface = true);

};
}

#endif
