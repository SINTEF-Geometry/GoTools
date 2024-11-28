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

#ifndef _REVENGPOINT_H
#define _REVENGPOINT_H

#include "GoTools/compositemodel/ftPointSet.h"
//#include "GoTools/compositemodel/RevEngRegion.h"
#include "GoTools/utils/Array.h"
#include "GoTools/utils/Point.h"
#include "GoTools/utils/DirectionCone.h"

namespace Go
{
  class RevEngRegion;
  
  /** RevEngPoint - A triangulation vertex enhanced with with information such as estimated 
      surface normal and curvature, classification flags and associated functionality. 
      A RevEngPoint knows which RevEngRegion it belongs to, if any. The reverse engineering
      functionality where RevEngPoint participates is run from RevEng.
   *
   */
  // Enumerations representing classification results

  /// Edge classification based on covariance. Not used
  enum
  {
   PCA_EDGE_UNDEF, PCA_NOT_EDGE, PCA_CLOSE_EDGE, PCA_EDGE
  };

  /// Edge classification based on curvature radius. Used
  enum
  {
   C1_EDGE_UNDEF, C1_NOT_EDGE, C1_CLOSE_EDGE, C1_EDGE
  };

  enum
  {
   C2_EDGE_UNDEF, C2_NOT_EDGE, C2_CLOSE_EDGE, C2_EDGE
  };

  /// Point classification based on covariance. Not used
  enum
  {
   PCA_UNDEF, PCA_PLANAR, PCA_LINEAR, PCA_OTHER
  };

  /// Point classification based on curvature. Used
  enum
  {
   C1_UNDEF, C1_PEAK, C1_RIDGE, C1_SRIDGE, C1_NONE, C1_FLAT, C1_MINSURF, C1_PIT, C1_VALLEY, C1_SVALLEY 
  };
    
  /** Subclass of ftSamplePoint; Enhanced point for use in reverse engineering
   */
  class RevEngPoint : public ftSamplePoint
  {

  public:
    /// Constructor
    RevEngPoint();

    /// Constuctor given point position
    /// \param bnd Boundary information. Inherited from ftSamplePoint. Not used
   RevEngPoint(Vector3D xyz, int bnd);

    /// Destructor
    virtual ~RevEngPoint();

    /// Compute normal based on the triangles meeting in the vertex
    /// \param lim Avoid very long triangle edges
    void computeTriangNormal(double lim);

    /// Set eigen vector and values in point
    void addCovarianceEigen(Point& eigen1, double lambda1, Point& eigen2,
			    double lambda2, Point& eigen3, double lambda3);

    /// Set normal and curvature information computed from expressing a group of points
    /// in the vicinity of the current in the coordinate system given by the eigen vector
    /// and approximate the points with a Bezier function. Surface normal minimum and
    /// maximum principal curvature vectors and values are computed as well as
    /// the distance from this point to the function and the average distance for all points
    /// in the group
    void addLocFuncInfo(Point& norm, Point& mincvec, double minc, Point& maxcvec,
		      double maxc, double currdist, double avdist);

    /// Enquire triangulation normal
    Point getTriangNormal()
    {
      Point zero(0.0,0.0,0.0);
      return normalcone_.greaterThanPi() ? zero :
	      normalcone_.centre();
    }

    /// Enquire the cone angle of a cone surrounding the surface normal of all
    /// triangles meeting in this vertex
    double getTriangAngle()
    {
      return normalcone_.greaterThanPi() ? 2.0*M_PI : normalcone_.angle();
    }

    /// Enquire the average length of triangle edges meeting in this vertex
    double getMeanEdgLen();

    /// Enquire the average length of triangle edges meeting in this vertex excluding
    /// edges with a length larger than maxlen
    double getMeanEdgLen(double maxlen);

    /// Enquire the normal computed from a local function approximation
    const Point& getLocFuncNormal()
    {
      return LocFuncnormal_;
    }

    /// Turn the normal computed from a local function approximation
    void turnLocFuncNorm()
      {
	LocFuncnormal_ *= -1;
      }


    /// Enquire first eigen vector
    const Point& getPCAEigen1()
    {
      return eigen1_;
    }

    /// Enquire normal given by covariance, e.g. the third eigen vector
     const Point& getPCANormal()
    {
      return eigen3_;
    }

    /// Turn direction of the coordinate system given by eigen vector
    void turnPCA()
    {
      eigen2_ *= -1;
      eigen3_ *= -1;
    }



    /// Traverse the triangulation and fetch points belonging to a given region, if given,
    /// within the sphere with radius radius from this vertex 
    /// \param min_nmb Minimum number of points required. The radius is
    /// increased if necessary to reach this number of points
    /// \param max_nmb When this number of points is reached, the search is stopped
   void fetchClosePoints2(double radius, int min_nmb, int max_nmb,
			   std::vector<RevEngPoint*>& nearpts,
			   RevEngRegion *region = 0);

    /// Fetch points belonging to the specified region connected through triangle edges,
    /// starting from current. Stop when max_nmb points are reached
    void fetchConnected(RevEngRegion *region, int max_nmb,
			std::vector<RevEngPoint*>& group);

     /// Fetch points marked with mark connected through triangle edges,
    /// starting from current. 
   void fetchConnectedMarked(int mark,
			      std::vector<RevEngPoint*>& group);

    /// Indicate that the points has been visited in a traversal
    void setVisited()
    {
      visited_ = 1;
    }

    /// A traversal is finished. Reset flag
    void unsetVisited()
    {
      visited_ = 0;
    }

    /// Check if a point has been visited in a traversal
    bool visited()
    {
      return (visited_ > 0);
    }

    /// Indicate that the point has been moved from one region to another. 
    void setMoved()
    {
      moved_ = 1;
    }

    /// Reset move flag
    void unsetMoved()
    {
      moved_ = 0;
    }

    /// Check if a ponit has been moved between regions
    bool moved()
    {
      return (moved_ > 0);
    }

    /// The eigen vectors and values are averaged between several computations.
    /// The number of computations are returned
    int nmbEigen()
    {
      return nmb_eigen_;
    }
    
    /// The results of the computation of normal and curvature by local funcations are
    /// averaged between several computations. The number of computations are returned
    int nmbLocFunc()
    {
      return nmb_locfunc_;
    }

    /// Enquire maximum principal curvature value
    double maxPrincipalCurvature()
    {
      return kmax_;
    }

    /// Enquire maximum principal curvature vector
    Point maxCurvatureVec()
    {
      return kvecmax_;
    }

    /// Enquire miniumum principal curvature value
    double minPrincipalCurvature()
    {
      return kmin_;
    }

    /// Enquire miniumum principal curvature vector
    Point minCurvatureVec()
    {
      return kvecmin_;
    }

    /// Enquire Gauss curvature
    double GaussCurvature()
    {
      return gausscurv_;
    }

    /// Enquire mean curvature
    double meanCurvature()
      {
	return meancurv_;
      }

    /// Gauss and mean curvature are smoothed to get less fragmented segmentation.
    /// Enquire previous value of Gauss curvature
    double GaussCurvature0()
    {
      return gausscurv0_;
    }

    /// Gauss and mean curvature are smoothed to get less fragmented segmentation.
    /// Enquire previous value of mean curvature
    double meanCurvature0()
      {
	return meancurv0_;
      }

    /// Set mean curvature. Used in smoothing of mean curvature
    void setMeanCurvature(double mean)
    {
      meancurv_ = mean;
    }

    /// Set Gauss curvature. Used in smoothing of Gauss curvature
    void setGaussCurvature(double Gauss)
    {
      gausscurv_ = Gauss;
    }

    /// Update curvature values. Used in iteration for smoothing Gauss and mean curvature
    void updateCurvature()
    {
      meancurv0_ = meancurv_;
      gausscurv0_ = gausscurv_;
    }

    double getCurvedness()
    {
      return curvedness_;
    }


    /// Set edge status in point using curvature radius
    void setEdgeClassification(int c1_edge)
    {
      edge_[0] = PCA_EDGE_UNDEF;
      edge_[1] = c1_edge;
      edge_[2] = C2_EDGE_UNDEF;
    }

    /// Set point classification using mean and Guass curvature
    void setClassification(int ctype)
    {
      surf_[0] = PCA_UNDEF;
      surf_[1] = ctype;
    }

    /// Check if a point is classified as edge. edge_class_type is expected to be C1_EDGE
    bool isEdge(int edge_class_type)
    {
      if (edge_class_type == 1)
	return (edge_[0] == PCA_EDGE);
      else if  (edge_class_type == 2)
	return (edge_[1] == C1_EDGE);
      else if  (edge_class_type == 3)
	return (edge_[2] == C2_EDGE);
      else
	return false;
    }

    /// Check if a point near being classified as edge. edge_class_type is expected to be C1_EDGE
    bool closeEdge(int edge_class_type)
    {
      if (edge_class_type == 1)
	return (edge_[0] >= PCA_CLOSE_EDGE);
      else if  (edge_class_type == 2)
	return (edge_[1] >= C1_CLOSE_EDGE);
      else if  (edge_class_type == 3)
	return (edge_[2] >= C2_CLOSE_EDGE);
      else
	return false;
    }

     /// Check if a point is different from an edge. edge_class_type is expected to be C1_EDGE
    bool notEdge(int edge_class_type)
    {
      if (edge_class_type == 1)
	return (edge_[0] <= PCA_NOT_EDGE);
      else if  (edge_class_type == 2)
	return (edge_[1] <= C1_NOT_EDGE);
      else if  (edge_class_type == 3)
	return (edge_[2] <= C2_NOT_EDGE);
      else
	return true;
    }

    /// Check if a point is classified as edge, but lacks other edge points
    /// in its vicinity
    bool isolatedEdge(int edge_class_type, int nmb, bool close);

    /// Unset edge information in point
    void setEdgeUndef()
    {
      edge_[0] = PCA_EDGE_UNDEF;
      edge_[1] = C1_EDGE_UNDEF;
      edge_[2] = C2_EDGE_UNDEF;
    }

    /// Get distance between point and the local function used to compute normal and
    /// curvature information
    double getPointDistance()
    {
      return ptdist_;
    }

    /// Get average distance between local point sets and the function used to compute normal 
    /// and curvature information
    double getAveragePointDistance()
    {
      return avdist_;
    }

    /// Enquire the point classification with respect to expected surface type
    int surfaceClassification(int classification_type) const;

    /// A point classified as edge is reclassified as near edge if the triangles meeting
    /// in the point don't indicate an edge
    void adjustWithTriangNorm(double anglim);

    /// Check if the point is assigned a region
    bool hasRegion()
    {
      return (region_ != 0);
    }

    /// Enquire the region assigned to the point
    RevEngRegion* region()
    {
      return region_;
    }

    /// Set region in point
    void setRegion(RevEngRegion* region)
    {
      region_ = region;
    }

    /// Unset region information
    void unsetRegion()
    {
      region_ = 0;
    }

    /// Unset accuracy and parameter information with respect to the surface
    /// approximating the points in the associated region
    void unsetSurfInfo()
    {
      sfdist_ = -1.0;
      sfang_ = -1.0;
      uv_ = Vector2D(0.0, 0.0);
    }

    /// Fetch regions different from the region of this point associated to the neighbouring
    /// points of this point
    void adjacentRegions(std::vector<RevEngRegion*>& adj) const;

    /// Check if any neighbour is associated the region reg
    bool nextToRegion(RevEngRegion *reg);

    /// Return classification of of point using Gauss and mean curvature 
    int C1_surf()
    {
      return surf_[1];
    }

    /// Mark point as outlier
    void setOutlier()
    {
      outlier_ = true;
    }

    /// Reset outlier information
    void unsetOutlier()
    {
      outlier_ = false;
    }

    /// Check if a point is classified as outlier
    bool isOutlier()
    {
      return outlier_;
    }

    /// The number of neighbours having the same classification as this point
    int nmbSameClassification(int classification_type) const;

    /// Record the distance between this point and the surface approximating the points
    /// of the associated region as well as the angular difference between surface normal
    /// and normal in point (least angle with respect to triangulation normal and local
    /// function normal)
    void setSurfaceDist(double dist, double ang)
    {
      sfdist_ = dist;
      sfang_ = ang;
    }

    /// Enquire distance to approximating surface
    double  getSurfaceDist()
    {
      return sfdist_;
    }
    
    /// Enquire distance and normal difference with respect to approximating surface
    void getSurfaceDist(double& dist, double& ang)
    {
      dist = sfdist_;
      ang = sfang_;
    }

    /// Check if any neighbouring vertex is associated the region reg
    bool isNeighbour(RevEngRegion* reg) const;

    /// Check if this vertex and pt are neighbours
    bool isNeighbour(RevEngPoint* pt) const;

    /// Fetch regions associated to neighbouring points that are associated a surface
    std::vector<RevEngRegion*> adjacentRegsWithSurf() const;

    /// Increase the number of moves between regions performed by this point
    void addMove()
    {
      nmb_move_++;
    }

    /// Enquire the number of moves between regions performed by this point
    /// (not complete information)
    int numMove()
    {
      return nmb_move_;
    }

    /// Fetch adjacent regions to this point and the closest point in those
    /// regions. Distant regions/points are excluded (distance depends on mean_edge_len)
    void getAdjInfo(double mean_edge_len, std::vector<RevEngRegion*>& adj_reg,
		    std::vector<RevEngPoint*>& adj_pt);

    /// Fetch regions associated to neighbouring points 
    std::vector<RevEngRegion*> getAdjacentRegions() const;

    
    /// Number of different regions associated to neighbouring points 
    int numAdjacentRegions() const;

    /// Include point in adjacent regions if the distance to this region is suffiently
    // small and the point normal information is consistent
    bool mergeWithAdjacent(double mean_edge_len);

    /// Set flag used in region growing and identification of connected groups
    void setMarkIx(int ix)
    {
      mark_ix_ = ix;
    }

    /// Unset flag used in region growing and identification of connected groups
    void unsetMarkIx()
    {
      mark_ix_ = -1;
    }

    /// Enquire flag used in region growing and identification of connected groups
    int getMarkIx()
    {
      return mark_ix_;
    }

    /// Store enhanced vertex information to file
    void store(std::ostream& os) const;

    /// Read enhanced vertex information from file. Return information about neighbouring vertices
    void read(std::istream& is, vector<int>& next_ix);
    
  private:
    /// Average length of edges meeting in this vertex
    double avedglen_;
    
    /// Eigenvectors of covariance matrix
    Point eigen1_, eigen2_, eigen3_;

    /// Eigenvalues of covariance matrix. lambda1_>= lambda2_ >= lambda3_
    double lambda1_, lambda2_, lambda3_;
    
    /// Normal and principal curvature vectors computed from local function approximation
    Point LocFuncnormal_, kvecmin_, kvecmax_;

    /// Principal curvatures computed from local function approximation
    double kmin_, kmax_;

    /// Distance to local function
    double ptdist_, avdist_;

    /// Number of point groups averaged to compute eigen vector/values, normal and
    /// curvature information
    int nmb_eigen_, nmb_locfunc_;

    /// Span of surface normals computed from triangulation
    DirectionCone normalcone_;  

    /// Previous and current version of mean curvature
    double meancurv0_, meancurv_;

    /// Previous and current version of Gauss curvature
    double gausscurv0_, gausscurv_;

    /// Curvedness computed from principal curvatures. Not used
    double curvedness_;   

    /// Parameters corresponding to classification results
    /// Results of edge classification. Sequence: Surface variation (PCA), curvature (C1),
    /// curveness (C2)
    int edge_[3];
    
    /// Results of surface classification. Sequence: PCA, curvature (C1), 
    int surf_[2];   // 
    
    /// Group (segment) of classified points
    RevEngRegion* region_;

    /// Outlier flag
    bool outlier_;

    /// Distance beteen point and surface approximating associated region
    double sfdist_;

    /// Angle between estimated normal in point and surface normal
    double sfang_;

    /// Number of times this point has been moved between region. Uncomplete information
    int nmb_move_;

    /// Flags mainly used in traversal
    mutable int visited_;
    mutable int moved_;
    mutable int mark_ix_;

    // Traverse the triangulation and fetch points within the sphere with
    // radius radius from this vertex 
    Point fetchClosePoints(double radius, int min_nmb, int max_nmb,
			   std::vector<Point>& nearpts);
    void getNearby(Vector3D xyz, double radius, int max_nmb,
		   std::vector<RevEngPoint*>& near,
		   RevEngRegion *region = 0);
  };
}

#endif
