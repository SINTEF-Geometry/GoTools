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

#ifndef _CLOSESTPOINTUTILS_H
#define _CLOSESTPOINTUTILS_H


#include <vector>
#include "GoTools/utils/Point.h"
#include "GoTools/geometry/GeomObject.h"


namespace Go
{

  /// Namespace for preprocessing structures on a surface model used to speed up closest point calculations.
  /// It consists of three classes:
  /// - BoundingBoxStructure, the overall class where one instance holds the enitre preprocess structure data
  /// - SurfaceData, one instance for each surface in the surface model
  /// - SubSurfaceBoundingBox, holding preprocessed data for one specific rectangular segment in a tensor mesh
  ///   sub structure of the parameter space on a surface (the underlying surface in case of a BoundedSurface)
  namespace boxStructuring
  {



    /// Class for preprocessed data on a specific surface in the model
    /// The surface parameter domain is split into a tensor mesh of sub segments. The splitting is neither
    /// required to be regular nor follow knot vector lines (when the underlying surface is a spline surface),
    /// but the current code creating the mesh structure chooses regular splittings for elementary
    /// surfaces and internal knot vector values for spline surfaces. For trimmed surfaces (BoundedSurface
    /// instances), the structure is created on the underlying surface, thus a segment might be either
    /// entirely inside the parameter domain, interely outside the parameter domain, or split by the boundary
    class SurfaceData
    {

    public:

      /// Constructor
      SurfaceData(shared_ptr<ParamSurface> surface)
      {
	surfaces_.push_back(surface);
      }

      /// Set the number of segments in the mesh structure in both parameter directions
      void setSegments(int segs_u, int segs_v)
      {
	segs_u_ = segs_u;
	segs_v_ = segs_v;
      }

      /// Get the number of segments in first parameter direction
      int segs_u() const
      {
	return segs_u_;
      }

      /// Get the number of segments in second parameter direction
      int segs_v() const
      {
	return segs_v_;
      }

      /// Set the index of this surface in the BoundingBoxStructure instance
      void setIndex(int index)
      {
	index_ = index;
      }

      /// Get the index of this surface in the BoundingBoxStructure instance
      int index() const
      {
	return index_;
      }

      /// Set number of surface copies
      void setSurfaceCopies(int nmb_copies)
      {
	if (nmb_copies < 1)
	  nmb_copies = 1;
	int old_nmb = surfaces_.size();
	surfaces_.resize(nmb_copies);
	for (int i = old_nmb; i < nmb_copies; ++i)
	  surfaces_[i] = shared_ptr<ParamSurface>(surfaces_[0]->clone());
      }

      /// Get a specific copy of the surface
      shared_ptr<ParamSurface> surface(int idx) const
	{
	  return surfaces_[idx];
	}

      /// Add internal surface point
      void add_inside_point(Point pt)
      {
	inside_points_.push_back(pt);
      }

      /// Get the internal surface points
      std::vector<Point> inside_points() const
	{
	  return inside_points_;
	}

    private:

      /// The index of this surface in the surface list in the overall BoundingBoxStructure instance
      int index_;

      /// The number of segments in first parameter direction
      int segs_u_;

      /// The number of segments in second parameter direction
      int segs_v_;

      /// The surface, might be cloned into copies to avoid evaluation errors when running multiple threads
      std::vector<shared_ptr<ParamSurface> > surfaces_;

      /// A set of points on the surface inside, but close to, the limiting curve loop, used to get an upper bound of the distance from a point to the surface
      /// Only used if the surface is a BoundedSurface
      std::vector<Point> inside_points_;

    };  // End class SurfaceData



    /// Class for the preprocessed data on a specific segment on the tensor mesh sub structure on the parameter
    /// space of a surface in the surface model
    class SubSurfaceBoundingBox
    {

    public:
      /// Constructor
      /// - surface_data is the surface data of the surface for this segment
      /// - pos_u and pos_v are the positions of this segment in the segment tensor mesh of the surface
      /// - box is the bounding box in the geometry space of the surface
      /// - par_domain holds the parameter domain of the segment (the bounding box in the parameter space)
      SubSurfaceBoundingBox(shared_ptr<SurfaceData> surface_data, int pos_u, int pos_v, BoundingBox box, shared_ptr<RectDomain> par_domain)
        : surface_data_(surface_data), domain_pos_u_(pos_u), domain_pos_v_(pos_v), box_(box), par_domain_(par_domain), domain_inside_boundary_(true)
      {
      }

      /// Set whether we know for sure this segment is entirely inside the parameter domain
      void setInside(bool inside)
      {
	domain_inside_boundary_ = inside;
      }

      /// Tell if the entire segment is known to be inside the parameter domain
      bool inside() const
      {
	return domain_inside_boundary_;
      }

      /// Tell if a specific point is guaranteed to be inside the parameter domain,
      /// by testing against an inside polygon of the entire parameter domain.
      /// Should always return false if the point is outside
      /// Should in most cases return true if the point is inside
      bool inside(double par_u, double par_v) const
      {
	if (domain_inside_boundary_)
	  return true;
	if (polygon_u_.size() == 0)
	  return false;

	// We go through the corner points in the polygon part inside the parameter domain
	// We do so by counting the quadrant changes of the polygon corners, viewed from the point
	// We count quadrants from 0 to 3, going clockwise and starting with 0 being the upper right
	// Quadrants 0 and 2 also include the axis aligned lines through the point
	int quadCount = 0;   // The total change of quadrants
	int n_pol = (int)polygon_u_.size();

	int previousQuad;
	for (int i = 0; i < n_pol; ++i)
	  {
	    double vec_u = polygon_u_[i] - par_u;
	    double vec_v = polygon_v_[i] - par_v;
	    int quad;
	    if (vec_u < 0)
	      quad = (vec_v > 0) ? 1 : 2;
	    else if (vec_u > 0)
	      quad = (vec_v < 0) ? 3 : 0;
	    else if (vec_v == 0.0)
	      return true;
	    else
	      quad = (vec_v < 0) ? 2 : 0;
	    if (i == 0)
	      {
		previousQuad = quad;
		continue;
	      }
	    if (quad == previousQuad)
	      continue;

	    int ch_quad = quad - previousQuad;
	    previousQuad = quad;
	    if (ch_quad == 3 || ch_quad == -1)
	      --quadCount;
	    else if (ch_quad == -3 || ch_quad == 1)
	      ++quadCount;
	    else
	      {
		// Crossing from one quadrant to the opposite, we must determin if we increase or decrease
		// the quadrant count by 2
		double cross_prod = vec_v * (polygon_u_[i-1] - par_u) - vec_u * (polygon_v_[i-1] - par_v);
		if (cross_prod == 0.0)
		  return true;  // Point is on the polygon
		quadCount += (cross_prod > 0) ? 2 : -2;
	      }
	  }

	if (quadCount != 0)
	  return quadCount > 0;

	// We end up with no quadrant change. To know if we are inside or outside, we must se if we move to the left or right
	// when going from the polygon start point to the end point, seen from the input point
	double end_pt_cross = (polygon_u_[0] - par_u) * (polygon_v_[n_pol - 1] - par_v) - (polygon_v_[0] - par_v) * (polygon_u_[n_pol - 1] - par_u);
	if (end_pt_cross != 0.0)
	  return end_pt_cross > 0.0;

	// We are left with the rare case when both the input point, polygon start point and end point lie on the same segment boundary line
	// We test by using the opposite to the input point
	double base_u = par_domain_->umax() + par_domain_->umin() - par_u;
	double base_v = par_domain_->vmax() + par_domain_->vmin() - par_v;
	return (polygon_u_[0] - base_u) * (polygon_v_[n_pol - 1] - base_v) > (polygon_v_[0] - base_v) * (polygon_u_[n_pol - 1] - base_u);
      }

      /// Get the segment position in the first direction in segment mesh of the surface
      int pos_u() const
      {
	return domain_pos_u_;
      }

      /// Get the segment position in the second direction in segment mesh of the surface
      int pos_v() const
      {
	return domain_pos_v_;
      }

      /// Get the geometry space bounding box of the image of the segment
      BoundingBox box() const
      {
	return box_;
      }

      /// Get the structure data of the surface
      shared_ptr<SurfaceData> surface_data() const
      {
	return surface_data_;
      }

      /// Get the domain of the segment in the surface parameter domain
      shared_ptr<RectDomain> par_domain() const
      {
	return par_domain_;
      }

      /// Add a point from the generated polygon inside the paramter domain of the surface
      /// The point is either inside or on the boudnary of the segment
      void add_polygon_corners(double par_u, double par_v)
      {
	polygon_u_.push_back(par_u);
	polygon_v_.push_back(par_v);
      }

      /// Tell if any of the generated polygon inside the parameter domain hits this segment
      bool has_polygon() const
      {
	return polygon_u_.size() > 0;
      }

      /// Get the number of polygon points
      int size_polygon() const
      {
	return (int)polygon_u_.size();
      }

      /// Remove the polygon information
      void remove_polygon()
      {
	polygon_u_.resize(0);
	polygon_v_.resize(0);
      }

    private:

      /// The structure data of the surface
      shared_ptr<SurfaceData> surface_data_;

      /// The segment position in the first direction in the segment mesh of the surface
      int domain_pos_u_;

      /// The segment position in the second direction in the segment mesh of the surface
      int domain_pos_v_;

      /// The geometry space bounding box of the image of the segment
      BoundingBox box_;

      /// The domain of the segment in the surface parameter domain
      shared_ptr<RectDomain> par_domain_;

      /// Variable telling if the entire segment is inside the parameter domain
      bool domain_inside_boundary_;

      /// The u-parameters of the points of the generated inside polygon of the parameter domain of the surface
      /// This segment only has the list of the points on the polygon part being inside the segment
      /// The first and last point will lie on the segment boundary, unless the entire polygon is a closed
      /// loop inside the segment. The polygon goes clockwise in the parameter domain.
      std::vector<double> polygon_u_;

      /// The v-paramters of the polygon corners
      std::vector<double> polygon_v_;

    };  // End class SubSurfaceBoundingBox



    /// Class for the preprocessed information of a surface model
    /// The gemetry space is split into voxels, axis aligned boxes
    /// of the same cubic shape (equal side length in all directions)
    /// The union of the voxels is rectangular, and holds all the surfaces
    class BoundingBoxStructure
    {

    public:

      /// Add information of a specific parameter sub segment of a surface to the structure
      void addBox(shared_ptr<SubSurfaceBoundingBox> box)
      {
	boxes_.push_back(box);
      }

      /// Add a surface from the surface model to the structure
      void addSurface(shared_ptr<SurfaceData> surface)
      {
	surface->setIndex((int)surfaces_.size());
	surfaces_.push_back(surface);
      }

      /// Get the numbder of segments in the structure
      int n_boxes() const
      {
	return (int)boxes_.size();
      }

      /// Get the numbder of surfaces in the structure
      int n_surfaces() const
      {
	return (int)surfaces_.size();
      }

      /// Get a specific segment in the structure
      shared_ptr<SubSurfaceBoundingBox> getBox(int i) const
      {
	return boxes_[i];
      }

      /// Get a specific surface in the structure
      shared_ptr<SurfaceData> getSurface(int i) const
      {
	return surfaces_[i];
      }

      /// Get the number of voxels in first coordinate direction
      int n_voxels_x() const
      {
	return n_voxels_x_;
      }

      /// Get the number of voxels in second coordinate direction
      int n_voxels_y() const
      {
	return n_voxels_y_;
      }

      /// Get the number of voxels in third coordinate direction
      int n_voxels_z() const
      {
	return n_voxels_z_;
      }

      /// Get the lower left corner of the entire voxel structure
      Point big_vox_low() const
      {
	return big_vox_low_;
      }

      /// Get the common length of all sides in the voxels
      double voxel_length() const
      {
	return voxel_length_;
      }

      /// Get all the segments that might hit a given voxel
      /// (those where the bounding boxes hit the voxel)
      std::vector<int> boxes_in_voxel(int i, int j, int k) const
      {
	return boxes_in_voxel_[i][j][k];
      }

      /// Set number of surface copies
      void setSurfaceCopies(int nmb_copies)
      {
	for (int i = 0; i < surfaces_.size(); ++i)
	  surfaces_[i]->setSurfaceCopies(nmb_copies);
      }

      /// Creates the voxel structure with information about the
      /// segment bounding boxes that hit each voxel
      /// bigbox is the entire bounding box of the surface structure,
      ///        the voxel structure must include bigbox
      /// volume_reduction gives the ratio of the volume of the bigbox structure.
      ///                  It will be a lower limit on the number of voxels
      void BuildVoxelStructure(BoundingBox bigbox, double volume_reduction)
      {
	// Get voxel dimensions
	Point diagonal = bigbox.high() - bigbox.low();
	double volume = diagonal[0] * diagonal[1] * diagonal[2];
	voxel_length_ = pow(volume/volume_reduction, 1.0/3.0);
	n_voxels_x_ = (int)(1.0 + diagonal[0] / voxel_length_);
	n_voxels_y_ = (int)(1.0 + diagonal[1] / voxel_length_);
	n_voxels_z_ = (int)(1.0 + diagonal[2] / voxel_length_);
	Point big_vox_center = bigbox.low() + diagonal * 0.5;
	big_vox_low_ = big_vox_center - Point((double)n_voxels_x_, (double)n_voxels_y_, (double)n_voxels_z_) * (0.5 * voxel_length_);

	// Add bounding boxes to voxel structure. Notice that a segment might
	// hit several voxels
	boxes_in_voxel_.resize(n_voxels_x_);
	for (int i = 0; i < n_voxels_x_; ++i)
	  {
	    boxes_in_voxel_[i].resize(n_voxels_y_);
	    for (int j = 0; j < n_voxels_y_; ++j)
	      boxes_in_voxel_[i][j].resize(n_voxels_z_);
	  }

	for (int i = 0; i < (int)boxes_.size(); ++i)
	  {
	    BoundingBox bb = boxes_[i]->box();
	    Point l_rel = (bb.low() - big_vox_low_) / voxel_length_;
	    Point h_rel = (bb.high() - big_vox_low_) / voxel_length_;
	    int l_x = (int)(l_rel[0]);
	    int l_y = (int)(l_rel[1]);
	    int l_z = (int)(l_rel[2]);
	    int h_x = (int)(h_rel[0]);
	    int h_y = (int)(h_rel[1]);
	    int h_z = (int)(h_rel[2]);
	    for (int jx = l_x; jx <= h_x; ++jx)
	      for (int jy = l_y; jy <= h_y; ++jy)
		for (int jz = l_z; jz <= h_z; ++jz)
		  boxes_in_voxel_[jx][jy][jz].push_back(i);
	  }
      }

      /// Test for closestPoint. Only used by the old code, closestVectorsOld()
      /// Will be removed if we know closestVectors() is safe
      bool closestPoint(int box_idx, bool any_tested, double best_dist, bool isInside, const Point& pt,
			double& clo_u, double& clo_v, Point& clo_pt, double& clo_dist,
			double epsilon, const RectDomain* domain_of_interest, double *seed) const
      {
	shared_ptr<SubSurfaceBoundingBox> surf_box = boxes_[box_idx];
	shared_ptr<SurfaceData> surf_data = surf_box->surface_data();
	shared_ptr<ParamSurface> paramSurf = surf_data->surface(0);
	shared_ptr<BoundedSurface> boundedSurf = dynamic_pointer_cast<BoundedSurface>(paramSurf);
	bool shall_test = true;
	if (isInside)
	  paramSurf = boundedSurf->underlyingSurface();
	else if (any_tested)
	  {
	    boundedSurf->underlyingSurface()->closestPoint(pt, clo_u, clo_v, clo_pt, clo_dist, epsilon, domain_of_interest, seed);
	    shall_test = clo_dist < best_dist;
	  }
	if (!shall_test)
	  return false;

	paramSurf->closestPoint(pt, clo_u, clo_v, clo_pt, clo_dist, epsilon, domain_of_interest, seed);
	return true;
      }

    private:

      /// The common length of all sides in the voxels
      double voxel_length_;

      /// The number of voxels in first coordinate direction
      int n_voxels_x_;

      /// The number of voxels in second coordinate direction
      int n_voxels_y_;

      /// The number of voxels in third coordinate direction
      int n_voxels_z_;

      /// The lower left corner of the entire voxel structure
      Point big_vox_low_;

      /// The parameter sub segments of the structure
      std::vector<shared_ptr<SubSurfaceBoundingBox> > boxes_;

      /// The surfaces of the structure
      std::vector<shared_ptr<SurfaceData> > surfaces_;

      /// The segments that hit each voxel
      std::vector<std::vector<std::vector<std::vector<int> > > > boxes_in_voxel_;

    };  // End class BoundingBoxStructure



  } // namespace Go::boxStructuring


  /// Create preprocessing data for the closest vector calculations on a surface.
  /// surfaces - a collection of the paramteric surfaces defining the surface model. Only instances of the ParamSurface subclass hierarchy are used
  /// par_len_el - a guiding for the side lengths of the segments in geometry space, used to determine the number of segments for elementary surfaces
  /// returns the preprocessing structures used as input for the closest point calculations
  shared_ptr<boxStructuring::BoundingBoxStructure> preProcessClosestVectors(const std::vector<shared_ptr<GeomObject> >& surfaces, double par_len_el);


  void closestPointSingleCalculation(int pt_idx, int start_idx, int skip,
				     const std::vector<float>& inPoints,
				     const std::vector<std::vector<double> >& rotationMatrix, const Point& translation,
				     const shared_ptr<boxStructuring::BoundingBoxStructure>& boxStructure,
				     std::vector<float>& result, std::vector<std::vector<int> >& lastBoxCall,
				     int return_type, int search_extend);

  /// Calculates the closest points of a point cloud to a surface model, after a SO(3)-rotation and translation is applied on the point clod.
  /// The method uses polygons inside the bounding curves on paramter domains to help determining if parameter pairs are inside the
  /// parameter domain. NB! The creation of the polygons is not yet proven to guarantee inside polygons, thus there is a theoretical
  /// risk of not getting the closest point in every case.
  /// pts            - The point cloud, of length 3N where N is the number of points, on format p[0][0], p[0][1], p[0][2], p[1][0] , ...
  /// structure      - the preprocessed structure used to improve the calculation speed. This also holds the surface model.
  /// rotationMatrix - An orthogonal 3x3 matrix describing the rotation to be applied in the point cloud before starting the calculations
  /// translation    - A translation vector to be added to the point cloud (after the orthogonal rotation) before starting the calculations
  /// return_type    - Tell whether the distances (0), signed distances (1) or closest points (2) should be returned
  /// start_idx      - Used for defining the subset of the points on which the calculation should be performed, see below
  /// skip           - Used for defining the subset of the points on which the calculation should be performed, see below
  /// max_idx        - Used for defining the subset of the points on which the calculation should be performed. The subset consists of the
  ///                  points with index    idx = skip + N * step   such that start_idx <= idx < max_idx. For calculations on the entire
  ///                  point cloud, use start_idx = 0, skip = 1 and max_idx >= number of points
  /// search_extend  - Used to define the number of segments to be added in each direction when defining the parameter subset on which
  ///                  the closest point functions should be performed. Will be removed.
  /// m_core         - Whether the calculations should be performed in parallell on multiple cores (only if OPENMP is included)
  /// returns   A vector of the distances to the closest points for the subset on which the calculations are performed.
  std::vector<float> closestPointCalculations(const std::vector<float>& pts, const shared_ptr<boxStructuring::BoundingBoxStructure>& structure,
					      const std::vector<std::vector<double> >& rotationMatrix, const Point& translation,
					      int return_type, int start_idx, int skip, int max_idx, int search_extend = 3, bool m_core = true);


  /// Calculates the closest points of a point cloud to a surface model, by not using the inside polygons in closestVectors()
  std::vector<float> closestVectorsOld(const std::vector<float>& inPoints, const shared_ptr<boxStructuring::BoundingBoxStructure>& boxStructure,
				       const std::vector<std::vector<double> >& rotationMatrix, const Point& translation,
				       int return_type, int start_idx, int skip, int max_idx, int search_extend = 3);




  /// Make closest point calculations on the entire set of a point cloud.
  /// Fot return_type == 0, the distances are returned
  /// Fot return_type == 1, the signed distances are returned
  /// Fot return_type == 2, the points are returned
  std::vector<float> closestPointCalculations(const std::vector<float>& pts, const shared_ptr<boxStructuring::BoundingBoxStructure>& structure,
					      const std::vector<std::vector<double> >& rotationMatrix, const Point& translation,
					      int return_type);

  /// Calculates the distance for each point in a point cloud to a given surface model, after a SO(3)-rotation and translation
  /// is applied on the point cloud.
  /// pts            - The point cloud, of length 3N where N is the number of points, on format p[0][0], p[0][1], p[0][2], p[1][0] , ...
  /// structure      - the preprocessed structure used to improve the calculation speed. This also holds the surface model.
  /// rotationMatrix - An orthogonal 3x3 matrix describing the rotation to be applied in the point cloud before starting the calculations
  /// translation    - A translation vector to be added to the point cloud (after the orthogonal rotation) before starting the calculations
  /// returns a vector of length pts.size()/3, holding the distances in the same order as the input points
  std::vector<float> closestDistances(const std::vector<float>& pts, const shared_ptr<boxStructuring::BoundingBoxStructure>& structure,
				      const std::vector<std::vector<double> >& rotationMatrix, const Point& translation);

  /// Calculates the signed distance for each point in a point cloud to a given surface model, after a SO(3)-rotation and translation
  /// is applied on the point cloud.
  /// For every point the signed distance
  /// * has absolute value equal to the distance between the point and the closest point on the surface model
  /// * is positive if the point is outside the model
  /// * is negative if the point is inside the model
  /// pts            - The point cloud, of length 3N where N is the number of points, on format p[0][0], p[0][1], p[0][2], p[1][0] , ...
  /// structure      - the preprocessed structure used to improve the calculation speed. This also holds the surface model.
  /// rotationMatrix - An orthogonal 3x3 matrix describing the rotation to be applied in the point cloud before starting the calculations
  /// translation    - A translation vector to be added to the point cloud (after the orthogonal rotation) before starting the calculations
  /// returns a vector of length pts.size()/3, holding the signed distances in the same order as the input points
  std::vector<float> closestSignedDistances(const std::vector<float>& pts, const shared_ptr<boxStructuring::BoundingBoxStructure>& structure,
					    const std::vector<std::vector<double> >& rotationMatrix, const Point& translation);

  /// Calculates the closest point for each point in a point cloud to a given surface model, after a SO(3)-rotation and translation
  /// is applied on the point cloud.
  /// pts            - The point cloud, of length 3N where N is the number of points, on format p[0][0], p[0][1], p[0][2], p[1][0] , ...
  /// structure      - the preprocessed structure used to improve the calculation speed. This also holds the surface model.
  /// rotationMatrix - An orthogonal 3x3 matrix describing the rotation to be applied in the point cloud before starting the calculations
  /// translation    - A translation vector to be added to the point cloud (after the orthogonal rotation) before starting the calculations
  /// returns a vector of same length as pts, holding the closest point coordinates in the same order as the input points
  std::vector<float> closestPoints(const std::vector<float>& pts, const shared_ptr<boxStructuring::BoundingBoxStructure>& structure,
				   const std::vector<std::vector<double> >& rotationMatrix, const Point& translation);

  /// Calculates the signed distance for each point in a point cloud to a given surface model, after a SO(3)-rotation and translation
  /// is applied on the point cloud.
  /// For every point the signed distance
  /// * has absolute value equal to the distance between the point and the closest point on the surface model
  /// * is positive if the point is outside the model
  /// * is negative if the point is inside the model
  /// pts            - The point cloud, of length 3N where N is the number of points, on format p[0][0], p[0][1], p[0][2], p[1][0] , ...
  /// structure      - the preprocessed structure used to improve the calculation speed. This also holds the surface model.
  /// rotationMatrix - An orthogonal 3x3 matrix describing the rotation to be applied in the point cloud before starting the calculations
  /// translation    - A translation vector to be added to the point cloud (after the orthogonal rotation) before starting the calculations
  /// returns a vector of length 4*pts.size()/3, holding the signed distances in the same order as the input points, as well as
  /// index of closest surface and the closest u & v parameters, stored as 4-tuples.
  std::vector<float> closestSignedDistanceSfParams(const std::vector<float>& inPoints,
                                                   const shared_ptr<boxStructuring::BoundingBoxStructure>& boxStructure,
                                                   const std::vector<std::vector<double> >& rotationMatrix, const Point& translation);

} // namespace Go


#endif // _CLOSESTPOINTUTILS_H
