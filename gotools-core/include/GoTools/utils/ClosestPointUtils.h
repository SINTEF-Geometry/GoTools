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

  namespace boxStructuring
  {

  class SurfaceData
  {

  public:

    SurfaceData(shared_ptr<ParamSurface> surface)
      : surface_(surface)
    {
    }

    void setSegments(int segs_u, int segs_v)
    {
      segs_u_ = segs_u;
      segs_v_ = segs_v;
    }

    int segs_u() const
    {
      return segs_u_;
    }

    int segs_v() const
    {
      return segs_v_;
    }

    void setIndex(int index)
    {
      index_ = index;
    }

    int index() const
    {
      return index_;
    }

    shared_ptr<ParamSurface> surface() const
    {
      return surface_;
    }

    void add_inside_point(Point pt)
    {
      inside_points_.push_back(pt);
    }

    std::vector<Point> inside_points() const
    {
      return inside_points_;
    }

  private:

    int index_;

    int segs_u_;

    int segs_v_;

    shared_ptr<ParamSurface> surface_;

    std::vector<Point> inside_points_;

  };  // End class SurfaceData

  class SubSurfaceBoundingBox
  {

  public:
    SubSurfaceBoundingBox(shared_ptr<SurfaceData> surface_data, int pos_u, int pos_v, BoundingBox box, shared_ptr<RectDomain> par_domain)
      : surface_data_(surface_data), domain_pos_u_(pos_u), domain_pos_v_(pos_v), box_(box), par_domain_(par_domain), domain_inside_boundary_(true)
    {
    }

    void setInside(bool inside)
    {
      domain_inside_boundary_ = inside;
    }

    bool inside() const
    {
      return domain_inside_boundary_;
    }

    BoundingBox box() const
    {
      return box_;
    }

    shared_ptr<SurfaceData> surface_data() const
    {
      return surface_data_;
    }

    int pos_u() const
    {
      return domain_pos_u_;
    }

    int pos_v() const
    {
      return domain_pos_v_;
    }

    shared_ptr<RectDomain> par_domain() const
    {
      return par_domain_;
    }

  private:

    BoundingBox box_;

    shared_ptr<RectDomain> par_domain_;

    shared_ptr<SurfaceData> surface_data_;

    int domain_pos_u_;

    int domain_pos_v_;

    bool domain_inside_boundary_;

  };  // End class SubSurfaceBoundingBox

  class BoundingBoxStructure
  {

  public:
    void addBox(shared_ptr<SubSurfaceBoundingBox> box)
    {
      boxes_.push_back(box);
    }

    void addSurface(shared_ptr<SurfaceData> surface)
    {
      surface->setIndex(surfaces_.size());
      surfaces_.push_back(surface);
    }

    int n_boxes() const
    {
      return boxes_.size();
    }

    int n_surfaces() const
    {
      return surfaces_.size();
    }

    shared_ptr<SubSurfaceBoundingBox> getBox(int i) const
    {
      return boxes_[i];
    }

    shared_ptr<SurfaceData> getSurface(int i) const
    {
      return surfaces_[i];
    }

    int n_voxels_x() const
    {
      return n_voxels_x_;
    }

    int n_voxels_y() const
    {
      return n_voxels_y_;
    }

    int n_voxels_z() const
    {
      return n_voxels_z_;
    }

    Point big_vox_low() const
    {
      return big_vox_low_;
    }

    double voxel_length() const
    {
      return voxel_length_;
    }

    std::vector<int> boxes_in_voxel(int i, int j, int k) const
    {
      return boxes_in_voxel_[i][j][k];
    }
    
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
	
      // Add bounding boxes to voxel structure
      boxes_in_voxel_.resize(n_voxels_x_);
      for (int i = 0; i < n_voxels_x_; ++i)
	{
	  boxes_in_voxel_[i].resize(n_voxels_y_);
	  for (int j = 0; j < n_voxels_y_; ++j)
	    boxes_in_voxel_[i][j].resize(n_voxels_z_);
	}

      for (int i = 0; i < boxes_.size(); ++i)
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

    // Do closestPoint test
    bool closestPoint(int box_idx, bool any_tested, double best_dist, bool isInside, const Point& pt,
		      double& clo_u, double& clo_v, Point& clo_pt, double& clo_dist,
		      double epsilon, const RectDomain* domain_of_interest, double *seed) const
    {
      shared_ptr<SubSurfaceBoundingBox> surf_box = boxes_[box_idx];
      shared_ptr<SurfaceData> surf_data = surf_box->surface_data();
      shared_ptr<ParamSurface> paramSurf = surf_data->surface();
      shared_ptr<BoundedSurface> boundedSurf = dynamic_pointer_cast<BoundedSurface>(paramSurf);
      bool shall_test = true;
      if (isInside)
	paramSurf = boundedSurf->underlyingSurface();
      else if (any_tested)
	{
	  boundedSurf->underlyingSurface()->closestPoint(pt, clo_u, clo_v, clo_pt, clo_dist, epsilon, domain_of_interest, seed);
	  /*
	  std::cout << " UL : surf = " << (surf_data->index()) << " seed = (" << seed[0] << ", " << seed[1]
		    << ") clo_par = (" << clo_u << ", " << clo_v << ")  clo_dist = " << clo_dist
		    << "  domain = [" << (domain_of_interest->umin()) << ", " << (domain_of_interest->umax()) << "]x["
		    << (domain_of_interest->vmin()) << ", " << (domain_of_interest->vmax()) << "]" << std::endl;
	  */
	  shall_test = clo_dist < best_dist;
	}
      if (!shall_test)
	return false;

      paramSurf->closestPoint(pt, clo_u, clo_v, clo_pt, clo_dist, epsilon, domain_of_interest, seed);
      /*
      std::cout << " TP : surf = " << (surf_data->index()) << " seed = (" << seed[0] << ", " << seed[1]
		<< ") clo_par = (" << clo_u << ", " << clo_v << ")  clo_dist = " << clo_dist
		<< "  domain = [" << (domain_of_interest->umin()) << ", " << (domain_of_interest->umax()) << "]x["
		<< (domain_of_interest->vmin()) << ", " << (domain_of_interest->vmax()) << "]" << std::endl;
      */
      return true;
    }

  private:

    double voxel_length_;

    int n_voxels_x_;

    int n_voxels_y_;

    int n_voxels_z_;

    Point big_vox_low_;

    std::vector<shared_ptr<SubSurfaceBoundingBox> > boxes_;

    std::vector<shared_ptr<SurfaceData> > surfaces_;

    std::vector<std::vector<std::vector<std::vector<int> > > > boxes_in_voxel_;

  };  // End class BoundingBoxStructure



  } // namespace Go::boxStructuring

  shared_ptr<boxStructuring::BoundingBoxStructure> preProcessClosestVectors(const std::vector<std::shared_ptr<GeomObject> >& surfaces, double par_len_el);


  std::vector<float> closestVectors(const std::vector<float>& inPoints, const shared_ptr<boxStructuring::BoundingBoxStructure>& boxStructure,
				    const std::vector<std::vector<double> >& rotationMatrix, const Point& translation,
				    int test_type, int start_idx, int skip, int max_idx, int search_extend = 3);


  std::vector<float> closestVectorsOld(const std::vector<float>& inPoints, const shared_ptr<boxStructuring::BoundingBoxStructure>& boxStructure,
				       const std::vector<std::vector<double> >& rotationMatrix, const Point& translation,
				       int test_type, int start_idx, int skip, int max_idx, int search_extend = 3);



} // namespace Go


#endif // _CLOSESTPOINTUTILS_H
