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

#ifndef _COONSPATCHVOLUMEGEN_H
#define _COONSPATCHVOLUMEGEN_H


#include <memory>
#include "GoTools/utils/Point.h"
#include "GoTools/geometry/SplineSurface.h"


namespace Go {


class SplineVolume;
class SplineSurface;
class SplineCurve;


/// This namespace contains functions used to create a Coons Patch volume

namespace CoonsPatchVolumeGen {


  /// Used internally in the geometry construction
  void get_corners(shared_ptr<Go::SplineSurface> surf, std::vector<Go::Point>& pts);
  /// Used internally in the geometry construction
  void push_corners(std::vector<Go::Point>& pts_to, const std::vector<Go::Point> pts_from,
		    int pos0, int pos1, int pos2, int pos3);
  /// Used internally in the geometry construction
  bool edge_curves_equal(shared_ptr<Go::SplineSurface> sf1, double par1, bool is_u_dir1,
			 shared_ptr<Go::SplineSurface> sf2, double par2, bool is_u_dir2,
			 double tol);


  /// Create a new SplineVolume representing the coons patch of six
  /// SplineSurfaces, the six faces of the volume.
  /// Here we assume that the faces are non-rational and lie in the
  /// same space, and all have the same B-spline basis
  /// (knot vectors and degrees are the same) along coinciding sides,
  /// i.e. for two surfaces with a common side, the two conciding
  /// edge curves have the same B-spline basis.
  /// Also, all B-spline bases are assumed to have knot interval
  /// [0.0, 1.0].
  /// \param surf_u_min The face/isosurface for u = 0.0
  /// \param surf_u_max The face/isosurface for u = 1.0
  /// \param surf_v_min The face/isosurface for v = 0.0
  /// \param surf_v_max The face/isosurface for v = 1.0
  /// \param surf_w_min The face/isosurface for w = 0.0
  /// \param surf_w_max The face/isosurface for u = 1.0
  /// \return a pointer to a newly created SplineVolume
  /// representing the coons patch.  The user assumes ownership of
  /// the object.
  SplineVolume* createCoonsPatchDirectly(const Go::SplineSurface* surf_u_min, const Go::SplineSurface* surf_u_max,
					 const Go::SplineSurface* surf_v_min, const Go::SplineSurface* surf_v_max,
					 const Go::SplineSurface* surf_w_min, const Go::SplineSurface* surf_w_max);



  /// Create a new SplineVolume representing the coons patch of six
  /// SplineSurfaces, the six faces of the volume.
  /// We test if the faces are non-rational, lie in the same space,
  /// and have conciding edge curves where needed. The B-spline bases
  /// of the surfaces might be reversed, swapped, have order raised,
  /// knots inserted and parameter interval rescaled to [0,1], in order to
  /// be able to use createCoonsPatchDirectly()
  /// \param surf_u_min One of faces/isosurfaces for the u-direction
  /// \param surf_u_max The other faces/isosurfaces for the u-direction
  /// \param surf_v_min One of faces/isosurfaces for the v-direction
  /// \param surf_v_max The other faces/isosurfaces for the v-direction
  /// \param surf_w_min One of faces/isosurfaces for the w-direction
  /// \param surf_w_max The other faces/isosurfaces for the w-direction
  /// \param tol Tolerance when identifiying equal curves and corners
  /// \return a pointer to a newly created SplineVolume
  /// representing the coons patch.  The user assumes ownership of
  /// the object
  SplineVolume* createCoonsPatch(const Go::SplineSurface* surf_u_min, const Go::SplineSurface* surf_u_max,
				 const Go::SplineSurface* surf_v_min, const Go::SplineSurface* surf_v_max,
				 const Go::SplineSurface* surf_w_min, const Go::SplineSurface* surf_w_max,
				 double tol = 1.0e-05);



} // namespace CoonsPatchVolumeGen


} // namespace Go


#endif // _COONSPATCHVOLUMEGEN_H
