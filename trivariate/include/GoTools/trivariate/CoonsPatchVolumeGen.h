//===========================================================================
//
// File : CoonsPatchVolumeGen.h
//
// Created: Fri Jun 19 09:03:16 2009
//
// Author: Kjell Fredrik Pettersen
//
// Revision: $Id: CoonsPatchVolumeGen.h,v 1.1 2009/06/24 10:42:07 kfp Exp $
//
// Description:
//
//===========================================================================



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


  void get_corners(std::shared_ptr<Go::SplineSurface> surf, std::vector<Go::Point>& pts);
  void push_corners(std::vector<Go::Point>& pts_to, const std::vector<Go::Point> pts_from,
		    int pos0, int pos1, int pos2, int pos3);
  bool edge_curves_equal(std::shared_ptr<Go::SplineSurface> sf1, double par1, bool is_u_dir1,
			 std::shared_ptr<Go::SplineSurface> sf2, double par2, bool is_u_dir2,
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
