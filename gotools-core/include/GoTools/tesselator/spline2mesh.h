#ifndef SPLINE2MESH_H_INCLUDED
#define SPLINE2MESH_H_INCLUDED






#include "GoTools/utils/Array.h"

#include <vector>


#include "GoTools/tesselator/2dpoly_for_s2m.h"
#include "GoTools/geometry/ParamSurface.h"
#include "GoTools/geometry/ParamCurve.h"






namespace Go
{
  
  // 081206: A version for more than one curve.
  void make_trimmed_mesh(shared_ptr<ParamSurface> srf, 
			 std::vector<shared_ptr<ParamCurve> >& crv_set,
			 std::vector< Vector3D > &vert,
			 std::vector< Vector2D > &vert_p,
			 std::vector< int > &bd,
			 std::vector< Vector3D > &norm,
			 std::vector<int> &mesh,
			 std::vector< Vector3D > &trim_curve,
			 std::vector< Vector3D > &trim_curve_p,
			 const int dn, const int dm,
			 double bd_res_ratio);

} // namespace Go

#endif
