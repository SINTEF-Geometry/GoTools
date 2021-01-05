//===========================================================================
//                                                                           
// File: LRSpline3DUtils.h                                                   
//                                                                           
// Created: Wed Mar  6 14:51:21 2013                                         
//                                                                           
// Author: Sverre Briseid <Sverre.Briseid@sintef.no>
//                                                                           
// Revision: $Id$
//                                                                           
// Description:
//                                                                           
//===========================================================================

#ifndef _LRSPLINE3DUTILS_H
#define _LRSPLINE3DUTILS_H

#include <array>
#include <functional>
#include <set>
#include <map>
#include <unordered_map>
#include <tuple>
#include <iostream> // @@ debug

#include "GoTools/utils/checks.h"
#include "GoTools/trivariate/SplineVolume.h"
#include "GoTools/lrsplines3D/LRSplineVolume.h"
#include "GoTools/lrsplines3D/Mesh3D.h"
#include "GoTools/lrsplines3D/LRBSpline3D.h"
#include "GoTools/lrsplines3D/Element3D.h"
#include "GoTools/lrsplines3D/LRSplineVolume.h"
#include "GoTools/lrsplines2D/BSplineUniLR.h"


namespace Go
{

  namespace LRSpline3DUtils
  {

      LRSplineVolume::ElementMap identify_elements_from_mesh(const Mesh3D& m);

      void update_elements_with_single_bspline(LRBSpline3D* b, 
					       LRSplineVolume::ElementMap& emap, 
					       const Mesh3D& mesh,
					       bool remove);

      // next_value and prev_value refers to the (x, y, z) orientation
      // (hence prev_value is the next value of next_value).
      int locate_interval(const Mesh3D& m, Direction3D d, double value, 
			  double next_value, double prev_value,
			  bool at_end);


      std::vector<int> set_uniform_meshlines(Direction3D d, Mesh3D& mesh);

      bool all_meshlines_uniform(Direction3D d, const Mesh3D& m);

      std::vector<int> knots_to_insert(const std::vector<int>& ref,
                                       const std::vector<int>& mults);

      double compute_alpha(int degree, const int*
                           const oldvec_ix,
                           const int* const newvec_ix,
                           const double* const kvals);

      std::tuple<std::vector<double>, std::vector<int> >
      insert_knots(const std::vector<int>& new_knots,
                                    std::unique_ptr<LRBSpline3D>& bfun,
                                    const Direction3D d,
                                    const double* const kvals);

      void tensor_split(std::unique_ptr<LRBSpline3D>& bfun,
                        const std::vector<int>& x_mults,
                        const std::vector<int>& y_mults,
                        const std::vector<int>& z_mults,
                        const Mesh3D& tensor_mesh,
			std::vector<std::unique_ptr<BSplineUniLR> >& bspline_vec1,
			std::vector<std::unique_ptr<BSplineUniLR> >& bspline_vec2,
			std::vector<std::unique_ptr<BSplineUniLR> >& bspline_vec3,
                        LRSplineVolume::BSplineMap& bmap);
      LRBSpline3D* 
      insert_basis_function(std::unique_ptr<LRBSpline3D>& b,
			    const Mesh3D& mesh, 
			    LRSplineVolume::BSplineMap& bmap);

      void iteratively_split(std::vector<std::unique_ptr<LRBSpline3D> >& bfuns,
			     const Mesh3D& mesh, 
			     std::vector<std::unique_ptr<BSplineUniLR> >& bspline_vec1,
			     std::vector<std::unique_ptr<BSplineUniLR> >& bspline_vec2,			
			     std::vector<std::unique_ptr<BSplineUniLR> >& bspline_vec3);
      void iteratively_split2(std::vector<LRBSpline3D*>& bfuns,
                              const Mesh3D& mesh,
                              LRSplineVolume::BSplineMap& bmap,
                              double* domain, 
			      std::vector<std::unique_ptr<BSplineUniLR> >& bspline_vec1,
			      std::vector<std::unique_ptr<BSplineUniLR> >& bspline_vec2,			
			      std::vector<std::unique_ptr<BSplineUniLR> >& bspline_vec3,
			      bool support = true);

      std::tuple<int, int, int, int, int, int>
      refine_mesh(Direction3D d, double fixed_val,
		  double start1, double end1,
		  double start2, double end2,
		  int mult, bool absolute,
		  int spline_degree, double knot_tol,
		  Mesh3D& mesh, 
		  std::vector<std::unique_ptr<BSplineUniLR> >& bsplines,
		  bool& refined);
#if 0
      void  zero_knot(Go::Direction3D d, double fixed_val,
		      double knot_tol,Go::Mesh3D& mesh,
		      std::vector<std::unique_ptr<BSplineUniLR> >& bsplines);
#endif      
     void get_affected_bsplines(const std::vector<LRSplineVolume::Refinement3D>& refs, 
				 const LRSplineVolume::ElementMap& emap,
				 double knot_tol, const Mesh3D& mesh,
				 std::vector<LRBSpline3D*>& affected);

      bool support_equal(const LRBSpline3D* b1, const LRBSpline3D* b2);


      SplineVolume* fullTensorProductVolume(const LRSplineVolume& lr_spline_vol);

      // Distribute given data points to elements
      void distributeDataPoints(LRSplineVolume* vol, std::vector<double>& points, 
                                bool add_distance_field = false, 
                                bool primary_points = true);
      
      // Evaluate all basis functions in a point
      void evalAllBSplines(const std::vector<LRBSpline3D*>& bsplines,
			   double upar, double vpar, double wpar,
			   bool u_at_end, bool v_at_end, bool w_at_end,
			   std::vector<double>& result);

      void evalAllBSplines2(const std::vector<LRBSpline3D*>& bsplines,
			   double upar, double vpar, double wpar,
			   bool u_at_end, bool v_at_end, bool w_at_end,
			    double* result, double* val);

      // Evaluate all basis functions in a point and multiplying with the
      // coefficient
      void evalAllBSplinePos(const std::vector<LRBSpline3D*>& bsplines,
			     double upar, double vpar, double wpar,
			     bool u_at_end, bool v_at_end, bool w_at_end,
			     std::vector<Point>& result);

      // Comparison functor which only takes into account the support of the LRBSpline3Ds,
      // not the coefficient or gamma.  This functor is sometimes useful when generating
      // STL structures that keep track of LRBSpline3Ds with unique support.
      //==============================================================================
      struct support_compare
      //==============================================================================
      {
        bool operator()(const LRBSpline3D* b1, const LRBSpline3D* b2) const
          {
            // to compare b1 and b2, compare the x-knotvectors.  If these are identical, compare
            // the y-knotvectors instead.
            const int tmp1 = compare_seq(b1->kvec(XDIR).begin(), b1->kvec(XDIR).end(),
                                         b2->kvec(XDIR).begin(), b2->kvec(XDIR).end());
            if (tmp1 != 0)
                return (tmp1 < 0);
            const int tmp2 = compare_seq(b1->kvec(YDIR).begin(), b1->kvec(YDIR).end(),
                                         b2->kvec(YDIR).begin(), b2->kvec(YDIR).end());
            if (tmp2 != 0)
                return (tmp2 < 0);
            const int tmp3 = compare_seq(b1->kvec(ZDIR).begin(), b1->kvec(ZDIR).end(),
                                         b2->kvec(ZDIR).begin(), b2->kvec(ZDIR).end());
            return (tmp3 < 0);

        }
      };

  } // end namespace LRSplineUtils

}; // end namespace Go

#endif // _LRSPLINE3DUTILS_H

