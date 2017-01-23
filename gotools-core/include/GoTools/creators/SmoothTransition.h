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

#ifndef _SMOOTHTRANSITION_H_
#define _SMOOTHTRANSITION_H_

#include <memory>

#include "GoTools/utils/Point.h"
#include "GoTools/creators/EvalCurveSet.h"
#include "GoTools/geometry/SplineCurve.h"
#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/geometry/SISLconversion.h"

namespace Go
{

    /// This abstract class provides an interface to a curve that can be evaluated.
    /// Given input of curve lying on two surfaces, we create a smooth transition
    /// surface by offsetting both surfaces, approximating intersecting curve,
    /// continue with projection of intersection curve on the two surfaces.
    /// Resulting projected curves and cross tangent curves may then be lofted.
    /// We're thus computing four points in eval().

class SmoothTransition : public EvalCurveSet
{
public:

  /// Constructor. 
  /// \param inters_crv intersection curve between two surfaces
  /// \param p_crv1 associated parameter curve in first surface
  /// \param p_crv2 associated parameter curve in second surface
  /// \param surf1 first surface
  /// \param surf2 second surface
  /// \param offset_dist1 offset distance related to first surface
  /// \param offset_dist2 offset distance related to second surface
  /// \param epsgeo approximation tolerance and tolerance used in intersection
  /// computations
  // By not using CurveOnSurface, object is more general.
    SmoothTransition(shared_ptr<const SplineCurve>& inters_crv,
		     shared_ptr<const SplineCurve>& p_crv1,
		     shared_ptr<const SplineCurve>& p_crv2,
		     shared_ptr<const ParamSurface> surf1,
		     shared_ptr<const ParamSurface> surf2,
		     double offset_dist1, double offset_dist2,
		     double epsgeo);

    // Inherited from EvalCurve

    /// Empty destructor.
    virtual ~SmoothTransition();
    virtual std::vector<Point> eval( double t); // offset1, offset1p, offset1_cross_tan,
    // offset2, offset2p, offset2_cross_tan
    virtual void eval(double t, int n, std::vector<std::vector<Point> >& der); // Exact value.
    virtual double start();
    virtual double end();
    virtual int dim(); // Dimension of space, i.e. 3.
    virtual bool approximationOK(double par, const std::vector<Point>& approxpos,
				 double tol1, double tol2);
    virtual int nmbCvs()
    { return 6; }

private:
    shared_ptr<const SplineCurve> inters_crv_;
    // Param curves serve as seed generators for closest point eval.
    shared_ptr<const SplineCurve> p_crv1_;
    shared_ptr<const SplineCurve> p_crv2_;
    shared_ptr<const ParamSurface> surf1_;
    shared_ptr<const ParamSurface> surf2_;
    shared_ptr<const SplineSurface> under_surf1_;
    shared_ptr<const SplineSurface> under_surf2_;
    double offset_dist1_; // In direction normal to surf1_.
    double offset_dist2_; // In direction normal to surf2_.
    const double epsgeo_;
    const double kinktol_;
    std::vector<double> tangent_lengths_; // We set required lengths on the tangents (0, 1, 3, 4).

    // Given space point, we project onto surface, returning parameter values.
    // If seed has size two, value is used in closest point evaluation.
    Point projectPoint(const Point& space_pt, const ParamSurface& surf,
		       std::vector<double>& seed, bool boundary_pt,
		       double epsgeo, double& dist);

    // Given input we compute the point and tangent std::vector of the cross tangent curve.
    // space_pt must be of size 2, local_pt of size derivs+1 (derivs not larger than 1).
    std::vector<Point> computeCrosstangentValues(std::vector<Point>& space_pt,
						 std::vector<Point>& local_pt, int derivs);

    // We try to guess parameter values of intersection between sf1 & sf2 & plane defined by
    // inters_cv_pt and it's tangent. Tangent is not needed as it is defined by normal in sfs.
    void guessParameterPoints(const Point& inters_cv_pt, double t,
			      const SplineCurve& inters_cv,
			      const ParamSurface& sf1, const ParamSurface& sf2,
			      const SplineCurve& p_inters_cv1, const SplineCurve& p_inters_cv2,
			      double offset_dist1, double offset_dist2,
			      Point& guess_pt1, Point& guess_pt2);

    // param_cv may be parametrized in the opposite direction, as given by pcv_turned.
    std::vector<double>
    getSuggestedSurfaceParameter(const SplineCurve& space_cv, double t,
				 const SplineCurve& param_cv,
				 bool pcv_turned);

    // Given a point in space, close to input surface point, we use partial derivatives in input
    // point to make a guess on parameter values of projection of space point.
    // surf_par_pt is of dimension 2, while space_pt shares dimension with surf.
    std::vector<double> getSuggestedSurfaceParameter(Point& surf_par_pt, const ParamSurface& surf,
						     Point& space_pt, double tolerance);

    void offsetIntersectionPoints(std::vector<Point>& ep, std::vector<Point>& eq,
				  std::vector<Point>& eoffp, std::vector<Point>& eoffq,
				  Point& eparp, Point& eparq,
				  std::vector<Point>& espine, std::vector<Point>& egeobb1,
				  std::vector<Point>& egeobb2, std::vector<Point>& ecrtan1,
				  std::vector<Point>& ecrtan2, std::vector<Point>& egeop,
				  std::vector<Point>& egeoq, std::vector<double>& curv_radis);

    void offsetIntersectionIterate(double arad1, double arad2, std::vector<Point>& epoint,
				   std::vector<Point>& epnt1, std::vector<Point>& epnt2,
				   Point& epar1, Point& epar2,
				   const SplineSurface& psurf1, const SplineSurface& psurf2,
				   double astep, double aepsge, std::vector<Point>& gpnt1,
				   std::vector<Point>& gpnt2, std::vector<Point>& goffpnt1,
				   std::vector<Point>& goffpnt2, Point& gpar1,
				   Point& gpar2p);


};

} // namespace Go

#endif // _SMOOTHTRANSITION_H_
