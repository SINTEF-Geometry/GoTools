//===========================================================================
//                                                                           
// File: SmoothTransition.h                                                
//                                                                           
// Created: Mon Dec  9 11:21:48 2002                                         
//                                                                           
// Author: Sverre Briseid <Sverre.Briseid@sintef.no>
//                                                                           
// Revision: $Id: SmoothTransition.h,v 1.5 2007-12-04 16:11:39 jbt Exp $
//                                                                           
// Description: Given input of curve lying on two surfaces, we create a smooth transition
//              surface by offsetting both surfaces, approximating intersecting curve,
//              continue with projection of intersection curve on the two surfaces.
//              Resulting projected curves and cross tangent curves may then be lofted.
//              We're thus computing four points in eval().
//                                                                           
//===========================================================================


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

    /// Constructor. By not using CurveOnSurface, object is more general.
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

    void blend_s1421(const SplineSurface* ps, double aoffset, int ider,
		     const Point& epar, int& ilfs, int& ilft,
		     std::vector<Point>& eoffpnt, std::vector<Point>& epnt, int* jstat);


};

} // namespace Go

#endif // _SMOOTHTRANSITION_H_
