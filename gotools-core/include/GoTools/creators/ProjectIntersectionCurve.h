
#ifndef _PROJECTCURVE_
#define _PROJECTCURVE_

#include <memory>

#include "GoTools/utils/Point.h"
#include "GoTools/creators/EvalCurve.h"
#include "GoTools/geometry/SplineCurve.h"
#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/utils/config.h"

namespace Go 

{

/// This class provides an interface to a curve that can be evaluated.
/// This evaluator based class computes the projected point in the
/// first surface of the offset point defined by input.
class ProjectIntersectionCurve : public EvalCurve
{
public:

    /// Constructor.
    /// \param inters_crv the space intersection curve between surf & other_surf.
    /// \param p_crv the corresponding parameter curve in surf.
    /// \param other_p_crv the corresponding parameter curve in other_surf.
    /// \param surf the first input surface.
    /// \param other_surf the second input surface.
    /// \param offset_dist the offset distance in surf.
    /// \param other_offset_dist the offset distance in other_surf.
    /// \param epsgeo the geometrical tolerance (for closes point evaluations).
    ProjectIntersectionCurve(shared_ptr<SplineCurve>& inters_crv,
			     shared_ptr<SplineCurve>& p_crv,
			     shared_ptr<SplineCurve>& other_p_crv,
			     shared_ptr<ParamSurface>& surf,
			     shared_ptr<ParamSurface>& other_surf,
			     double offset_dist, double other_offset_dist,
			     double epsgeo);

    /// Destructor.
    virtual ~ProjectIntersectionCurve();

    /// The evaluator part of the class, returns the projected offset point in
    /// the first surface.
    /// \param t the parameter in which to evaluate.
    virtual Point eval(double t) const;

    /// The evaluator part of the class, returns the projected offset point in
    /// the first surface.  For n == 1 the tangent in the projected offset curve
    /// is also computed.
    /// \param t the parameter in which to evaluate.
    /// \param n the number of derivatives to compute (at most 1).
    /// \param der the evaluated point.  Size of array is 'n + 1'.
    virtual void eval(double t, int n, Point der[]) const;

    /// Start parameter of curve.
    /// \return the start parameter.
    virtual double start() const;
    /// End parameter of curve.
    /// \return the end parameter.
    virtual double end() const;

    /// Dimension of inters_crv_.
    /// \return the dimension of the evaluator point.
    virtual int dim() const;

    /// Whether the evaluated point in par is close enough to approxpos.
    /// \param par the parameter in which to evaluate.
    /// \param approxpos postition to check for accuracy.
    /// \param tol1 currently not used.
    /// \param tol2 currently not used.
    /// \return whether approxpos is within satisfactory accuracy (i.e. epsgeo_).
    virtual bool approximationOK(double par, Point approxpos,
				 double tol1, double tol2) const; 

private:
    const shared_ptr<SplineCurve> inters_crv_;
    //Param curves serve as seed generators for closest point evaluations.
    const shared_ptr<SplineCurve> p_crv_;
    const shared_ptr<SplineCurve> other_p_crv_;
    const shared_ptr<ParamSurface> surf_;
    const shared_ptr<ParamSurface> other_surf_;
    const double offset_dist_; // In direction normal to surf_.
    const double other_offset_dist_; // In direction normal to other_surf_.
    const double epsgeo_;

    /// Compute the parameter point in t.  If pcv_turned the returned point takes this into
    /// account.
    /// \param space_cv the curve defining the parametrization.
    /// \param t the parameter in which to evaluate.
    /// \param_cv the curve in which to evaluate.
    /// \param pcv_turned whether the parametrization of param_cv is the opposite of space_cv.
    /// return the corresponding parameter point in param_cv.
    std::vector<double>
    getSuggestedSurfaceParameter(const SplineCurve& space_cv, double t,
				 const SplineCurve& param_cv,
				 bool pcv_turned) const;

};

} // namespace Go

#endif //_PROJECTCURVE_
