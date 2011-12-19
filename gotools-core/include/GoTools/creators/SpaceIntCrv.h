#ifndef _SPACEINTCRV_
#define _SPACEINTCRV_

#include "GoTools/creators/EvalCurve.h"
#include "GoTools/geometry/CurveOnSurface.h"
#include "GoTools/utils/Point.h"
#include "GoTools/utils/config.h"

namespace Go
{
  /// This class represents an intersection curve approximated by two
  /// curve on surface instances
  class SpaceIntCrv : public EvalCurve
  {
  public:
    /// Constructor
    SpaceIntCrv(shared_ptr<ParamCurve> init_crv, int pardir,
		std::vector<shared_ptr<CurveOnSurface> >& sfcv1, 
		std::vector<double> start1, 
		std::vector<double> end1,
		std::vector<shared_ptr<CurveOnSurface> >& sfcv2,
		std::vector<double> start2, std::vector<double> end2,
		std::vector<bool> opposite, bool same_orient);

    /// Destructor.
    virtual ~SpaceIntCrv();

    /// Evaluate the curves.
    /// \param t parameter in which to evaluate.
    /// \return the evaluated point for the curve.
    virtual Point eval(double t) const;

    /// Evaluate the curve derivatives.
    /// \param t parameter in which to evaluate.
    /// \param n number of derivatives to compute.
    /// \param der the evaluated points up to the n'th derivative for the curve.
    virtual void eval(double t, int n, Point der[]) const; // n = order of diff

    /// Start parameter of domain.
    /// \return start parameter of the spline space.
    virtual double start() const;

    /// End parameter of domain.
    /// \return end parameter of the spline space.
    virtual double end() const;

    /// The geometric dimension of the spline curves.
    virtual int dim() const;

    /// Whether the approximation is within tolerances in input parameter.
    /// \param par parameter in which to evaluate.
    /// \param approxpos whether the input point are within tolerance from the
    ///                  evaluated points (as given by eval()).
    /// \param tol1 tolerance used to decide approximation accuracy.
    /// \param tol2 tolerance used to decide approximation accuracy.
    /// \return whether the approximation is within tolerances in input parameter.
    virtual bool approximationOK(double par, Point approxpos,
				 double tol1, double tol2) const;

  private:
    shared_ptr<ParamCurve> init_crv_;
    std::vector<shared_ptr<CurveOnSurface> > sfcv1_;
    std::vector<shared_ptr<CurveOnSurface> > sfcv2_;
    std::vector<double> start1_, end1_, start2_, end2_;
    std::vector<double> segment_;
    std::vector<bool> opposite_;
    bool same_orient_;

    void evaluate(double t, int n, Point result[]) const;
  };

} // namespace Go

#endif //_SPACEINTCRV_
