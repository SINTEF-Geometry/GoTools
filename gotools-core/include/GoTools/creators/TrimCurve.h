#ifndef _TRIMCURVE_
#define _TRIMCURVE_

#include <memory>

#include "GoTools/utils/Point.h"
#include "GoTools/creators/EvalCurveSet.h"
#include "GoTools/geometry/SplineCurve.h"
#include "GoTools/geometry/SplineSurface.h"


namespace Go
{

  class CurveOnSurface;
/// This class represents the curve obtained by projecting a 
/// given 3D curve onto a given part of a given 3D surface.
/// Used to improve the accuracy of an already existing trimming curve.

class TrimCurve : public EvalCurveSet
{
public:

  /// Constructor given the CurveOnSurface curve representing the trim curve
  TrimCurve(CurveOnSurface* bd_crv);

  /// Constructor given the CurveOnSurface curve representing the trim curve
  /// limited in the parameter domain
  TrimCurve(CurveOnSurface* bd_crv, double start, double end);

  /// Constructor given the CurveOnSurface curve representing the trim curve
  /// limited in the geometry space
  TrimCurve(Point startpt, Point endpt, CurveOnSurface* bd_crv);

    /// virtual destructor ensures safe inheritance
    virtual ~TrimCurve();
    
    // Inherited from EvalCurveSet
    std::vector<Point> eval( double t);

    // Inherited from EvalCurveSet
    virtual void eval(double t, int n, std::vector<std::vector<Point> >& der);

    // Inherited from EvalCurveSet
    virtual double start();

    // Inherited from EvalCurveSet
    virtual double end();

    /// Inherited from EvalCurveSet::dim().  
    virtual int dim();

    /// Inherited from EvalCurveSet::approximationOK().  For this class, the specified tolerances
    /// are not used; the internally stored 'epsgeo' value is used as tolerance (this value was
    /// specified in the constructor).
    /// \param par the parameter at which to check the curve
    /// \param approxpos the position we want to check whether or not the curve
    ///                  approximates for parameter 'par'.
    /// \param tol1 unused
    /// \param tol2 unused
    /// \return 'true' if the curve approximates the point at the parameter, 'false'
    ///         otherwise.
    virtual bool approximationOK(double par, const std::vector<Point>& approxpos,
				 double tol1, double tol2); 

    /// The number of curves in the curve set.
    /// \return the number of curves in the curve set.
    virtual int nmbCvs();

private:
    CurveOnSurface* sfcv_;
    Point startpt_;
    Point endpt_;
    double start_;
    double end_;

    void evaluate(double t, int n, std::vector<Point>& result);
};


}

#endif //_TRIMCURVE_
