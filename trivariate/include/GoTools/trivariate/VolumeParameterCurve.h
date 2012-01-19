#ifndef _VOLUMEPARAMETERCURVE_H_
#define _VOLUMEPARAMETERCURVE_H_

#include <memory>
#include "GoTools/utils/Point.h"
#include "GoTools/creators/EvalCurve.h"


namespace Go
{

  class ParamVolume;
  class ParamCurve;

  /// \brief An evaluator based curve representing the parameter domain curve
  /// in a given volume which represents the same curve as a given space curve.
  /// Project a geometry curve into the parameter domain of a volume.

  class VolumeParameterCurve : public EvalCurve
{
public:
  /// Constructor
  /// \param vol volume in which the curve lies
  /// \param crv space curve representing the same curve
  VolumeParameterCurve(shared_ptr<ParamVolume> vol, 
		       shared_ptr<ParamCurve> crv);

  /// Constructor
  /// \param vol volume in which the curve lies
  /// \param crv space curve representing the same curve
  /// \param par1 startparameter of the final curve par1 >= crv->startparam()
  /// \param par2 endparameter of the final curve par1 <= crv->endparam()
  VolumeParameterCurve(shared_ptr<ParamVolume> vol, 
		       shared_ptr<ParamCurve> crv,
		       shared_ptr<Point> par1,
		       shared_ptr<Point> par2);

  /// virtual destructor ensures save inheritance
  virtual ~VolumeParameterCurve();

  /// Evaluate a point on the curve for a given parameter
  /// \param t the parameter for which to evaluate the curve.
  /// \return the evaluated point
  virtual Point eval( double t) const;

  /// Evaluate a point and a certain number of derivatives 
  /// on the curve for a given parameter.
  /// \param t the parameter for which to evaluate the curve.
  /// \param n the number of derivatives (0 or more)
  /// \retval der pointer to an array of Points where the 
  ///         result will be written.  The position will be stored
  ///         first, then the first derivative (tangent), then the
  ///         second, etc..
  ///         \b NB: For most (all) derived classes of 'EvalCurve', 
  ///         the implementation actually only supports the computation of 
  ///         one derivative, i.e. if n > 1, only one derivative will be 
  ///         computed anyway.
  virtual void eval( double t, int n, Point der[]) const; // n = order of diff

  /// Get the start parameter of the curve.
  /// \return the start parameter of the curve.
  virtual double start() const;
  
  /// Get the end parameter of the curve.
  /// \return  the end parameter of the curve.
  virtual double end() const;

  /// Get the dimension of the space in which the curve lies.
  /// \return the space dimension of the curve.
  virtual int dim() const;

  /// Check if the curve, evaluated at a given parameter, approximates
  /// a given position within a given tolerance.
  /// \param par the parameter at which to check the curve
  /// \param approxpos the position we want to check whether or not the curve
  ///                  approximates for parameter 'par'.
  /// \param tol1 approximation tolerance.
  /// \param tol2 another approximation tolerance (its use is defined by some of
  ///             the derived classes.
  /// \return 'true' if the curve approximates the point at the parameter, 'false'
  ///         otherwise.
  virtual bool approximationOK(double par, Point approxpos,
			       double tol1, double tol2) const;

 private:
  shared_ptr<ParamVolume> vol_;
  shared_ptr<ParamCurve> crv_;
  shared_ptr<Point> par1_;
  shared_ptr<Point> par2_;
};

}

#endif
