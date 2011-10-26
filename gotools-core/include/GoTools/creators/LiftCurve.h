//===========================================================================
//                                                                           
// File: LiftCurve.C                                                    
//                                                                           
// Created: Sat Mar 29 09:00:22 2003                                         
//                                                                           
// Author: Sverre Briseid <Sverre.Briseid@sintef.no>
//                                                                           
// Revision: $Id: LiftCurve.h,v 1.2 2007-12-04 16:11:39 jbt Exp $
//                                                                           
// Description:
//                                                                           
//===========================================================================

#ifndef _LIFTCURVE_
#define _LIFTCURVE_



#include "GoTools/utils/Point.h"
#include "GoTools/creators/EvalCurve.h"
#include "GoTools/geometry/SplineCurve.h"
#include "GoTools/geometry/SplineSurface.h"

#include <memory>

namespace Go
{

/// This class represents a "lift curve", generated from taking
/// a surface and a 2D curve and then evaluation the surface 
/// at the parameter values obtained by evaluating the 2D curve.

class LiftCurve : public EvalCurve
{
public:

  /// Constructor, taking a 2D parameter curve and a surface.
  /// \param parameter_crv the 2D parameter curve that will be 'lifted'
  /// \param surf the surface on which the resulting, 'lifted' curve will
  ///             lie
  /// \param epsgeo geometrical tolerance used when running the 'approximationOK'
  ///               function.
  LiftCurve(std::shared_ptr<Go::ParamCurve>& parameter_crv,
	      std::shared_ptr<Go::ParamSurface>& surf,
	      double epsgeo);

  /// virtual destructor enables safe inheritance
  virtual ~LiftCurve();

  // Inherited from EvalCurve
  virtual Point eval( double t) const;

  // Inherited from EvalCurve
  virtual void eval(double t, int n, Point der[]) const;

  // Inherited from EvalCurve
  virtual double start() const;

  // Inherited from EvalCurve
  virtual double end() const;

  /// Dimension of the lifted curve (i.e. 3).
  virtual int dim() const;

  /// Inherited from EvalCurve::approximationOK().  
  /// \param par the parameter at which to check the curve
  /// \param approxpos the position we want to check whether or not the curve
  ///                  approximates for parameter 'par'.
  /// \param tol1 unused
  /// \param tol2 unused
  /// \return 'true' if the curve approximates the point at the parameter
  ///         (within the tolerance given in the constructor, 'epsgeo'). 'false'
  ///         otherwise.
  virtual bool approximationOK(double par, Point approxpos,
			       double tol1, double tol2) const;

 private:
  const std::shared_ptr<Go::ParamCurve> parameter_crv_;
  const std::shared_ptr<Go::ParamSurface> surf_;
  const double epsgeo_;

};


} // namespace Go

#endif //_LIFTCURVE_
