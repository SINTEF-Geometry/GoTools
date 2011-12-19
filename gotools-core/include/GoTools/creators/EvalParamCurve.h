//===========================================================================
//                                                                           
// File: EvalParamCurve.C                                                    
//                                                                           
// Created: September 2010
//                                                                           
// Author: Vibeke Skytt
//                                                                           
// Revision: $Id: EvalParamCurve.h
//                                                                           
// Description:
//                                                                           
//===========================================================================

#ifndef _EVALPARAMCURVE_
#define _EVALPARAMCURVE_



#include "GoTools/utils/Point.h"
#include "GoTools/creators/EvalCurve.h"
#include "GoTools/geometry/SplineCurve.h"
#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/utils/config.h"

#include <memory>

namespace Go
{

/// This class represents a interface to a ParamCurve to make
/// it fit into the evaluator based curve concept

class EvalParamCurve : public EvalCurve
{
public:

  /// Constructor, taking a parametric curve
  EvalParamCurve(shared_ptr<Go::ParamCurve>& crv);

  /// virtual destructor enables safe inheritance
  virtual ~EvalParamCurve();

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

  // Debug
  virtual void write(std::ostream& out) const;
	

 private:
  const shared_ptr<Go::ParamCurve> crv_;

};


} // namespace Go

#endif //_EVALPARAMCURVE_
