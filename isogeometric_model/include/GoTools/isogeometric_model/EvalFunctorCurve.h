//===========================================================================
//
// File : EvalFunctorCurve.h
//
// Created: Mon Oct 11 15:47:47 2010
//
// Author: Kjell Fredrik Pettersen
//
// Revision: $Id:$
//
// Description:
//
//===========================================================================



#ifndef __EVALFUNCTORCURVE_H
#define __EVALFUNCTORCURVE_H


#include <memory>
#include "GoTools/creators/EvalCurve.h"
#include "GoTools/geometry/SplineCurve.h"
#include "GoTools/isogeometric_model/BdCondFunctor.h"


namespace Go
{


  class EvalFunctorCurve : public EvalCurve
  {

  public:

    // Constructor
    EvalFunctorCurve(BdCondFunctor* fbd, shared_ptr<SplineCurve> geo_curve, int dimension);

    // Destructor
    virtual ~EvalFunctorCurve();

    // Inherited functions from EvalCurve

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

    // The boundary condition functor
    BdCondFunctor* functor_;

    // The spline curve in the geometry space
    shared_ptr<SplineCurve> geo_curve_;

    // The dimension of the solution space
    int dim_;

  };    // Class EvalFunctorCurve


} // namespace Go

#endif    // #ifndef __EVALFUNCTORCURVE_H
