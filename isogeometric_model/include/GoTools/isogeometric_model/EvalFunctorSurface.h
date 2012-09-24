//===========================================================================
//                                                                           
// File: EvalFunctorSurface.h                                                
//                                                                           
// Created: Wed Sep 19 14:56:05 2012                                         
//                                                                           
// Author: Sverre Briseid <Sverre.Briseid@sintef.no>
//                                                                           
// Revision: $Id$
//                                                                           
// Description:
//                                                                           
//===========================================================================

#ifndef _EVALFUNCTORSURFACE_H
#define _EVALFUNCTORSURFACE_H


#include <memory>
#include "GoTools/creators/EvalSurface.h"
#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/isogeometric_model/BdCondFunctor.h"


namespace Go
{


  class EvalFunctorSurface : public EvalSurface
  {

  public:

    // Constructor
    EvalFunctorSurface(BdCondFunctor* fbd, shared_ptr<SplineSurface> geo_surface, int dimension);

    // Destructor
    virtual ~EvalFunctorSurface();

    // Inherited functions from EvalSurface

    /// Evaluate a point on the surface for a given parameter
    /// \param t the parameter for which to evaluate the surface.
    /// \return the evaluated point
      virtual Point eval( double u, double v) const;

    /// Evaluate a point and a certain number of derivatives 
    /// on the surface for a given parameter.
    /// \param t the parameter for which to evaluate the surface.
    /// \param n the number of derivatives (0 or more)
    /// \retval der pointer to an array of Points where the 
    ///         result will be written.  The position will be stored
    ///         first, then the first derivative (tangent), then the
    ///         second, etc..
    ///         \b NB: For most (all) derived classes of 'EvalSurface', 
    ///         the implementation actually only supports the computation of 
    ///         one derivative, i.e. if n > 1, only one derivative will be 
    ///         computed anyway.
      virtual void eval( double u, double v, int n, Point der[]) const; // n = order of diff

    /// Get the start parameter of the surface.
    /// \return the start parameter of the surface.
    virtual double start_u() const;
    virtual double start_v() const;
  
    /// Get the end parameter of the surface.
    /// \return  the end parameter of the surface.
    virtual double end_u() const;
    virtual double end_v() const;

    /// Get the dimension of the space in which the surface lies.
    /// \return the space dimension of the surface.
    virtual int dim() const;

    /// Check if the surface, evaluated at a given parameter, approximates
    /// a given position within a given tolerance.
    /// \param par the parameter at which to check the surface
    /// \param approxpos the position we want to check whether or not the surface
    ///                  approximates for parameter 'par'.
    /// \param tol1 approximation tolerance.
    /// \param tol2 another approximation tolerance (its use is defined by some of
    ///             the derived classes.
    /// \return 'true' if the surface approximates the point at the parameter, 'false'
    ///         otherwise.
      virtual bool approximationOK(double par_u, double par_v, Point approxpos,
				 double tol1, double tol2) const;
    

  private:

    // The boundary condition functor
    BdCondFunctor* functor_;

    // The spline surface in the geometry space
    shared_ptr<SplineSurface> geo_surface_;

    // The dimension of the solution space
    int dim_;

  };    // Class EvalFunctorSurface


} // namespace Go


#endif // _EVALFUNCTORSURFACE_H

