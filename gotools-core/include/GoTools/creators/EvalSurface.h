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

#ifndef _EVALSURFACE_H
#define _EVALSURFACE_H


#include "GoTools/utils/Point.h"
#include <iostream>



namespace Go
{

/// This is the abstract base class of a curve type that can be evaluated.
/// Representing the exact geometry, typically used when iteratively
/// approximating the curve.

class EvalSurface
{
public:
  /// virtual destructor ensures save inheritance
  virtual ~EvalSurface();

  /// Evaluate a point on the surface for a given parameter
  /// \param t the parameter for which to evaluate the curve.
  /// \return the evaluated point
    virtual Point eval( double u, double v) const = 0;

  /// Evaluate a point and a certain number of derivatives 
  /// on the surface for a given parameter.
  /// \param t the parameter for which to evaluate the curve.
  /// \param n the number of derivatives (0 or more)
  /// \retval der pointer to an array of Points where the 
  ///         result will be written.  The position will be stored
  ///         first, then the first derivative (tangent), then the
  ///         second, etc..
  ///         \b NB: For most (all) derived classes of 'EvalSurface', 
  ///         the implementation actually only supports the computation of 
  ///         one derivative, i.e. if n > 1, only one derivative will be 
  ///         computed anyway.
  virtual void eval( double u, double v, int n, Point der[]) const = 0; // n = order of diff

  /// Get the start parameter of the curve.
  /// \return the start parameter of the curve.
  virtual double start_u() const =0;
  virtual double start_v() const =0;
  
  /// Get the end parameter of the curve.
  /// \return  the end parameter of the curve.
  virtual double end_u() const =0;
  virtual double end_v() const =0;

  /// Get the dimension of the space in which the curve lies.
  /// \return the space dimension of the curve.
  virtual int dim() const = 0;

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
  virtual bool approximationOK(double par_u, double par_v, Point approxpos,
			       double tol1, double tol2) const = 0;

  // Debug
  virtual void write(std::ostream& out) const;

};

}

#endif // _EVALSURFACE_H

