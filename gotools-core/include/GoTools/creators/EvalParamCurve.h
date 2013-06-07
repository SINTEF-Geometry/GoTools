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
