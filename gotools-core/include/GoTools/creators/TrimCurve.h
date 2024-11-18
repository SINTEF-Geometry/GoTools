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
  /// limited in the geometry space and in parameter space. This can avoid
  /// troubles in the case of a closed input curve. It is left to the user
  /// to ensure consistence
  TrimCurve(double startpar, Point startpt, double endpar, Point endpt,
	    CurveOnSurface* bd_crv);

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
