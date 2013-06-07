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

#ifndef _INTCRVEVALUATOR_
#define _INTCRVEVALUATOR_

#include "GoTools/creators/EvalCurveSet.h"
#include "GoTools/geometry/CurveOnSurface.h"
#include "GoTools/utils/Point.h"

namespace Go
{
  /// This class represents an intersection curve approximated by two
  /// curve on surface instances
  class IntCrvEvaluator : public EvalCurveSet
  {
  public:
    /// Constructor
    IntCrvEvaluator(shared_ptr<CurveOnSurface> sfcv1,
		    double start1, double end1,
		    shared_ptr<CurveOnSurface> sfcv2,
		    double start2, double end2,
		    bool same_orientation,
		    int keep_crv = 0);

    /// Destructor.
    virtual ~IntCrvEvaluator();

    /// Evaluate the curves.
    /// \param t parameter in which to evaluate.
    /// \return the evaluated points for the curve set.
    virtual std::vector<Point> eval(double t);

    /// Evaluate the curve derivatives.
    /// \param t parameter in which to evaluate.
    /// \param n number of derivatives to compute.
    /// \param der the evaluated points up to the n'th derivative for the curve set.
    virtual void eval(double t, int n, std::vector<std::vector<Point> >& der); // n = order of diff

    /// Start parameter of domain.
    /// \return start parameter of the spline space.
    virtual double start();

    /// End parameter of domain.
    /// \return end parameter of the spline space.
    virtual double end();

    /// The geometric dimension of the spline curves.
    /// \return geometric dimension of the total space.
    virtual int dim();

    /// Whether the approximation is within tolerances in input parameter.
    /// \param par parameter in which to evaluate.
    /// \param approxpos whether the input points are within tolerance from the
    ///                  evaluated points (as given by eval()).
    /// \param tol1 tolerance used to decide approximation accuracy.
    /// \param tol2 tolerance used to decide approximation accuracy.
    /// \return whether the approximation is within tolerances in input parameter.
    virtual bool approximationOK(double par, const std::vector<Point>& approxpos,
				 double tol1, double tol2);

    /// The number of curves in the curve set.
    /// \return the number of curves in the curve set.
    virtual int nmbCvs();

    /// Reset intermediate error. New iterationa
    virtual void resetErr()
    {
      max_err_ = 0.0;
    }

    /// Get maximum approximation error
    double getMaxErr()
    {
      return std::max(max_dist_, max_err_);
    }

  private:
    shared_ptr<CurveOnSurface> sfcv1_;
    double start1_;
    double end1_;
    shared_ptr<CurveOnSurface> sfcv2_;
    double start2_;
    double end2_;
    bool same_orientation_;
    int keep_crv_;
    double max_dist_;
    double max_err_;

    void evaluate(double t, int n, std::vector<Point>& result);
  };

} // namespace Go

#endif //_INTCRVEVALUATOR_
