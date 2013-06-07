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

#ifndef _SPACEINTCRV_
#define _SPACEINTCRV_

#include "GoTools/creators/EvalCurve.h"
#include "GoTools/geometry/CurveOnSurface.h"
#include "GoTools/utils/Point.h"
#include "GoTools/utils/config.h"

namespace Go
{
  /// This class represents an intersection curve approximated by two
  /// curve on surface instances
  class SpaceIntCrv : public EvalCurve
  {
  public:
    /// Constructor
    /// \param init_crv approximation to true intersection curve
    /// \param pardir assumes monotonicity of intersection curve in this
    /// parameter direction
    /// \param sfcv1 set of curves joining into the intersection curve representation 
    /// related to first surface
    /// \param start1 startparameter of the curves sfcv1. May restrict the extent of each curve.
    /// \param end1 endparameter of the curves sfcv1. May restrict the extent of each curve.
    /// \param sfcv2 set of curves joining into intersection curve representation related to second surface
    /// \param start2 startparameter of the curves sfcv2. May restrict the extent of each curve.
    /// \param end2 endparameter of the curves sfcv2. May restrict the extent of each curve.
    /// \param opposite pairwise comparisment of orientation from the curves in sfcv1 and sfcv2
    /// \param same_orient whether the curves in sfcv1 have the same orientation
    /// as the underlying surface
    SpaceIntCrv(shared_ptr<ParamCurve> init_crv, int pardir,
		std::vector<shared_ptr<CurveOnSurface> >& sfcv1, 
		std::vector<double> start1, 
		std::vector<double> end1,
		std::vector<shared_ptr<CurveOnSurface> >& sfcv2,
		std::vector<double> start2, std::vector<double> end2,
		std::vector<bool> opposite, bool same_orient);

    /// Destructor.
    virtual ~SpaceIntCrv();

    /// Evaluate the curves.
    /// \param t parameter in which to evaluate.
    /// \return the evaluated point for the curve.
    virtual Point eval(double t) const;

    /// Evaluate the curve derivatives.
    /// \param t parameter in which to evaluate.
    /// \param n number of derivatives to compute.
    /// \param der the evaluated points up to the n'th derivative for the curve.
    virtual void eval(double t, int n, Point der[]) const; // n = order of diff

    /// Start parameter of domain.
    /// \return start parameter of the spline space.
    virtual double start() const;

    /// End parameter of domain.
    /// \return end parameter of the spline space.
    virtual double end() const;

    /// The geometric dimension of the spline curves.
    virtual int dim() const;

    /// Whether the approximation is within tolerances in input parameter.
    /// \param par parameter in which to evaluate.
    /// \param approxpos whether the input point are within tolerance from the
    ///                  evaluated points (as given by eval()).
    /// \param tol1 tolerance used to decide approximation accuracy.
    /// \param tol2 tolerance used to decide approximation accuracy.
    /// \return whether the approximation is within tolerances in input parameter.
    virtual bool approximationOK(double par, Point approxpos,
				 double tol1, double tol2) const;

  private:
    shared_ptr<ParamCurve> init_crv_;
    std::vector<shared_ptr<CurveOnSurface> > sfcv1_;
    std::vector<shared_ptr<CurveOnSurface> > sfcv2_;
    std::vector<double> start1_, end1_, start2_, end2_;
    std::vector<double> segment_;
    std::vector<bool> opposite_;
    bool same_orient_;

    void evaluate(double t, int n, Point result[]) const;
  };

} // namespace Go

#endif //_SPACEINTCRV_
