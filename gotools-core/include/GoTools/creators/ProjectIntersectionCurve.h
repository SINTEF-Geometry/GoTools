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

#ifndef _PROJECTCURVE_
#define _PROJECTCURVE_

#include <memory>

#include "GoTools/utils/Point.h"
#include "GoTools/creators/EvalCurve.h"
#include "GoTools/geometry/SplineCurve.h"
#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/utils/config.h"

namespace Go 

{

/// This class provides an interface to a curve that can be evaluated.
/// This evaluator based class computes the projected point in the
/// first surface of the offset point defined by input.
class ProjectIntersectionCurve : public EvalCurve
{
public:

    /// Constructor.
    /// \param inters_crv the space intersection curve between surf & other_surf.
    /// \param p_crv the corresponding parameter curve in surf.
    /// \param other_p_crv the corresponding parameter curve in other_surf.
    /// \param surf the first input surface.
    /// \param other_surf the second input surface.
    /// \param offset_dist the offset distance in surf.
    /// \param other_offset_dist the offset distance in other_surf.
    /// \param epsgeo the geometrical tolerance (for closes point evaluations).
    ProjectIntersectionCurve(shared_ptr<SplineCurve>& inters_crv,
			     shared_ptr<SplineCurve>& p_crv,
			     shared_ptr<SplineCurve>& other_p_crv,
			     shared_ptr<ParamSurface>& surf,
			     shared_ptr<ParamSurface>& other_surf,
			     double offset_dist, double other_offset_dist,
			     double epsgeo);

    /// Destructor.
    virtual ~ProjectIntersectionCurve();

    /// The evaluator part of the class, returns the projected offset point in
    /// the first surface.
    /// \param t the parameter in which to evaluate.
    virtual Point eval(double t) const;

    /// The evaluator part of the class, returns the projected offset point in
    /// the first surface.  For n == 1 the tangent in the projected offset curve
    /// is also computed.
    /// \param t the parameter in which to evaluate.
    /// \param n the number of derivatives to compute (at most 1).
    /// \param der the evaluated point.  Size of array is 'n + 1'.
    virtual void eval(double t, int n, Point der[]) const;

    /// Start parameter of curve.
    /// \return the start parameter.
    virtual double start() const;
    /// End parameter of curve.
    /// \return the end parameter.
    virtual double end() const;

    /// Dimension of inters_crv_.
    /// \return the dimension of the evaluator point.
    virtual int dim() const;

    /// Whether the evaluated point in par is close enough to approxpos.
    /// \param par the parameter in which to evaluate.
    /// \param approxpos postition to check for accuracy.
    /// \param tol1 currently not used.
    /// \param tol2 currently not used.
    /// \return whether approxpos is within satisfactory accuracy (i.e. epsgeo_).
    virtual bool approximationOK(double par, Point approxpos,
				 double tol1, double tol2) const; 

private:
    const shared_ptr<SplineCurve> inters_crv_;
    //Param curves serve as seed generators for closest point evaluations.
    const shared_ptr<SplineCurve> p_crv_;
    const shared_ptr<SplineCurve> other_p_crv_;
    const shared_ptr<ParamSurface> surf_;
    const shared_ptr<ParamSurface> other_surf_;
    const double offset_dist_; // In direction normal to surf_.
    const double other_offset_dist_; // In direction normal to other_surf_.
    const double epsgeo_;

    /// Compute the parameter point in t.  If pcv_turned the returned point takes this into
    /// account.
    /// \param space_cv the curve defining the parametrization.
    /// \param t the parameter in which to evaluate.
    /// \param_cv the curve in which to evaluate.
    /// \param pcv_turned whether the parametrization of param_cv is the opposite of space_cv.
    /// return the corresponding parameter point in param_cv.
    std::vector<double>
    getSuggestedSurfaceParameter(const SplineCurve& space_cv, double t,
				 const SplineCurve& param_cv,
				 bool pcv_turned) const;

};

} // namespace Go

#endif //_PROJECTCURVE_
