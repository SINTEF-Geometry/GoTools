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
#include  "GoTools/utils/config.h"


namespace Go
{

/// This class represents the curve obtained by projecting a 
/// given 3D curve onto a given part of a given 3D surface.

class ProjectCurve : public EvalCurve
{
public:

    /// Constructor, taking one 3D curve and one 3D surface.  The user 
    /// may additionally provide explicit values for the start and end points
    /// of the curve (supposedly on the surface) as well as specifying
    /// the parameter domain of interest of the surface.
    /// \param space_crv the 3D space curve that will be projected onto a surface
    /// \param surf the 3D surface that we will project the space curve onto.
    ///             Assuming that the surface is k-regular.
    /// \param start_par_pt explicit position of start point (can be a zero pointer,
    ///                     in which case the start point will be evaluated by projection,
    ///                     just like any other point).
    /// \param end_par_pt explicit position of end point (can be a zero pointer,
    ///                   in which case the start point will be evaluated by projection,
    ///                   just like any other point).
    /// \param epsgeo1 geometric tolerance to use when projecting curve onto surface, and 
    ///               when using the approximationOK() function, max dist from exact proj.
    /// \param epsgeo2 geometric tolerance when using the approximationOK() function, max dist from space_crv. Negative = ignored.
    /// \param domain_of_interest if the user wants to limit the surface to a certain 
    ///                           parametric domain, it can be specified here.
    ProjectCurve(shared_ptr<Go::ParamCurve>& space_crv,
		 shared_ptr<Go::ParamSurface>& surf,
		 shared_ptr<Go::Point>& start_par_pt,
		 shared_ptr<Go::Point>& end_par_pt,
		 double epsgeo1,
// 		 double epsgeo2,
		 const RectDomain* domain_of_interest = NULL);

    /// virtual destructor ensures safe inheritance
    virtual ~ProjectCurve();
    
    // Inherited from EvalCurve
    virtual Go::Point eval( double t) const;

    /// Evaluate point, given seed for the closest point iteration involved
    Go::Point eval( double t, Go::Point seed) const;

    // Inherited from EvalCurve
    virtual void eval(double t, int n, Go::Point der[]) const;

    // Inherited from EvalCurve
    virtual double start() const;

    // Inherited from EvalCurve
    virtual double end() const;

    /// Inherited from EvalCurve::dim().  For this class, the returned dimension will be that
    /// of the surface parameter domain, ie. 2, NOT that of the space curve.
    virtual int dim() const;

    /// Inherited from EvalCurve::approximationOK().  For this class, the specified tolerances
    /// are not used; the internally stored 'epsgeo' value is used as tolerance (this value was
    /// specified in the constructor).
    /// \param par the parameter at which to check the curve
    /// \param approxpos the position we want to check whether or not the curve
    ///                  approximates for parameter 'par'.
    /// \param tol1 unused
    /// \param tol2 unused
    /// \return 'true' if the curve approximates the point at the parameter, 'false'
    ///         otherwise.
    virtual bool approximationOK(double par, Go::Point approxpos,
				 double tol1, double tol2) const; 

private:
    const shared_ptr<Go::ParamCurve> space_crv_;
    const shared_ptr<Go::ParamSurface> surf_;
    const shared_ptr<Go::Point> start_par_pt_; // When projecting end pts may be of special interest.
    const shared_ptr<Go::Point> end_par_pt_;
    // We store our found proj pts, useful for seed creation.
//     vector<pair<double, Point> > found_pts_; // (cv_par, sf_par)
    const double epsgeo1_; // Max dist from exact proj.
//     const double epsgeo2_; // Max dist from space_crv. Negative = ignored.
    const RectDomain* domain_of_interest_;
    bool closed_dir_u_; // Closed cfs and projection requires extra caution.
    bool closed_dir_v_;
    double umin_, umax_, vmin_, vmax_;  // Parameter domain of surface

    // Simple function, interpolates end pts.
    std::vector<double> createSeed(double tpar) const;

    void surfaceClosed(const SplineSurface& sf, bool& dir_u, bool& dir_v) const;

    // Assuming surface is closed, make sure we landed on the right
    // side.
    void placeBorderPoint(double t,
			  double& clo_u, double& clo_v) const;

    // If close to the boundary we should check distance when
    // projecting towards boundary.
    // True if inside domain_fraction*(tmax-tmin) in one or both dirs.
    bool closeToSurfaceBoundary(double upar, double vpar,
				double domain_fraction = 0.01) const;

    // If a boundary point is closer to space_pt we snap to that param.
    void snapIfBoundaryIsCloser(Go::Point space_pt,
				double& upar, double& vpar, double& dist) const;

};


}

#endif //_PROJECTCURVE_
