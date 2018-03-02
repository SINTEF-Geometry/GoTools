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


#include "GoTools/creators/OffsetCurveNormalDir.h"

#include "GoTools/creators/CoonsPatchGen.h"
#include "GoTools/geometry/LineCloud.h"
#include "GoTools/utils/config.h"
#include "GoTools/geometry/ParamCurve.h"

#include <fstream>

using namespace std;
using namespace Go;


//===========================================================================
OffsetCurveNormalDir::OffsetCurveNormalDir(shared_ptr<Go::ParamCurve>& parameter_crv,
                                           shared_ptr<Go::ParamCurve>& space_crv,
                                           shared_ptr<Go::ParamSurface>& surf,
                                           double epsgeo,
                                           double offset_dist)
    : parameter_crv_(parameter_crv), space_crv_(space_crv), surf_(surf),
      epsgeo_(epsgeo), offset_dist_(offset_dist)
//===========================================================================
{
    // Test input
    ALWAYS_ERROR_IF(((parameter_crv_.get() == nullptr) && (space_crv_.get() == nullptr)) ||
                    surf_.get() == nullptr,
                    "Missing input data.");
    ALWAYS_ERROR_IF(parameter_crv_.get() != nullptr && parameter_crv_->dimension() != 2,
                    "Dimension mismatch.");
    ALWAYS_ERROR_IF(space_crv_.get() != nullptr && space_crv_->dimension() != 3,
                    "Dimension mismatch.");
    ALWAYS_ERROR_IF(surf_->dimension() != 3,
                    "Dimension mismatch.");
}


//===========================================================================
OffsetCurveNormalDir::~OffsetCurveNormalDir()
//===========================================================================
{
}


//===========================================================================
Point OffsetCurveNormalDir::eval(double t) const
//===========================================================================
{
    Point par_pt(2);
    if (parameter_crv_.get() != nullptr)
    {
        par_pt = parameter_crv_->ParamCurve::point(t);
    }
    else
    {
        // We project onto the surface.
        Point clo_pt;
        double clo_u, clo_v, clo_dist;
        Point space_pt = space_crv_->point(t);
        const double epsilon = 1.0e-08;
        surf_->closestPoint(space_pt, clo_u, clo_v, clo_pt, clo_dist, epsilon);
        par_pt[0] = clo_u;
        par_pt[1] = clo_v;
    }

    Point space_pt = surf_->ParamSurface::point(par_pt[0], par_pt[1]);
    Point sf_normal;
    surf_->normal(sf_normal, par_pt[0], par_pt[1]);
    Point offset_pt = space_pt + offset_dist_*sf_normal;

    return offset_pt;
}


//===========================================================================
void OffsetCurveNormalDir::eval(double t, int n, Point der[]) const
//===========================================================================
{
    MESSAGE_IF(n > 1, "Only one derivative will be computed");

    if (n == 0)
    {
        der[0] = eval(t);
    }
    else
    {
        Point par_pt(2);
        if (parameter_crv_.get() != nullptr)
        {
            par_pt = parameter_crv_->ParamCurve::point(t);
            if (n > 0)
            {
                THROW("Derivatives not yet supported!");
            }
        }
        else
        {
            // We project onto the surface.
            Point clo_pt;
            double clo_u, clo_v, clo_dist;
            vector<Point> space_pt = space_crv_->point(t, 1);
            const double epsilon = 1.0e-08;
            surf_->closestPoint(space_pt[0], clo_u, clo_v, clo_pt, clo_dist, epsilon);
            par_pt[0] = clo_u;
            par_pt[1] = clo_v;
            if (n > 0)
            {
                der[1] = space_pt[1];
                // @@sbr201802 We use the tangent of the input curve, which is only correct for planes ...
                // For our current usage it should suffice as the boundary loop will be given as external curves.
                if (surf_->instanceType() != Class_Plane)
                {
                    MESSAGE("Derivatives not yet supported (will only be ok for planes)!");
                }
            }
        }

        Point space_pt = surf_->ParamSurface::point(par_pt[0], par_pt[1]);
        Point sf_normal;
        surf_->normal(sf_normal, par_pt[0], par_pt[1]);
        Point offset_pt = space_pt + offset_dist_*sf_normal;

        der[0] = offset_pt;
    }
}


//===========================================================================
double OffsetCurveNormalDir::start() const
//===========================================================================
{
    if (parameter_crv_.get() != nullptr)
    {
        return parameter_crv_->startparam();
    }
    else
    {
        return space_crv_->startparam();
    }
}


//===========================================================================
double OffsetCurveNormalDir::end() const
//===========================================================================
{
    if (parameter_crv_.get() != nullptr)
    {
        return parameter_crv_->endparam();
    }
    else
    {
        return space_crv_->endparam();
    }
}


//===========================================================================
int OffsetCurveNormalDir::dim() const
//===========================================================================
{
    return 3; // Dimension of lifted curve, not that of the parameter curve.
}


//===========================================================================
bool OffsetCurveNormalDir::approximationOK(double par, Point approxpos,
                                           double tol1, double tol2) const
//===========================================================================
{
    // Only first tolerance is used.
  Point pos = eval(par);

  double dist = pos.dist(approxpos); // Both points are in space.

  // @@sbr Currently no input tolerance is used, only the epsgeo_.
  return (dist < epsgeo_);
}
