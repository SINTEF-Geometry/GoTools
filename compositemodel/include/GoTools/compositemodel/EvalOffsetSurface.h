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

#ifndef _EVALOFFSETSURFACE_H
#define _EVALOFFSETSURFACE_H


#include <memory>
#include "GoTools/compositemodel/ftFaceBase.h"
#include "GoTools/geometry/ParamSurface.h"
#include "GoTools/creators/EvalSurface.h"
#include "GoTools/creators/HermiteGrid2D.h"


namespace Go
{

    class BaseSurface
    {
        // point();
        // normal();
        // asSplineSurface();
        // containingDomain();
        // dimension();
        // closestPoint(); // Not needed at the moment.
    };

    class EvalOffsetSurface : public EvalSurface
    {

    public:

        // Constructor
        // @@sbr201612 To be replaced with a ChartSurface at a later stage.
        EvalOffsetSurface(shared_ptr<ftFaceBase> base_sf,
                          double offset_dist, double epsgeo);

        // Destructor
        virtual ~EvalOffsetSurface();

        // Inherited functions from EvalSurface

        /// Evaluate a point on the surface for a given parameter
        /// \param t the parameter for which to evaluate the surface.
        /// \return the evaluated point
        virtual Point eval( double u, double v) const;

        virtual void eval( double u, double v, int n, Point der[]) const; // n = order of diff

        /// Get the start parameter of the curve.
        /// \return the start parameter of the curve.
        virtual double start_u() const;
        virtual double start_v() const;
  
        /// Get the end parameter of the curve.
        /// \return  the end parameter of the curve.
        virtual double end_u() const;
        virtual double end_v() const;

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
        virtual bool approximationOK(double par_u, double par_v, Point approxpos,
                                     double tol1, double tol2) const;

#if 0
        // Debug
        virtual void write(std::ostream& out) const;
#endif

        void gridSelfIntersections(const HermiteGrid2D& grid,
                                   std::vector<int>& grid_self_intersections,
                                   std::vector<double>& radius_of_curv) const;

        void gridKinks(const HermiteGrid2D& grid,
                       const std::vector<shared_ptr<SplineCurve> >& kink_cvs_2d,
                       std::vector<int>& grid_kinks) const;

        // Project the kink curves onto the approximaitng spline_sf_.
        // We also return the parameter curves and surfaces in the original surfaces.
        std::vector<shared_ptr<SplineCurve> >
        getProjKinkCurves(std::vector<pair<shared_ptr<ParamCurve>, shared_ptr<ParamCurve> > >& par_cvs,
                          std::vector<pair<shared_ptr<ParamSurface>, shared_ptr<ParamSurface> > >& sfs);

    private:

        shared_ptr<ftFaceBase> base_sf_;
        shared_ptr<SplineSurface> spline_sf_; // The guide surface defining the domain and
                                                    // iso-lines of the offset surface. An approximation
                                                    // of the underlying surface, not to be used for
                                                    // evaluations.
        double offset_dist_;
        double epsgeo_;

        //@@sbr201612 Should we support reparametrization inside this class?

        // Given a parameter point in the guide surface, find the corresponding surface in the surface
        // set as well as the corresponding parameter values.
        ParamSurface* findLocalSurface(double u, double v,
                                       double& local_u, double& local_v) const;

        // 
        std::vector<shared_ptr<ParamCurve> >
        get3DKinkCurves(std::vector<pair<shared_ptr<ParamCurve>, shared_ptr<ParamCurve> > >& par_cvs,
                        std::vector<pair<shared_ptr<ParamSurface>, shared_ptr<ParamSurface> > >& sfs);
        
    };    // Class EvalOffsetSurface

}

#endif // _EVALOFFSETSURFACE_H

