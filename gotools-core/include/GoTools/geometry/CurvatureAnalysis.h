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

#ifndef _CURVATUREANALYSIS_H
#define _CURVATUREANALYSIS_H


#include "GoTools/geometry/ParamSurface.h"
#include "GoTools/utils/Point.h"

/// \file CurvatureAnalysis.h
/// Curvature analysis related to surfaces.
/// Functions computing fundamental forms and curvature.

namespace Go
{

/// Curvature analysis related to surfaces.
/// Functions computing fundamental forms and curvature.
namespace CurvatureAnalysis
{

    /// Computes the coefficients of the first fundamental form.
    /// The returned values are stored 
    /// { E F G [Eu Fu Gu Ev Fv Gv
    ///   [Euu Fuu Guu Euv Fuv Guv Evv Fvv Gvv [...] ] ] }.
    /// Only one derivative is implemented so far.
    /// \param sf reference to the concerned surface
    /// \param u first parameter of point where we want to compute the first 
    ///          fundamental form
    /// \param v second parameter of point where we want to compute the first
    ///          fundamental form
    /// \param derivs number of (partial) derivatives of the first fundamental
    ///               form that we want to include in the result.  So far, 
    ///               only the computation of the first derivative is actually 
    ///               implemented.
    /// \param form the computed values, on the format described above
    void computeFirstFundamentalForm(const ParamSurface& sf,
				     double u, double v, int derivs,
				     std::vector<double>& form);

    /// Computes the coefficients of the first and second fundamental
    /// forms.
    /// The returned values are stored 
    /// { E F G } in form1 and { e f g } in form2.
    /// This function cannot compute derivatives.
    /// \param sf reference to the concerned surface
    /// \param u first parameter of point where we want to compute the fundamental
    ///          forms.
    /// \param v second parameter of point where we want to compute the fundamental
    ///          forms
    /// \param form1 the values associated with the first fundamental form will be 
    ///              returned here.
    /// \param form2 the values associated with the second fundamental form will be
    ///              returned here.
    void computeSecondFundamentalForm(const ParamSurface& sf,
				      double u, double v,
				      double form1[3],
				      double form2[3]);

    /// Computes the Gaussian (K) and mean (H) curvatures.
    /// \param sf reference to the concerned surface
    /// \param u first parameter of point where we want to carry out computation
    /// \param v second parameter of point where weh want to carry out computation
    /// \param K value of Gaussian curvature returned here
    /// \param H value of mean curvature returned here
    void curvatures(const ParamSurface& sf,
		    double u, double v,
		    double& K, double& H);

    /// Computes the principal curvatures of a surface, sf, in the parameter
    /// value (u,v) along with the corresponding parameter directions.
    void principalCurvatures(const ParamSurface& sf,
			     double u, double v,
			     double& k1, Point& d1,
			     double& k2, Point& d2);

    /// Estimate the minimum curvature radius of the surface sf. 
    /// \param sf the given surface
    /// \param tolerance influences the density of the search
    /// \param mincurve the computed minimum curvature radius
    /// \param pos_u parameter in the first parameter direction corresponding
    /// to the minimum curvature radius
    /// \param pos_v parameter in the second parameter direction corresponding
    /// to the minimum curvature radius
    /// \param degenerate_eps maximum length of a degenenerate boundary, used
    /// to check whether special treatment is required
    /// \param curv_tol default parameter used to limit the search
    void minimalCurvatureRadius(const ParamSurface& sf,
				double tolerance,
				double& mincurv,
				double& pos_u,
				double& pos_v,
				double degenerate_eps,
				double curv_tol = 1.0e-3);


    /// Estimate the minium curvature radius of the surface sf in the parameter
    /// domain [start_u,end_u]x[start_v,end_v] by evaluating the curvature radius
    /// in a grid and fetch the smallest one. The grid density is computed from
    /// the given tolerance
    void evaluateMinCurvatureRadius(const ParamSurface& sf,
				    double start_u, double end_u, double start_v, double end_v,
				    double tolerance,
				    std::vector<double>& param_u, std::vector<double>& param_v,
				    std::vector<std::vector<double> >& curvs,
				    double& mincurv,
				    double& minpos_u,
				    double& minpos_v,
				    bool initialize);


} //namespace CurvatureAnalysis


} // namespace Go

#endif // _CURVATUREANALYSIS_H

