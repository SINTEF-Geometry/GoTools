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

#ifndef __VOLUMEINTERPOLATOR_H
#define __VOLUMEINTERPOLATOR_H

#include "GoTools/geometry/BsplineBasis.h"
#include "GoTools/utils/Point.h"

namespace Go
{
  class SplineVolume;

/// This namespace contains functions used to interpolate a set of points
  namespace VolumeInterpolator
  {
    /// Interpolate a set of regular, parameterized interpolation
    /// points. The parameterization is assumed to correspond to the
    /// given B-spline basises (select interpolation points in the
    /// Greville points). The function throws if the input parameters
    /// are inconsistent
    /// \param basis_u spline basis in the first parameter direction
    /// \param basis_v spline basis in the second parameter direction
    /// \param basis_w spline basis in the third parameter direction
    /// \param par_u parameter values in 1. parameter direction corresponding
    /// to point
    /// \param par_v parameter values in 2. parameter direction corresponding
    /// to point
    /// \param par_w parameter values in 2. parameter direction corresponding
    /// to point
    /// \param points the regular point set
    /// \param dimension dimension of geometry space
    /// \param rational whether or not a rational surface is expected
    /// \param weights the weights of the rational volume, used only if 
    /// rational==true
    SplineVolume* regularInterpolation(const BsplineBasis& basis_u,
				       const BsplineBasis& basis_v,
				       const BsplineBasis& basis_w,
				       std::vector<double>& par_u,
				       std::vector<double>& par_v,
				       std::vector<double>& par_w,
				       std::vector<double>& points,
				       int dimension,
				       bool rational,
				       std::vector<double>& weights);

  };    // namespace VolumeInterpolator


} // namespace Go


#endif    // #ifndef __VOLUMEINTPEROLATOR_H

