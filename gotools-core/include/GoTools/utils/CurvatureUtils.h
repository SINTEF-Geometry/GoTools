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

#ifndef _CURVATUREUTILS_
#define _CURVATUREUTILS_

#include "GoTools/utils/Point.h"
#include <vector>

namespace Go
{
    /** Help functions in related to curvature.
     */

    /// Given position, first and second derivative
    /// of a curve passing through a point, compute
    /// the unit tangent, curvature vector and curvature 
    /// radius of this curve.
    double curvatureRadius(const std::vector<Point>& der,
			   std::vector<Point>& unitder);

    /// Computes the step length along a curve based on radius of curvature
    /// at a point on the curve, and an absolute tolerance.
    double stepLenFromRadius(double radius, double aepsge);

    /// To create the tangent length for interpolating a
    /// circular arc with an almost equi-oscillating Hermit qubic.
    double tanLenFromRadius(double radius, double angle);

    /// Given position, first and second derivative in both ends of
    /// an Hermite segment, compute parameter interval and tangent lengths
    /// in order to stay close to a circular segment.
    void getHermiteData(const std::vector<Point>& der1,
			const std::vector<Point>& der2, 
			double& parint, double& len1, double& len2);

} // End of namespace Go


#endif

