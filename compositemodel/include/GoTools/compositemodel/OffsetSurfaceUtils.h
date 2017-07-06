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

#ifndef _OFFSETSURFACEUTILS_H
#define _OFFSETSURFACEUTILS_H

#include "GoTools/geometry/SplineSurface.h"

#include <vector>

enum OffsetSurfaceStatus {
    OFFSET_OK = 0,
    OFFSET_FAILED = 1, // For a lack of more precise status.
    SELF_INTERSECTING_INTERIOR = 2, // Failure to handle. Should have been removed by routine.
    SELF_INTERSECTING_BOUNDARY = 3, // Currently not supported.
    TOLERANCE_ERROR = 4,
    NOT_FOUR_CORNERS = 5,
    NON_ISO_KINK_CURVE = 6
};

namespace Go
{
    

namespace OffsetSurfaceUtils
{

    /// Compute the offset surface for a set of parametric surfaces (possibly 1 surface only).
    /// The surface set must have 4 corners.
    /// \param param_sfs input parametric surfaces
    /// \param offset_dist distance for the offset surface.
    /// \param epsgeo geometric tolerance for the offset surface.
    /// \return status of the offset function. 0 => success, everything else a failure.
    OffsetSurfaceStatus offsetSurfaceSet(const std::vector<shared_ptr<ParamSurface> >& param_sfs,
                                         double offset_dist, double epsgeo,
                                         shared_ptr<SplineSurface>& offset_sf);

} // namespace OffsetSurfaceUtils

} // namespace Go

#endif // _OFFSETSURFACEUTILS_H

