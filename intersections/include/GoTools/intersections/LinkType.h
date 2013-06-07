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

#ifndef _LINKTYPE_H
#define _LINKTYPE_H


namespace Go {


/// A classification of IntersectionLink objects. The classification
/// is strictly based on how it was made within the algorithm.

enum LinkType {
    /// 0
    LINK_UNDEFINED = 0,
    /// 1
    SIMPLE_CONE,
    /// 2
    SIMPLE_MONOTONE,
    /// 3
    SIMPLE_IMPLICIT,
    /// 4
    SIMPLE_TWO_POINTS,
    /// 5
    COINCIDENCE_CVPT,
    /// 6
    COINCIDENCE_CVCV,
    /// 7
    COINCIDENCE_SFCV,
    /// 8
    COINCIDENCE_SFPT,
    /// 9
    MICRO_CVCV,
    /// 10
    MICRO_SFPT,
    /// 11
    MICRO_SFCV,
    /// 12
    MICRO_SFSF,
    /// 13
    MICRO_PAR1FUNC,
    /// 14
    MICRO_PAR2FUNC,
    /// 15
    LINEAR_CVCV,
    /// 16
    LINEAR_SFCV,
    /// 17
    LINEAR_SFSF,
    /// 18
    MERGED_UNDEFINED,
    /// 19
    MERGED_SIMPLE_CONE,
    /// 20
    MERGED_SIMPLE_MONOTONE,
    /// 21
    MERGED_SIMPLE_CONE_MONOTONE,
    /// 22
    MERGED_COINCIDENCE_SFCV_SFCV,
    /// 23
    DEG_TRIANGLE,
    /// 24
    POST_ITERATE,
    /// 25
    BRANCH_CONNECTION,
    /// 26
    COMPLEX_SFSF,
    /// 27
    SPLIT_LINK,
    /// 28
    REPAIRED_MISSING_LINK,
    /// 29
    INSIDE_OUTSIDE_SINGULARITY_BOX

};


} // namespace Go


#endif // _LINKTYPE_H

