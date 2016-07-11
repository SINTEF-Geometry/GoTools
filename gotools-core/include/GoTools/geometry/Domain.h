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

#ifndef _DOMAIN_H
#define _DOMAIN_H


#include "GoTools/utils/config.h"
#include "GoTools/utils/Array.h"


namespace Go
{
    /** Abstract base class representing a 2D parameter domain.
     *
     */

class GO_API Domain
{
 public:
    /// virtual destructor ensures safe inheritance
    virtual ~Domain() {}

    /// check whether a given parameter pair is located inside the domain
    /// \param point the (u,v)-pair that we want to test.
    /// \param tolerance the tolerance used (ruling what to do when 'point' is
    ///        located very near the edge of the domain).
    /// \return 'true' if 'point' is inside the domain, or within 'tolerance' 
    ///         from being inside the domain, 'false' otherwise'.
    virtual bool isInDomain(const Array<double, 2>& point, 
			    double tolerance) const = 0;

    /// Query whether a given parameter pair is inside the domain or
    /// not.
    /// \param point array containing the parameter pair
    /// \param tolerance the tolerance to be used.  In order to be considered
    ///                  'inside', the point must be located inside one of the
    ///                  defining CurveLoop s, as well as being at a distance
    ///                  more than 'tolerance' from any point on that CurveLoop.
    /// \return '1' if the point is found to be inside the domain, 
    ///         '2' if it is found to be on the boundary '0' otherwise.
    virtual int isInDomain2(const Array<double, 2>& point, 
			    double tolerance) const = 0;

    /// check whether a given parameter pair is located on the domain boundary
    /// \param point the (u,v)-pair that we want to test
    /// \param tolerance the tolerance used (how 'far' from the boundary our
    ///        (u,v) pair can be and still be considered 'on' the boundary.
    /// \return 'true' if the point is considered to be on the boundary (within
    ///         'tolerance', 'false' otherwise.
    virtual bool isOnBoundary(const Array<double, 2>& point, 
			      double tolerance) const = 0;

    /// Find the (u, v) point in the Domain that is closest (using Euclidean distance
    /// in R^2) to a given (u, v) point.  If the given point is in the domain, then 
    /// the answer is obviously the same point.
    /// \param point the (u,v) parameter pair that we want to find the closest 
    ///              parameter pair to \em inside Domain.
    /// \param clo_pt the resulting closest parameter point.
    /// \param tolerance the tolerance used in defining whether the given point is 
    ///        already inside the domain.
    virtual void closestInDomain(const Array<double, 2>& point,
				 Array<double, 2>& clo_pt,
				 double tolerance) const = 0;
    
    /// Find the (u, v) point on the boundary of the Domain that is closest
    /// (using Euclidean distance in R^2) to a given (u, v) point.  If the
    /// point is already considered \em on the boundary, then the answer is obviously
    /// the same point.
    /// \param point the (u,v) parameter pair that we want to find to closest parameter
    ///              pair to \em on the Domain border.
    /// \param clo_bd_pt the resulting closest border point.
    /// \param tolerance the tolerance used in defining whether the given point is
    ////                 already inside the domain.
    virtual void closestOnBoundary(const Array<double, 2>& point,
				   Array<double, 2>& clo_bd_pt,
				   double tolerance) const = 0;

};


} // namespace Go


#endif // _DOMAIN_H

