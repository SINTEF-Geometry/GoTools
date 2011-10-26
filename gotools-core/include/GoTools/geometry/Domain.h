//===========================================================================
//                                                                           
// File: Domain.h                                                          
//                                                                           
// Created: Fri Sep  1 14:46:32 2000                                         
//                                                                           
// Author: Atgeirr F Rasmussen <atgeirr@sintef.no>
//                                                                           
// Revision: $Id: Domain.h,v 1.11 2007-12-04 16:12:01 jbt Exp $
//                                                                           
// Description:
//                                                                           
//===========================================================================

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

