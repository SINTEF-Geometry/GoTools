//===========================================================================
//                                                                           
// File: RectDomain.h                                               
//                                                                           
// Created: Fri Sep  1 14:56:17 2000                                         
//                                                                           
// Author: Atgeirr F Rasmussen <atgeirr@sintef.no>
//                                                                           
// Revision: $Id: RectDomain.h,v 1.17 2009-05-13 07:30:49 vsk Exp $
//                                                                           
// Description:
//                                                                           
//===========================================================================

#ifndef _RECTDOMAIN_H
#define _RECTDOMAIN_H

#include "GoTools/utils/Array.h"
#include "GoTools/geometry/Domain.h"
#include "GoTools/utils/config.h"

namespace Go
{

  /** Represents a rectangular parameter domain
   *
   */

class GO_API RectDomain : public Domain
{
public:
    /// Constructs an uninitialized domain
    RectDomain() {}

    /// Constructs and defines a rectangular domain by two of its
    /// opposing corners.
    /// \param corner1 the first corner
    /// \param corner2 the opposing corner to 'corner1'
    RectDomain(const Array<double, 2>& corner1, 
	       const Array<double, 2>& corner2);

    /// Virtual destructor, enables safe inheritance.
    virtual ~RectDomain();
    
    // check whether a given parameter pair is located inside the domain.
    // DOXYGEN documentation can be found in the base class header Domain.h
    virtual bool isInDomain(const Array<double, 2>& point, 
			    double tolerance) const;

    // check whether a gien parameter pair is located on the Domain boundary. 
    // DOXYGEN documentation can be found in the base class header Domain.h
    virtual bool isOnBoundary(const Array<double, 2>& point, 
			      double tolerance) const;

    bool isOnCorner(const Array<double, 2>& point, 
		    double tolerance) const;

    int whichBoundary(const Array<double, 2>& point1, const Array<double, 2>& point2, 
		      double tolerance) const;

    /// Find the (u, v) point in the Domain that is closest (using Euclidean distance
    /// in R^2) to a given (u, v) point.  If the given point is in the domain, then 
    /// the answer is obviously the same point.
    // DOXYGEN documentation can be found in the base class header Domain.h
    virtual void closestInDomain(const Array<double, 2>& point,
				 Array<double, 2>& clo_pt,
				 double tolerance) const;

    /// Find the (u, v) point on the boundary of the Domain that is closest
    /// (using Euclidean distance in R^2) to a given (u, v) point.  If the
    /// point is already considered \em on the boundary, then the answer is obviously
    /// the same point.
    // DOXYGEN documentation can be found in the base class header Domain.h
    virtual void closestOnBoundary(const Array<double, 2>& point,
				   Array<double, 2>& clo_bd_pt,
				   double tolerance) const;

    /// Expand the RectDomain just enough to cover the RectDomain given as argument.
    /// \param rd the RectDomain that we want to be covered by 'this' RectDomain.
    void addUnionWith(const RectDomain& rd);

    /// Set 'this' RectDomain to be the intersection of its current extent and that of
    /// the argument RectDomain.
    /// \param rd the RectDomain that we want to intersect with 'this' one.
    void intersectWith(const RectDomain& rd);

    /// Get the RectDomain's smallest value for the first parameter
    /// \return the RectDomain's smallest value for the first parameter
    double umin() const { return ll_[0]; }

    /// Get the RectDomain's largest value for the first parameter
    /// \return the RectDomain's largest value for the first parameter
    double umax() const { return ur_[0]; }

    /// Get the RectDomain's smallest value for the second parameter
    /// \return the RectDomain's smallest value for the second parameter
    double vmin() const { return ll_[1]; }
    
    /// Get the RectDomain's largest value for the second parameter
    /// \return the RectDomain's largest value for the second parameter
    double vmax() const { return ur_[1]; }

    /// Length of diagonal
    double diagLength()
    {
      return ll_.dist(ur_);
    }

    /// Get the 'lower left' corner of this RectDomain.
    /// \return a 2D array containing the 'lower left' corner of this RectDomain
    Array<double, 2> lowerLeft()  const { return ll_; }

    /// Get the 'upper right' corner of this RectDomain
    /// \return a 2D array containing the 'upper right' corner of this RectDomain
    Array<double, 2> upperRight() const { return ur_; }

private:
    // We store the lower left and upper right points
    Array<double, 2> ll_;
    Array<double, 2> ur_;
};


} // namespace Go

#endif // _RECTDOMAIN_H

