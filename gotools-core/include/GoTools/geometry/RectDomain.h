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

    // check whether a given parameter pair is located inside the domain.
    // DOXYGEN documentation can be found in the base class header Domain.h
    virtual int isInDomain2(const Array<double, 2>& point, 
			    double tolerance) const;

    // check whether a gien parameter pair is located on the Domain boundary. 
    // DOXYGEN documentation can be found in the base class header Domain.h
    virtual bool isOnBoundary(const Array<double, 2>& point, 
			      double tolerance) const;

    /// bd = -1 : Not a boundary point
    /// bd = 0: umin
    /// bd = 1: umax
    /// bd = 2: vmin
    /// bd = 3: vmax
    /// bd2 = -1 : Not a corner point
    /// otherwise: as bd
    bool isOnBoundary(const Array<double, 2>& point, 
		      double tolerance, int& bd, int& bd2) const;

    /// Check if a given parameter pair lies on a corner in the domain within
    /// the given tolerance
    bool isOnCorner(const Array<double, 2>& point, 
		    double tolerance) const;

    /// Given two parameter pairs, check if they specify a domain boundary
    /// return value: -1=no boundary, 0=umin, 1=umax, 2=vmin, 3=vmax
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

    /// Check if two domains overlap, boundary overlap within tolerance
    /// included
    bool overlap(const RectDomain& rd, double tol);

    /// Check if two domains overlap, boundary overlap within tolerance
    /// dependent on parameter direction included
    bool overlap(const RectDomain& rd, double tol1, double tol2);

    /// Get the 'lower left' corner of this RectDomain.
    /// \return a 2D array containing the 'lower left' corner of this RectDomain
    Array<double, 2> lowerLeft()  const { return ll_; }

    /// Get the 'upper right' corner of this RectDomain
    /// \return a 2D array containing the 'upper right' corner of this RectDomain
    Array<double, 2> upperRight() const { return ur_; }

    /// Translate box
    void move(Array<double, 2> vec)
    {
      ll_ += vec;
      ur_ += vec;
    }

private:
    // We store the lower left and upper right points
    Array<double, 2> ll_;
    Array<double, 2> ur_;
};


} // namespace Go

#endif // _RECTDOMAIN_H

