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

#ifndef _RECTTRIDOMAIN_H
#define _RECTTRIDOMAIN_H

#include "GoTools/utils/Array.h"
#include "GoTools/geometry/Domain.h"
#include "GoTools/utils/config.h"

namespace Go
{

  /** Represents a rectangular parameter domain
   *
   */

class GO_API RectTriDomain 
{
public:
    /// Constructs an uninitialized domain
    RectTriDomain() {}

    /// Constructs and defines a rectangular domain by two of its
    /// opposing corners.
    /// \param corner1 the first corner
    /// \param corner2 the opposing corner to 'corner1'
    RectTriDomain(const Array<double, 3>& corner1, 
	       const Array<double, 3>& corner2);

    /// destructor
    ~RectTriDomain();
    
    // check whether a given parameter pair is located inside the domain.
    // DOXYGEN documentation can be found in the base class header Domain.h
    bool isInDomain(const Array<double, 3>& point, 
			    double tolerance) const;

    // check whether a given parameter pair is located inside the domain.
    // DOXYGEN documentation can be found in the base class header Domain.h
    int isInDomain2(const Array<double, 3>& point, 
			    double tolerance) const;

    // check whether a gien parameter pair is located on the Domain boundary. 
    // DOXYGEN documentation can be found in the base class header Domain.h
    bool isOnBoundary(const Array<double, 3>& point, 
			      double tolerance) const;

    /// Find the (u, v) point in the Domain that is closest (using Euclidean distance
    /// in R^2) to a given (u, v) point.  If the given point is in the domain, then 
    /// the answer is obviously the same point.
    // DOXYGEN documentation can be found in the base class header Domain.h
    void closestInDomain(const Array<double, 3>& point,
				 Array<double, 3>& clo_pt,
				 double tolerance) const;

    /// Find the (u, v) point on the boundary of the Domain that is closest
    /// (using Euclidean distance in R^2) to a given (u, v) point.  If the
    /// point is already considered \em on the boundary, then the answer is obviously
    /// the same point.
    // DOXYGEN documentation can be found in the base class header Domain.h
    void closestOnBoundary(const Array<double, 3>& point,
				   Array<double, 3>& clo_bd_pt,
				   double tolerance) const;

    /// Expand the RectTriDomain just enough to cover the RectTriDomain given as argument.
    /// \param rd the RectTriDomain that we want to be covered by 'this' RectTriDomain.
    void addUnionWith(const RectTriDomain& rd);

    /// Set 'this' RectTriDomain to be the intersection of its current extent and that of
    /// the argument RectTriDomain.
    /// \param rd the RectTriDomain that we want to intersect with 'this' one.
    void intersectWith(const RectTriDomain& rd);

    /// Get the RectTriDomain's smallest value for the first parameter
    /// \return the RectTriDomain's smallest value for the first parameter
    double umin() const { return ll_[0]; }

    /// Get the RectTriDomain's largest value for the first parameter
    /// \return the RectTriDomain's largest value for the first parameter
    double umax() const { return ur_[0]; }

    /// Get the RectTriDomain's smallest value for the second parameter
    /// \return the RectTriDomain's smallest value for the second parameter
    double vmin() const { return ll_[1]; }
    
    /// Get the RectTriDomain's largest value for the second parameter
    /// \return the RectTriDomain's largest value for the second parameter
    double vmax() const { return ur_[1]; }

    /// Get the RectTriDomain's smallest value for the third parameter
    /// \return the RectTriDomain's smallest value for the third parameter
    double wmin() const { return ll_[1]; }
    
    /// Get the RectTriDomain's largest value for the third parameter
    /// \return the RectTriDomain's largest value for the third parameter
    double wmax() const { return ur_[1]; }

     /// Length of diagonal
    double diagLength()
    {
      return ll_.dist(ur_);
    }

    /// Check if two domains overlap, boundary overlap within tolerance
    /// included
    bool overlap(const RectTriDomain& rd, double tol);

    /// Get the 'lower left' corner of this RectTriDomain.
    /// \return a 2D array containing the 'lower left' corner of this RectTriDomain
    Array<double, 3> lowerLeft()  const { return ll_; }

    /// Get the 'upper right' corner of this RectTriDomain
    /// \return a 2D array containing the 'upper right' corner of this RectTriDomain
    Array<double, 3> upperRight() const { return ur_; }

    /// Translate box
    void move(Array<double, 3> vec)
    {
      ll_ += vec;
      ur_ += vec;
    }

private:
    // We store the lower left and upper right points
    Array<double, 3> ll_;
    Array<double, 3> ur_;
};


} // namespace Go

#endif // _RECTTRIDOMAIN_H

