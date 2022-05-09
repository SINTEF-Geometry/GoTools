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

#ifndef LR_SPLINE3DMBA_H
#define LR_SPLINE3DMBA_H

#include "GoTools/lrsplines3D/LRSplineVolume.h"
#include "GoTools/lrsplines3D/LRBSpline3D.h"
#include "GoTools/utils/Point.h"

namespace Go
{

  /// Update coefficients of LR B-spline volume using multi level B-spline
  /// approximation (MBA)
  namespace LRSpline3DMBA
  {
    // Update LRSplineVolume according to data points stored in the surface elements
    // using the 3DMBA algorithm
    /// Compute distances between points stored with the volume and
    /// update LRSplineVolume according to data points stored in the volume elements
    /// using the MBA algorithm
    /// \param vol volume to update
    void MBADistAndUpdate(LRSplineVolume *vol);
    
    /// Compute distances between points stored with the volume and
    /// update LRSplineVolume according to data points stored in the volume elements
    /// using the MBA algorithm
    /// Parallel version
    /// \param vol volume to update
    /// \param eps tolerance, not used
    /// \param delta equal to 0.0, not used
    void MBADistAndUpdate_omp(LRSplineVolume *vol, double eps, double delta);
    
    /// Update LRSplineVolume according to data points stored in the volume elements
    /// using the MBA algorithm
    /// \param vol volume to update
    void MBAUpdate(LRSplineVolume *vol);
    //void MBAUpdate_omp(LRSplineVolume *vol);
    
    /// Update LRSplineVolume according to data points stored in specified elements
    /// using the MBA algorithm
    /// \param srf surface to update
   /// \param elems elements where the surface are to be updated
   /// \param elems2 elements influenced by shared B-spline support with elements in elems
    void MBAUpdate(LRSplineVolume *vol, std::vector<Element3D*>& elems,
		   std::vector<Element3D*>& elems2);

    /// Help function to 3DMBAUpdate
    void
      add_contribution(int dim,
    		       std::map<const LRBSpline3D*, Array<double,2> >& target,
    		       const LRBSpline3D* bspline, double nom[], double denom);
    void 
      add_contribution2(int dim,
			std::map<const LRBSpline3D*, Array<double,4> >& target, 
			const LRBSpline3D* bspline, double nom[], double denom);

  }; // end namespace LRSpline3DMBA

}; // end namespace Go

#endif
 
