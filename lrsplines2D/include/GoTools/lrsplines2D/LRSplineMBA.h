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

#ifndef LR_SPLINEMBA_H
#define LR_SPLINEMBA_H

#include "GoTools/lrsplines2D/LRSplineSurface.h"
#include "GoTools/lrsplines2D/LRBSpline2D.h"
#include "GoTools/utils/Point.h"

namespace Go
{

  namespace LRSplineMBA
  {
    // Update LRSplineSurface according to data points stored in the surface elements
    // using the MBA algorithm
    void MBADistAndUpdate(LRSplineSurface *srf);
    void MBADistAndUpdate_omp(LRSplineSurface *srf);
    void MBAUpdate(LRSplineSurface *srf);
    void MBAUpdate_omp(LRSplineSurface *srf);
    void MBAUpdate(LRSplineSurface *srf, std::vector<Element2D*>& elems,
		   std::vector<Element2D*>& elems2);

    // Help function to MBAUpdate
    void 
      add_contribution(int dim,
		       std::map<const LRBSpline2D*, Array<double,2> >& target, 
		       const LRBSpline2D* bspline, double nom[], double denom);

  }; // end namespace LRSplineMBA

}; // end namespace Go

#endif
 
