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

#ifndef BSPLINEUNIUTILS_H
#define BSPLINEUNIUTILS_H

#include <array>

#include "GoTools/lrsplines2D/BSplineUniLR.h"

namespace Go
{

  /// Utility functionality for keeping track of univariate B-splines (BSplineUniLR)
  ///  in LR spline context.
  namespace BSplineUniUtils
  {
    bool identify_bsplineuni(const BSplineUniLR* bspline, 
			     std::vector<std::unique_ptr<BSplineUniLR> >& bspline_vec,
			     int& ix);
    bool identify_bsplineuni(std::vector<int>::const_iterator start,
			     std::vector<int>::const_iterator end,
			     std::vector<std::unique_ptr<BSplineUniLR> >& bspline_vec,
			     int& ix);

    void insert_univariate(std::vector<std::unique_ptr<BSplineUniLR> >& bspline_vec,
			   BSplineUniLR *bspline, int& ix);

    int 
      last_overlapping_bsplineuni(int knot_ix,
				  std::vector<std::unique_ptr<BSplineUniLR> >& bspline_vec);

    /* bool bsplineuni_range(std::vector<std::unique_ptr<BSplineUniLR> >& bspline_vec, */
    /* 			  int start_ix, int end_ix, int& first, int& last); */

  }; // end namespace BSplineUniUtils


}; // end namespace Go

#endif
