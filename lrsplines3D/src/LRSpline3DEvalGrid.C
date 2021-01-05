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

#include "GoTools/lrsplines3D/LRSpline3DEvalGrid.h"



//==============================================================================
namespace Go
//==============================================================================
{


LRSpline3DEvalGrid::LRSpline3DEvalGrid()
  : dim_(0)
{
}


LRSpline3DEvalGrid::LRSpline3DEvalGrid(LRSplineVolume& lr_spline)
    : dim_(lr_spline.dimension())
{
    assert(!lr_spline.rational());

    orig_dom_ = lr_spline.parameterSpan();

    order_u_ = 1 + lr_spline.degree(XDIR);
    order_v_ = 1 + lr_spline.degree(YDIR);
    order_w_ = 1 + lr_spline.degree(ZDIR);

    // We run through the sf and extract the elements.
    auto iter = lr_spline.elementsBegin();
    while (iter != lr_spline.elementsEnd())
    {
	elements_.push_back(*iter->second);
	++iter;
    }

    mesh_ = lr_spline.mesh(); 

}
/*
void LRSpline3DEvalGrid::testCoefComputation()
{

    std::vector<double> tmpCoefs(dim_*order_u_*order_v_*order_w_);
    std::vector<double> coefs;
    std::vector<double> tmpPoints(dim_*order_u_*order_v_*order_w_);
    double *p =&tmpPoints[0];
    int i=0;
//    for(auto it=orig.elements_begin(); it!=orig.elements_end(); it++, i++) {
    for(auto it=elements_begin(); it!=elements_end(); it++, i++) {
      param_float_type ll_x, ll_y, ll_z, ur_x, ur_y, ur_z;
      low(*it, ll_x, ll_y, ll_z);
      high(*it, ur_x, ur_y, ur_z);
      double ll[3];
      ll[0] = ll_x;
      ll[1] = ll_y;
      ll[2] = ll_z;
      double ur[3];
      ur[0] = ur_x;
      ur[1] = ur_y;
      ur[2] = ur_z;
      int index = i;
      //knotIntervals.push_back(KnotInterval(ll, ur, index));
      evaluateGrid(*it, &tmpPoints[0]);
      //computeBezCoefs(dim_, &tmpPoints[0], order_u_, order_v_, order_w_, &tmpCoefs[0]);
      //coefficients.insert(coefficients.end(), tmpCoefs.begin(), tmpCoefs.end());
      coefs.insert(coefs.end(), tmpCoefs.begin(), tmpCoefs.end());
      // element_boundaries.push_back(ll_x);
      // element_boundaries.push_back(ll_y);
      // element_boundaries.push_back(ur_x);
      // element_boundaries.push_back(ur_y);
    }
    //    m_knotSearcher.reset(new KnotIntervalSearcher(knotIntervals, glm::ivec2(128, 128)));
    // m_patchesU = 10;
    // m_patchesV = 10;
    // init(&element_boundaries[0], &coefficients[0]);
    std::cout << "Done with testCoefComputation()." << std::endl;
}
*/


} // end namespace Go
