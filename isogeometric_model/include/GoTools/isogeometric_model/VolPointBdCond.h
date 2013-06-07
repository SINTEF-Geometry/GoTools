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

#ifndef __VOLPOINTBDCOND_H
#define __VOLPOINTBDCOND_H


#include <vector>
#include "GoTools/utils/Point.h"
#include "GoTools/isogeometric_model/BlockPointBdCond.h"



namespace Go
{

  class VolSolution;

  // This class represents a pointwise boundary condition of Dirichlet type
  class VolPointBdCond : public BlockPointBdCond
  {
  public:
    // Constructor
    VolPointBdCond(int face_nmb, double param[], Point& condition_value);

    // Destructor
    ~VolPointBdCond();

    // Get the enumeration of affected surface coefficients
    virtual void 
      getCoefficientsEnumeration(std::vector<int>& local_enumeration) const;

    // Get the value of the condition
    virtual Point getConditionValue() const;

    // Get the parameter corresponding to the condition value
    double* getParam() const;

    // Get enumeration of affected surface coefficients and the factors in the
    // equation for requiring the current boundary to interpolate the condition
    virtual void 
      getInterpolationFactors(std::vector<std::pair<int,double> >& factors) const;

    // Get face number
    int faceNumber() const;

  private:
    // The boundary surface it corresponds to
    // facenmb_ = 0: the boundary corresponding to the minimum parameter in the first parameter
    //                 direction (u_min)
    // facenmb_ = 1: the boundary corresponding to the maximum parameter in the first parameter
    //                 direction (u_max)
    // facenmb_ = 2: the boundary corresponding to the minimum parameter in the second parameter
    //                 direction (v_min)
    // facenmb_ = 3: the boundary corresponding to the maximum parameter in the second parameter
    //                 direction (v_max)
    // facenmb_ = 4: the boundary corresponding to the minimum parameter in the third parameter
    //                 direction (w_min)
    // facenmb_ = 5: the boundary corresponding to the maximum parameter in the third parameter
    //                 direction (w_max)
    int facenmb_; 
  
    // The value of the boundary condition of the solution in the given point
    Point value_;

    // The parameter value corresponding to this pointwise condition
    double param_[2];

    // Pointer to the block solution to which this boundary condition belongs
    VolSolution* parent_;

  };   // end class VolPointBdCond

} // end namespace Go


#endif    // #ifndef __VOLPOINTBDCOND_H
