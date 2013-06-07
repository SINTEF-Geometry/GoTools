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

#ifndef _CONSTRAINTDEFINITIONS_H
#define _CONSTRAINTDEFINITIONS_H

#include <vector>


namespace Go
{

    /// Struct defining linear side constraints between control points in
    /// a surface.
    /// The elements of factor_ correspond to a linear combination of
    /// control points on the left side of the equation, whilst
    /// constant_term_ denotes the right side of the equation.
    typedef struct sideConstraint
    {
      /// Dimension of coefficient (max 3).
      int dim_;  
      /// For each coefficient involved in the constraint, the 
      /// index of the coefficient is given and the factor corresponding
      /// to the coefficient in the equation.
      std::vector<std::pair<int, double> > factor_;
      /// The constant term in the current equation.
      double constant_term_[3];  
    } sideConstraint;


    /// Struct defining linear side constraints between control points in
    /// a set of surfaces.
    /// The elements of factor_ correspond to a linear combination of
    /// control points on the left side of the equation, whilst
    /// constant_term_ denotes the right side of the equation.
    typedef struct sideConstraintSet
    {
      /// Dimension of coefficient (max 3).
      int dim_;  
      /// For each coefficient involved in the constraint, the index
      /// of the surface, the index of the coefficient and the 
      /// factor corresponding to the coefficient in the equation is given.
      std::vector<std::pair<std::pair<int,int>, double> > factor_; 
      /// The constant term in the current equation, on the right side of the
      ///  equation. For dim < 3 not all elements are used.
      double constant_term_[3];
      /// Default constructor.
      sideConstraintSet()
      {
	dim_ = 3;
	constant_term_[0] = constant_term_[1] = constant_term_[2] = 0.0;
      }

      /// Constructor.
      /// \param dim dimension of geometric space.
      sideConstraintSet(int dim)
      {
	dim_ = dim;
	constant_term_[0] = constant_term_[1] = constant_term_[2] = 0.0;
      }
    } sideConstraintSet;

}

#endif // _CONSTRAINTDEFINITIONS_H

