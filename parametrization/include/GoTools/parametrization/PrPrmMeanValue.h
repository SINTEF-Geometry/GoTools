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

#ifndef PRPRMMEANVALUE_H
#define PRPRMMEANVALUE_H

#include "GoTools/parametrization/PrParametrizeInt.h"

/*<PrPrmMeanValue-syntax: */

/** PrPrmMeanValue -  Implement the mean value parametrization
 * by implementing the virtual function makeWeights.
 * See preprint "Mean Value Coordinates", M.S. Floater.
 * This method is now THE recommended method.
 * The mapped points of the mean value parameterization
 * depend infinitely smoothly on the original data points,
 * unlike the shape-preserving parameterization.
 * However, in many examples that have been tested,
 * the two methods give similar results.
 * (MF. August 2002).
 */
class PrPrmMeanValue : public PrParametrizeInt
{
protected:

  virtual bool makeWeights(int i);

public:
  /// Default constructor
  PrPrmMeanValue();
  /// Empty destructor
  virtual ~PrPrmMeanValue();

};


/*>PrPrmMeanValue-syntax: */

/*Class:PrPrmMeanValue

Name:              PrPrmMeanValue
Syntax:	           @PrPrmMeanValue-syntax
Keywords:
Description:       Implement the mean value parametrization
                   by implementing the virtual function makeWeights.
                   See preprint "Mean Value Coordinates", M.S. Floater.
                   This method is now THE recommended method.
                   The mapped points of the mean value parameterization
                   depend infinitely smoothly on the original data points,
                   unlike the shape-preserving parameterization.
                   However, in many examples that have been tested,
                   the two methods give similar results.
                   (MF. August 2002).
Member functions:
                   "attach(PrOrganizedPoints& graph)" --\\
                   Set the planar graph.

Constructors:
Files:
Example:
See also:
Developed by:      SINTEF Applied Mathematics, Oslo, Norway
Author:	           Michael Floater, SINTEF
Date:              Apr. 2002
*/

#endif // PRPRMMEANVALUE_H
