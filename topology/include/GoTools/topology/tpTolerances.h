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

#ifndef _TPTOLERANCES_H_
#define  _TPTOLERANCES_H_

namespace Go
{

/** tpTolerances -  Tolerances used in adjacency analysis for faces
 */
struct tpTolerances
{
public:
  /// Tolerance for when two faces are assumed to be C0 continous
    double gap;
  /// Tolerance for when two faces are assumed to be neighbours
    double neighbour;
  /// Tolerance for when two adjacent faces are assumed to be C1 continous
    double kink;
  /// Tolerance for when two adjacent faces are assumed to have an 
  /// intentially smooth connection
    double bend;
    tpTolerances(double g, double n, double k, double b)
	: gap(g), neighbour(n), kink(k), bend(b)
    {}
    tpTolerances(const tpTolerances& tol)
	: gap(tol.gap), neighbour(tol.neighbour), kink(tol.kink),
          bend(tol.bend)
    {}
};

} // namespace Go

#endif //  _TPTOLERANCES_H_
