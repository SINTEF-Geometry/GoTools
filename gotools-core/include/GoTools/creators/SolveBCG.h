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

#ifndef _SOLVEBCG_H
#define _SOLVEBCG_H


#include "GoTools/creators/SolveCG.h"

namespace Go
{

  /// Using the biconjugate gradient method to solve sparse and
  /// positive definite matrix systems.

class SolveBCG : public SolveCG
{

public:

  /// Constructor.
  /// \param conv_type in the range 1 to 4. Convergence estimate. The selection
  /// of preconditioner is based on this number
  SolveBCG(int conv_type, bool symm);

  /// Destructor.
  virtual ~SolveBCG();

  /// Prepare for preconditioning.
  virtual void precond(double relaxfac);

  /// Solve the equation system.
  virtual int solve(double *ex, double *eb, int nn);


private:

  double omega_;        // Relaxation parameter.
  int conv_type_; // We operate with 4 different types of convergence estimates.
  bool symm_;

  void precondRILU(double relaxfac); // @@sbr Not sure RILU is what we will be using...

  // Helper funtion, return largest element in array (i.e. the sup-norm).
  double max_abs_element(double* x, int nn);

};

} // end namespace Go


#endif // _SOLVEBCG_H

