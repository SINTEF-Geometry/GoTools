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

#include "GoTools/utils/QRFactorization.h"

using namespace Go;
using std::vector;
#include <iostream>
#include <ostream>

int main(int argc, char *argv[])
{
  vector<double> A(12);
  vector<double> b(4);

  A[0] = A[2] = A[3] = A[6] = A[9] = b[0] = b[1] = b[2] = b[3] = 1.0;
  A[1] = A[4] = 3.0;
  A[7] = A[10] = -1.0;
  A[5] = 7.0;
  A[8] = -4.0;
  A[11] = 2.0;

  int m = 4, n = 3, dim = 1;

  vector<double> Q, R, x;
  QRFactorization::QRDecomp(A, n, m, Q, R);

  std::cout << "Q: " << std::endl;
  for (int kb=0; kb<m; ++kb)
    {
      for (int ka=0; ka<m; ++ka)
	std::cout << Q[kb*m+ka] << " ";
      std::cout << std::endl;
    }

  std::cout << "R: " << std::endl;
  for (int kb=0; kb<m; ++kb)
    {
      for (int ka=0; ka<n; ++ka)
	std::cout << R[kb*m+ka] << " ";
      std::cout << std::endl;
    }

  QRFactorization::QRSolve(Q, R, n, m, b, dim, x);
  std::cout << "x: " << std::endl;
  for (int kb=0; kb<dim; ++kb)
    {
      for (int ka=0; ka<n; ++ka)
	std::cout << x[kb*dim+ka] << " ";
      std::cout << std::endl;
    }
}


