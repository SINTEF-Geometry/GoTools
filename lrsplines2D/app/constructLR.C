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

#include "GoTools/lrsplines2D/LRSplineSurface.h"

#include <iostream>
#include <vector>

using namespace Go;
using namespace std;

 //------------------------------------------------------------------------------
 vector<double> make_regular_kvec(int degree, int num_coefs, double parmin, double parmax)
 //------------------------------------------------------------------------------
 {
   assert (parmin < parmax);

   vector<double> result;
   const int num_intervals = num_coefs - degree;
   const double interval_size = (parmax - parmin) / num_intervals;

   for (int i = 0; i != degree+1; ++i) result.push_back(parmin);
   for (int i = 1; i != num_intervals; ++i) result.push_back(parmin + i * interval_size);
   for (int i = 0; i != degree+1; ++i) result.push_back(parmax);

   return result;
 }

LRSplineSurface copy(const LRSplineSurface& cp )
{
  LRSplineSurface result = cp;
  return result;
} 

int main(int argc, char *argv[])
{
  int deg_x = 2;
  int deg_y = 2;
  int patches_x = 2;
  int patches_y = 2;
  
  const vector<double> kvec_x = make_regular_kvec(deg_x, deg_x + patches_x, 0.0, 1.0);
  const vector<double> kvec_y = make_regular_kvec(deg_y, deg_y + patches_y, 0.0, 1.0);
  /*
  for (int ix=0; ix != kvec_x.size(); ++ix) {
    cout << kvec_x[ix] << " ";
  } cout << endl;

  for (int ix=0; ix != kvec_y.size(); ++ix) {
    cout << kvec_y[ix] << " ";
  } cout << endl;
  */
  double knot_tol = 0.001;

  LRSplineSurface result(deg_x, deg_y, deg_x + patches_x, deg_y + patches_y, 1, &kvec_x[0], &kvec_y[0], knot_tol);

  //result.write(cout);

  LRSplineSurface result2 = copy(result);

  cout << "Object starts here" << endl;
  result2.write(cout); 

  return 0;
}
