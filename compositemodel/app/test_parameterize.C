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

#include "GoTools/compositemodel/PointSetApp.h"
#include <fstream>
#include <stdlib.h> // For atof()

using namespace std;
using namespace Go;

int main( int argc, char* argv[] )
{
  if (argc != 3) {
    std::cout << "Input parameters : Input file, Output file"  << std::endl;
    exit(-1);
  }

  // Read input arguments
  std::ifstream is(argv[1]);
  ALWAYS_ERROR_IF(is.bad(), "Input file not found or file corrupt");
  std::ofstream os(argv[2]);

  // Read input file
  int numpnts = 0;
  int numtrs = 0;
  is >> numpnts;
  is >> numtrs;
   
  vector<double> xyz(3*numpnts);
  vector<int> tri(3*numtrs);

  int ki;
  for (ki=0; ki<numpnts; ++ki)
    is >> xyz[3*ki] >> xyz[3*ki+1] >> xyz[3*ki+2];
  for (ki=0; ki<numtrs; ++ki)
    is >> tri[3*ki] >> tri[3*ki+1] >> tri[3*ki+2];

  // Parameterize
  vector<double> param;
  PointSetApp::parameterizeTriang(&xyz[0], numpnts, &tri[0], numtrs,
				  param);

  streamsize prev = os.precision(15);
  os << numpnts << std::endl;
  for (ki=0; ki<numpnts; ++ki)
    {
      os << param[2*ki] << " " << param[2*ki+1] << " ";
      os << xyz[3*ki] << " " << xyz[3*ki+1] << " " << xyz[3*ki+2];
      os << std::endl;
    }
}

