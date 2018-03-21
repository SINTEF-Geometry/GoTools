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

#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/geometry/CurveOnSurface.h"
#include "GoTools/geometry/BoundedUtils.h"
#include <fstream>
#include <algorithm>


using namespace Go;
using namespace std;


int main(int argc, char** argv)
{
  if (argc != 4)
    {
      std::cout << "Input parameters: infile (2 surfaces, g2), tolerance, outfile" << std:: endl;
      exit(1);
    }

    // Read the surfaces from file
    std::ifstream input(argv[1]);
    if (input.bad()) {
	std::cerr << "File error (no file or corrupt file specified)."
		  << std::endl;
	return 1;
    }
    
    double tol = atof(argv[2]);

    std::ofstream out(argv[3]);

    shared_ptr<ObjectHeader> header(new ObjectHeader());
    shared_ptr<ParamSurface> surf1(new SplineSurface());
    shared_ptr<ParamSurface> surf2(new SplineSurface());

    header->read(input);
    if (header->classType() != Class_SplineSurface)
      {
	std::cout << "Object one is not a spline surface" << std::endl;
	exit(1);
      }
    surf1->read(input);

    header->read(input);
    if (header->classType() != Class_SplineSurface)
      {
	std::cout << "Object two is not a spline surface" << std::endl;
	exit(1);
      }
    surf2->read(input);


    vector<shared_ptr<CurveOnSurface> > int_seg1, int_seg2;
    BoundedUtils::getIntersectionCurve(surf1, surf2, int_seg1, int_seg2, tol);

    for (size_t ki=0; ki<int_seg1.size(); ++ki)
      {
	int_seg1[ki]->spaceCurve()->writeStandardHeader(out);
	int_seg1[ki]->spaceCurve()->write(out);
      }
}
