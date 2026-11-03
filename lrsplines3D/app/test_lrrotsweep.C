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

#include <iostream>
#include <fstream>
#include "GoTools/geometry/SplineCurve.h"
#include "GoTools/lrsplines2D/LRSplineSurface.h"
#include "GoTools/lrsplines3D/LRSplineVolume.h"
#include "GoTools/lrsplines3D/SweepVolumeLR.h"
#include "GoTools/geometry/ObjectHeader.h"

using namespace std;
using namespace Go;


int main(int argc, char* argv[] )
{

    ALWAYS_ERROR_IF(argc != 10, "Usage: " << argv[0]
		    << " surfaceinfile angle volumeoutfile point axis" << endl);

    // Open input surface file
    ifstream is_surf(argv[1]);
    ALWAYS_ERROR_IF(is_surf.bad(), "Bad or no surface input filename");

    double angle = atof(argv[2]);

    // Open output volume file
    ofstream os(argv[3]);
    ALWAYS_ERROR_IF(os.bad(), "Bad output filename");

    int dim = 3;
    Point pt(dim), axis(dim);
    for (int ki=0; ki<dim; ++ki)
      pt[ki] = atof(argv[4+ki]);
    for (int ki=0; ki<dim; ++ki)
      axis[ki] = atof(argv[7+ki]);

    // Read surface from file
    LRSplineSurface surf;
    ObjectHeader head;
    is_surf >> head >> surf;

    LRSplineVolume* vol = SweepVolumeLR::rotationalSweptVolume(surf, angle, pt, axis);

    vol->writeStandardHeader(os);
    vol->write(os);

    const Array<double,6> span = vol->parameterSpan();
    std::cout << "Span: " << span[0] << " " << span[1]  << " " << span[2]  << " " << span[3]  << " " << span[4]  << " " << span[5]  << std::endl;



    std::cout << "Num eval: " << std::endl;
    int num;
    std::cin >> num;
    std::ofstream ofpt("pos.g2");
    ofpt << "400 1 0 4 255 0 0 255" << std::endl;
    ofpt << num << std::endl;
    std::cout << "Give " << num << " trivariate parameters" << std::endl;
    for (int ki=0; ki<num; ++ki)
      {
	double u, v, w;
	std::cin >> u >> v >> w;
	Point pp;
	vol->point(pp, u, v, w);
	ofpt << pp << std::endl;
      }
}

