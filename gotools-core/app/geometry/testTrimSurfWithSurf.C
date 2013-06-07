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

#include "GoTools/geometry/BoundedSurface.h"
#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/geometry/BoundedUtils.h"

#include <fstream>
#include <vector>


using namespace Go;
using std::vector;


int main(int argc, char** argv)
{
    if (argc != 5) {
	MESSAGE("Usage: surf1 surf2 epsge outfile");
	return 1;
    }

    std::ifstream infile1(argv[1]);
    std::ifstream infile2(argv[2]);
    double epsge = atof(argv[3]);
    std::ofstream outfile(argv[4]);

    ObjectHeader header;
    shared_ptr<SplineSurface> surf1(new SplineSurface());
    shared_ptr<SplineSurface> surf2(new SplineSurface());
    header.read(infile1);
    surf1->read(infile1);
    header.read(infile2);
    surf2->read(infile2);

    vector<shared_ptr<BoundedSurface> > trimmed_sfs =
	BoundedUtils::trimSurfWithSurf(surf1, surf2, epsge);

    for (size_t ki = 0; ki < trimmed_sfs.size(); ++ki) {
	trimmed_sfs[ki]->writeStandardHeader(outfile);
	trimmed_sfs[ki]->write(outfile);
    }

    return 0;
}
