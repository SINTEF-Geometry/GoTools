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

#include "GoTools/geometry/SISLconversion.h"
#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/utils/Array.h"
#include "GoTools/geometry/ClassType.h"
#include "GoTools/geometry/sisl_file_io.h"
#include "sislP.h"

#include <fstream>
#include <stdlib.h>
#include <stdio.h>



int main(int argc, char** argv)
{
  if (argc != 3)
    {
      std::cout << "Usage; infile oufile" << std::endl;
      return -1;
    }

  std::ifstream infile(argv[1]);
  FILE* to_file = fopen(argv[2], "w" );

  int type, dummy, nmb_col_elem;
  std::vector<shared_ptr<Go::GeomObject> > objs;
  while (infile >> type)
    {
      infile >> dummy;
      infile >> dummy;
      infile >> nmb_col_elem;
      for (int ki = 0; ki < nmb_col_elem; ++ki)
	  infile >> dummy;

    if (type == 200)
      {
	shared_ptr<Go::SplineSurface> surf(new Go::SplineSurface());
	surf->read(infile);
	objs.push_back(surf);
      }
    else if (type == 100)
      {
	shared_ptr<Go::SplineCurve> crv(new Go::SplineCurve());
	crv->read(infile);
	objs.push_back(crv);
      }
    else
      {
	std::cout << "Expecting spline curves & surfaces only!" << std::endl;
      }
    }

  // We write to file the number of objects.
  fprintf(to_file,"$ Number of objects in file\n");
  fprintf(to_file,"%d\n\n",int(objs.size()));

  for (size_t ki = 0; ki < objs.size(); ++ki)
    {
      /* open the file with the name to_file */
      if (objs[ki]->instanceType() == Go::Class_SplineSurface)
	{
	  SISLSurf* srf=GoSurf2SISL
	    (*(dynamic_pointer_cast<Go::SplineSurface, Go::GeomObject>(objs[ki])));
	  surface_to_file(to_file, srf);
	}
      else if (objs[ki]->instanceType() == Go::Class_SplineCurve)
	{
	  SISLCurve* cv=Curve2SISL
	    (*(dynamic_pointer_cast<Go::SplineCurve, Go::GeomObject>(objs[ki])));
	  curve_to_file(to_file, cv);
	}
      fprintf(to_file,"\n");
    }

  fclose(to_file);
}
