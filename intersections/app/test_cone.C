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

#include <fstream>
#include <iomanip>
#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/intersections/ParamSurfaceInt.h"
#include "GoTools/intersections/SplineSurfaceInt.h"
#include "GoTools/intersections/SfSfIntersector.h"
#include "GoTools/geometry/SISLconversion.h"
#include "sislP.h"


using namespace std;
using namespace Go;


int main(int argc, char** argv)
{
 


  ObjectHeader header;

  // Read the surface from file
  ifstream input1(argv[1]);
  if (input1.bad()) {
    return 1;
  }
  
  header.read(input1);
  shared_ptr<SplineSurface> surf(new SplineSurface());
  surf->read(input1);
  input1.close();

  DirectionCone cone = surf->normalCone();
  printf("Greater than pi : %d \n",cone.greaterThanPi());
  printf("Angle : %13.7f \n",cone.angle());

  int kstat = 0;
  SISLSurf *sf = GoSurf2SISL(*surf.get());
  s1990(sf, 0.000001, &kstat);
  printf("Sisl igtpi : %d, angle %13.7f \n",sf->pdir->igtpi,sf->pdir->aang);
  
  freeSurf(sf);
}
