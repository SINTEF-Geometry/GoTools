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
 
  if (argc != 4) {
    cout << "Usage: test_SfSfIntersector FileSf1 FileSf2 aepsge"
	 << endl;
    return 0;
  }


  ObjectHeader header;

  // Read the first curve from file
  ifstream input1(argv[1]);
  if (input1.bad()) {
    cerr << "File #1 error (no file or corrupt file specified)."
	 << std::endl;
    return 1;
  }
  header.read(input1);
  shared_ptr<SplineSurface> surf1(new SplineSurface());
  surf1->read(input1);
  input1.close();
    
  // Read the second curve from file
  ifstream input2(argv[2]);
  if (input2.bad()) {
    cerr << "File #2 error (no file or corrupt file specified)."
	 << std::endl;
    return 1;
  }
  header.read(input2);
  shared_ptr<SplineSurface> surf2(new SplineSurface());
  surf2->read(input2);
  input2.close();

  double aepsge;
  aepsge = atof(argv[3]);

  SISLSurf *sf1 = GoSurf2SISL(*surf1.get());
  SISLSurf *sf2 = GoSurf2SISL(*surf2.get());

  int kstat = 0;
  int kpnt, kcrv;
  double *spar1=0, *spar2=0;
  SISLIntcurve **ucrv=0;
  s1859(sf1, sf2, 0.0, aepsge, &kpnt, &spar1, &spar2, &kcrv,
	&ucrv, &kstat);

  printf("kstat = %d \n",kstat);
  printf("Number of points: %d \n",kpnt);
  printf("Number of curves: %d \n",kcrv);

  freeSurf(sf1);
  freeSurf(sf2);
  free(spar1);
  free(spar2);
  freeIntcrvlist(ucrv, kcrv);
  return 0;
}
