//==========================================================================
//                                                                          
// File: test_CvCvIntersector.C
//
// Created:
//                                                                          
// Author: B. Spjelkavik <bsp@sintef.no>
//                                                                          
// Revision: $Id: test_SfSfsisl.C,v 1.3 2005-09-22 15:01:29 oan Exp $
//                                                                          
// Description: Test of class CvCvIntersector
//                                                                          
//==========================================================================

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
