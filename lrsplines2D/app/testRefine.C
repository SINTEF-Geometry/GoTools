//===========================================================================
//                                                                           
// File: createLR
//                                                                           
// Created: 20.10.12
//                                                                           
// Author: Vibeke Skytt
//                                                                           
// Revision: 
//                                                                           
// Description: Create LRSplineSurface from SplineSurface
//                                                                           
//===========================================================================


#include "GoTools/lrsplines2D/LRSplineSurface.h"
#include "GoTools/lrsplines2D/LRBSpline2D.h"
#include "GoTools/lrsplines2D/LRSplinePlotUtils.h"
#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/geometry/SplineCurve.h"
#include "GoTools/geometry/CurveLoop.h"
#include "GoTools/geometry/ObjectHeader.h"

#include <iostream>
#include <fstream>
#include <string.h>

using namespace Go;

int main(int argc, char *argv[])
{
  if (argc != 4) {
    std::cout << "Usage: lrspline_in (.g2) refinement_in lrspline_out.g2 " << std::endl;
    return -1;
  }

  std::ifstream filein(argv[1]);
  std::ifstream filein2(argv[2]);
  std::ofstream fileout(argv[3]);

  ObjectHeader header;
  header.read(filein);
  LRSplineSurface lrsf;
  lrsf.read(filein);
  
  int nmb_refs;
  filein2 >> nmb_refs;
  for (int ki=0; ki<nmb_refs; ++ki)
    {
      double parval, start, end;
      int dir;
      int mult;

      filein2 >> parval;
      filein2 >> start;
      filein2 >> end;
      filein2 >> dir;
      filein2 >> mult;
      lrsf.refine((dir==0) ? XFIXED : YFIXED, parval, start, end, mult);
    }
  puts("Writing lr-spline to file.");
  if (lrsf.dimension() == 1)
    lrsf.to3D();
  lrsf.writeStandardHeader(fileout);
  lrsf.write(fileout);

  return 0;
}
