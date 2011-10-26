//===========================================================================
//
// File : make_cone.C
//
// Created: Tue Dec 16 12:16:06 2008
//
// Author: Kjell Fredrik Pettersen
//
// Revision: $Id:$
//
// Description:
//
//===========================================================================


#include <iostream>
#include <fstream>
#include "GoTools/geometry/SplineCurve.h"
#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/geometry/SweepSurfaceCreator.h"
#include "GoTools/utils/Point.h"
#include "GoTools/geometry/ObjectHeader.h"

using namespace Go;
using namespace std;

int main(int argc, char** argv)
{

  if (argc != 11)
    {
      cout << "Usage: " << argv[0] << " top_x top_y top_z bot_x bot_y bot_z axis_x axis_y axis_z outfile" << endl;
      exit(-1);
    }

  Point top (atoi(argv[1]), atoi(argv[2]), atoi(argv[3]));
  Point bot (atoi(argv[4]), atoi(argv[5]), atoi(argv[6]));
  Point axis (atoi(argv[7]), atoi(argv[8]), atoi(argv[9]));

  SplineCurve* line = new SplineCurve(bot, top);
  SplineSurface* srf = SweepSurfaceCreator::rotationalSweptSurface(*line,
								   2 * M_PI,
								   top,
								   axis);

  ofstream outfile(argv[10]);
  if (outfile.bad())
    {
      cout << "Error when creating output file " << argv[1] << endl;
      exit(-1);
    }

  srf->writeStandardHeader(outfile);
  srf->write(outfile);

}
