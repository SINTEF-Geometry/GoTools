//===========================================================================
//
// File : testConstParam.C
//
// Created: Tue Nov 11 14:13:52 2008
//
// Author: Kjell Fredrik Pettersen
//
// Revision: $Id: testConstParam.C,v 1.1 2008-11-12 16:10:08 kfp Exp $
//
// Description:
//
//===========================================================================


#include <fstream>
#include "GoTools/trivariate/SplineVolume.h"
#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/utils/errormacros.h"
#include "GoTools/geometry/BsplineBasis.h"

using namespace Go;
using namespace std;

int main(int argc, char* argv[] )
{
  if (argc != 5 && argc != 6)
    {
      cout << "Usage: " << argv[0] << "infile outfile direction(0,1, or 2) #surfaces(>=2)" << endl;
      cout << "   Or: " << argv[0] << "infile outfile direction(0,1, or 2) 1 parameter" << endl;
      exit(-1);
    }
 
  int nmbsurfs = atoi(argv[4]);
  if (argc != 5 + (nmbsurfs == 1))
    {
      cout << "Usage: " << argv[0] << "infile outfile direction(0,1, or 2) #surfaces(>=2)" << endl;
      cout << "   Or: " << argv[0] << "infile outfile direction(0,1, or 2) 1 parameter" << endl;
      exit(-1);
    }

  // Open input surface file
  ifstream is(argv[1]);
  ALWAYS_ERROR_IF(is.bad(), "Bad or no input filename");

  ObjectHeader head;
  is >> head;

  // Read volume from file
  SplineVolume vol;
  is >> vol;

  // Open outfile
  ofstream os(argv[2]);
  ALWAYS_ERROR_IF(os.bad(), "Bad or no output filename");

  int pardir = atoi(argv[3]);

  double param;
  double step;
  if (nmbsurfs == 1)
    {
      param = atof(argv[5]);
      step = 0.0;
    }
  else
    {
      const BsplineBasis& b = vol.basis(pardir);
      param = b.startparam();
      step = (b.endparam() - param) / (double)(nmbsurfs - 1);
    }

  for (int i = 0; i < nmbsurfs; ++i)
    {
      SplineSurface* ss;
      ss = vol.constParamSurface(param, pardir);
      ss->writeStandardHeader(os);
      ss->write(os);
      param += step;
    }

}





