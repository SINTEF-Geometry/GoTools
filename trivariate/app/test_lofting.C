//===========================================================================
//
// File : test_lofting.C
//
// Created: Thu Dec 18 13:42:39 2008
//
// Author: Kjell Fredrik Pettersen
//
// Revision: $Id: test_lofting.C,v 1.1 2008-12-19 09:54:29 kfp Exp $
//
// Description:
//
//===========================================================================



#include <iostream>
#include <fstream>
#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/trivariate/SplineVolume.h"
#include "GoTools/trivariate/LoftVolumeCreator.h"
#include "GoTools/geometry/ObjectHeader.h"

using std::endl;
using std::ifstream;
using std::ofstream;
using std::vector;
using std::shared_ptr;
using namespace Go;


int main(int argc, char* argv[] )
{

  ALWAYS_ERROR_IF(argc != 3, "Usage: " << argv[0]
		  << " surfacesinfile volumeoutfile" << endl);

  // Open input surfaces file
  ifstream is(argv[1]);
  ALWAYS_ERROR_IF(is.bad(), "Bad or no surface input filename");

  // Open output volume file
  ofstream os(argv[2]);
  ALWAYS_ERROR_IF(os.bad(), "Bad output filename");

  vector<shared_ptr<SplineSurface> > surfaces;

  int nmb_sfs;
  is >> nmb_sfs;

  ObjectHeader head;
  for (int i = 0; i < nmb_sfs; ++i)
    {
      shared_ptr<SplineSurface> srf(new SplineSurface());
      is >> head;
      srf->read(is);
      surfaces.push_back(srf);
    }

  SplineVolume* vol = LoftVolumeCreator::loftVolume(surfaces.begin(), (int)surfaces.size());
  vol->writeStandardHeader(os);
  vol->write(os);
}
