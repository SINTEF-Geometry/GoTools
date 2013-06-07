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
#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/trivariate/SplineVolume.h"
#include "GoTools/trivariate/LoftVolumeCreator.h"
#include "GoTools/geometry/ObjectHeader.h"

using std::endl;
using std::ifstream;
using std::ofstream;
using std::vector;
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
