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
#include "GoTools/trivariate/SplineVolume.h"
#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/utils/errormacros.h"
#include "GoTools/trivariate/CoonsPatchVolumeGen.h"

using namespace Go;
using namespace std;

int main(int argc, char* argv[] )
{
  if (argc != 3 && argc != 4)
    {
      cout << "Usage: " << argv[0] << " volumeinfile [instructionsinfile] outfile" << endl;
      exit(-1);
    }


  // Open input volume file
  ifstream is(argv[1]);
  ALWAYS_ERROR_IF(is.bad(), "Bad or no volume input filename");

  // Open output volume file
  ofstream os(argv[argc-1]);
  ALWAYS_ERROR_IF(os.bad(), "Bad volume output filename");

  // Read volume from file
  ObjectHeader head;
  SplineVolume vol;
  is >> head >> vol;

  vector<SplineSurface*> faces, old_faces;
  faces.push_back(vol.constParamSurface(vol.startparam(0), 0));
  faces.push_back(vol.constParamSurface(vol.endparam(0), 0));
  faces.push_back(vol.constParamSurface(vol.startparam(1), 1));
  faces.push_back(vol.constParamSurface(vol.endparam(1), 1));
  faces.push_back(vol.constParamSurface(vol.startparam(2), 2));
  faces.push_back(vol.constParamSurface(vol.endparam(2), 2));

  if (argc == 4)
    {
      // Open instructions input file
      ifstream instrf(argv[2]);
      ALWAYS_ERROR_IF(instrf.bad(), "Bad or no instructions input filename");

      // Read instructions. The following code is the only existing documentation
      int cmd;
      instrf >> cmd;
      while (cmd != 0)
	{
	  if (cmd == 1) // Permute u,v,w-directions
	    {
	      int perm[3];
	      instrf >> perm[0] >> perm[1] >> perm[2];
	      old_faces.resize(6);
	      for (int i = 0; i < 6; ++i) old_faces[i] = faces[i];
	      for (int i = 0; i < 6; ++i) faces[i] = old_faces[(perm[i>>1]<<1) | (i&1)];
	    }

	  else if (cmd == 2) // Swap surfaces in same direction
	    {
	      int pos;
	      instrf >> pos;
	      pos <<= 1;
	      SplineSurface *ss = faces[pos];
	      faces[pos] = faces[pos+1];
	      faces[pos+1] = ss;
	    }

	  else if (cmd == 3) // Reverse parameter direction
	    {
	      int pos, dir;
	      instrf >> pos >> dir;
	      faces[pos]->reverseParameterDirection(dir == 0);
	    }

	  else if (cmd == 4) // Swap parameter direction
	    {
	      int pos;
	      instrf >> pos;
	      faces[pos]->swapParameterDirection();
	    }

	  else if (cmd == 5) // Rescale knot intervals
	    {
	      int pos;
	      double begin_u, end_u, begin_v, end_v;
	      instrf >> pos >> begin_u >> end_u >> begin_v >> end_v;
	      faces[pos]->setParameterDomain(begin_u, end_u, begin_v, end_v);
	    }

	  else if (cmd == 6) // Insert new knot value
	    {
	      int pos, dir;
	      double knot_val;
	      instrf >> pos >> dir >> knot_val;
	      if (dir==0)
		faces[pos]->insertKnot_u(knot_val);
	      else
		faces[pos]->insertKnot_v(knot_val);
	    }

	  else if (cmd == 7) // Raise order
	    {
	      int pos, raise_u, raise_v;
	      instrf >> pos >> raise_u >> raise_v;
	      faces[pos]->raiseOrder(raise_u, raise_v);
	    }

	  instrf >> cmd;
	}
    }


  SplineVolume* volnew = Go::CoonsPatchVolumeGen::createCoonsPatch(faces[0], faces[1], faces[2], faces[3], faces[4], faces[5]);

  /*
  vol.setParameterDomain(0.0, 1.0, 0.0, 1.0, 0.0, 1.0);
  SplineVolume* volnew = Go::CoonsPatchVolumeGen::createCoonsPatchDirectly(vol.constParamSurface(vol.startparam(0), 0),
									   vol.constParamSurface(vol.endparam(0), 0),
									   vol.constParamSurface(vol.startparam(1), 1),
									   vol.constParamSurface(vol.endparam(1), 1),
									   vol.constParamSurface(vol.startparam(2), 2),
									   vol.constParamSurface(vol.endparam(2), 2));
  */

  volnew->writeStandardHeader(os);
  volnew->write(os);
}





