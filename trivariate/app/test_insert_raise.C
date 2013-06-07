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

#include "GoTools/trivariate/SplineVolume.h"
#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/utils/errormacros.h"
#include <fstream>
#include <sstream>


using namespace Go;
using namespace std;


int main(int argc, char* argv[] )
{
    ALWAYS_ERROR_IF(argc < 3, "Usage: " << argv[0]
		    << " volumeinfile volumeoutfile [ -i pardir knot1 ... kbotN ] [ -r raise_u raise_v raise_w ]" << endl);

    // Open input volume file
    ifstream is(argv[1]);
    ALWAYS_ERROR_IF(is.bad(), "Bad or no input filename");

    // Open output volume swap file 1
    ofstream os(argv[2]);
    ALWAYS_ERROR_IF(os.bad(), "Bad output filename");

    // Read volume from file
    SplineVolume vol;
    ObjectHeader head;
    is >> head >> vol;

    vector<int> knotIntPardir;
    vector<vector<double> >  knotInt;
    int raise_u = 0;
    int raise_v = 0;
    int raise_w = 0;

    int argPos = 3;
    while (argPos < argc)
      {
	string arg = argv[argPos++];
	if (arg != "-r" && arg != "-i")
	  continue;

	if (arg == "-i")
	  {
	    knotIntPardir.push_back(atoi(argv[argPos++]));
	    vector<double> kInt;
	    while (argPos < argc && string(argv[argPos]) != "-i" && string(argv[argPos]) != "-r")
	      kInt.push_back(atof(argv[argPos++]));
	    knotInt.push_back(kInt);
	  }
	else if (arg == "-r")
	  {
	    raise_u = atoi(argv[argPos++]);
	    raise_v = atoi(argv[argPos++]);
	    raise_w = atoi(argv[argPos++]);
	  }
      }

    for (size_t i = 0; i < knotInt.size(); ++i)
      if (knotInt[i].size() > 0)
	{
	  if (knotInt[i].size() == 1)
	    {
	      cout << "Single knot insrtion. Dir = " << knotIntPardir[i] << " Knot = " << knotInt[i][0] << endl;
	      vol.insertKnot(knotIntPardir[i], knotInt[i][0]);
	    }
	  else
	    {
	      cout << "Multiple knots insrtion. Dir = " << knotIntPardir[i] << " Knots =";
	      for (size_t j = 0; j < knotInt[i].size(); ++j)
		cout << " " << knotInt[i][j];
	      cout << endl;
	      vol.insertKnot(knotIntPardir[i], knotInt[i]);
	    }
	}

    cout << "Raise order " << raise_u << "," << raise_v << "," << raise_w << endl;
    vol.raiseOrder(raise_u, raise_v, raise_w);

    vol.writeStandardHeader(os);
    vol.write(os);
}
