//===========================================================================
//
// File : test_clone_swap.C
//
// Created: Fri Nov 21 12:07:04 2008
//
// Author: Kjell Fredrik Pettersen
//
// Revision: $Id: test_insert_raise.C,v 1.1 2008-11-24 10:36:11 kfp Exp $
//
// Description:
//
//===========================================================================


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
