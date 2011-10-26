//===========================================================================
//
// File : test_derivVolume.C
//
// Created: Wed Nov 26 13:01:31 2008
//
// Author: Kjell Fredrik Pettersen
//
// Revision: $Id: test_derivVolume.C,v 1.2 2009-01-02 11:37:22 vsk Exp $
//
// Description:
//
//===========================================================================



#include <fstream>
#include "GoTools/trivariate/SplineVolume.h"
#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/utils/errormacros.h"

using namespace Go;
using namespace std;



int main(int argc, char* argv[] )
{

    ALWAYS_ERROR_IF(argc != 8 && argc != 10, "Usage: " << argv[0]
		    << " volumeinfile paru parv parw deru derv derw (swap pardir1, swap pardir2)" << endl);

    // Open input volume file
    ifstream is(argv[1]);
    ALWAYS_ERROR_IF(is.bad(), "Bad or no input filename");

    // Read volume from file
    SplineVolume vol;
    ObjectHeader head;
    is >> head >> vol;

    double pu = atof(argv[2]);
    double pv = atof(argv[3]);
    double pw = atof(argv[4]);
    int deru = atoi(argv[5]);
    int derv = atoi(argv[6]);
    int derw = atoi(argv[7]);

    int dtot = deru + derv + derw;
    int dvwtot = derv + derw;

    int dir1=-1, dir2=-1;
    if (argc > 8)
    {
	dir1 = atoi(argv[8]);
	dir2 = atoi(argv[9]);
    }

    vector<Point> pts;
    pts.resize(((dtot+3)*(dtot+2)*(dtot+1))/6);
    vol.point(pts, pu, pv, pw, dtot);

    cout << "From point evaluation in original surface: " << pts[((dtot+2)*(dtot+1)*dtot)/6 + ((dvwtot+1)*dvwtot)/2 + derw] << endl;

  if (dir1 >= 0 && dir2 >=0 && dir1 != dir2)
      vol.swapParameterDirection(dir1, dir2);

    SplineVolume* dervol = vol.derivVolume(deru, derv, derw);
  if (dir1 >= 0 && dir2 >=0 && dir1 != dir2)
      dervol->swapParameterDirection(dir2, dir1);

    Point pderiv;
    dervol->point(pderiv, pu, pv, pw);

    cout << "From point evaluation in derivative surface: " << pderiv << endl;

    ofstream os("data/derivOut.g2");
    dervol->writeStandardHeader(os);
    dervol->write(os);

}
