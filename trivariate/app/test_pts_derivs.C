//===========================================================================
//
// File : test_pts_derivs.C
//
// Created: Mon Nov 17 12:37:37 2008
//
// Author: Kjell Fredrik Pettersen
//
// Revision: $Id: test_pts_derivs.C,v 1.2 2008-11-18 13:37:14 kfp Exp $
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
    ALWAYS_ERROR_IF(argc < 4, "Usage: " << argv[0]
		    << " inputvol inputpoints deriveDepth" << endl);

    // Open input surface file
    ifstream is(argv[1]);
    ALWAYS_ERROR_IF(is.bad(), "Bad or no input filename");

    ObjectHeader head;
    is >> head;

    // Read surface from file
    SplineVolume vol;
    is >> vol;

    // Get points
    ifstream pts(argv[2]);
    ALWAYS_ERROR_IF(pts.bad(), "Bad or no input filename");
    int n;
    pts >> n;
    vector<double> pt(n*3);
    for (int i = 0; i < n; ++i) {
	pts >> pt[3*i] >> pt[3*i+1] >> pt[3*i+2];
    }

    int derivs = atoi(argv[3]);

    std::vector<Point> derivPts;
    derivPts.resize(((derivs+3)*(derivs+2)*(derivs+1))/6);

    for (int j = 0; j < n; ++j)
      {
	vol.point(derivPts, pt[3*j],  pt[3*j+1],  pt[3*j+2], derivs);
	for (size_t i = 0; i < derivPts.size(); ++i)
	  cout << derivPts[i][0] << " " << derivPts[i][1] << " " << derivPts[i][2] << " ";
	cout << endl;
	cout << endl;
      }

}





