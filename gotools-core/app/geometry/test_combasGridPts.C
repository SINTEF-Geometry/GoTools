//===========================================================================
//
// File : test_combasGridPts.C
//
// Created: Fri Jan 16 12:58:14 2009
//
// Author: Kjell Fredrik Pettersen
//
// Revision: $Id: test_combasGridPts.C,v 1.1 2009-01-16 12:18:00 kfp Exp $
//
// Description:
//
//===========================================================================


#include <fstream>
#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/utils/errormacros.h"

using namespace Go;
using namespace std;


int main(int argc, char* argv[] )
{
    ALWAYS_ERROR_IF(argc != 4, "Usage: " << argv[0]
		    << " surfaceinfile numpts_u numpts_v" << endl);

    // Open input surface file
    ifstream is(argv[1]);
    ALWAYS_ERROR_IF(is.bad(), "Bad or no input filename");

    ObjectHeader head;
    is >> head;

    // Read surface from file
    SplineSurface sf;
    is >> sf;

    vector<vector<double> > pars(2);

    for (int i = 0; i < 2; ++i)
      {
	double current_par = i == 0 ? sf.startparam_u() : sf.startparam_v();
	double end_par = i == 0 ? sf.endparam_u() : sf.endparam_v();
	int nmb = atoi(argv[2+i]);
	double step = (end_par - current_par)/(double)(nmb-1);
	for (int j = 0; j < nmb; ++j, current_par += step) pars[i].push_back(step);
      }

    vector<vector<double> > pts_only, pts_with_derivs, derivs_u, derivs_v;

    sf.computeBasisGrid(pars[0], pars[1],
			 pts_only);
    sf.computeBasisGrid(pars[0], pars[1],
			 pts_with_derivs, derivs_u, derivs_v);

    if (pts_only.size() != pts_with_derivs.size())
      {
	cout << "Test 1: Wrong base size" << endl;
	exit(-1);
      }
    for (size_t i = 0; i < pts_only.size(); ++i)
      {
	if (pts_only[i].size() != pts_with_derivs[i].size())
	  {
	    cout << "Test 1: Wrong base size at position " << i << endl;
	    exit(-1);
	  }
	for (size_t j = 0; j < pts_only[i].size(); ++j)
	  if (pts_only[i][j] != pts_with_derivs[i][j])
	    {
	      cout << "Test 1: Different entries at (" << i << "," << j << "): " << pts_only[i][j] << " and " << pts_with_derivs[i][j] << endl;
	      exit(-1);
	    }
      }

    cout << "Test 1 complete" << endl;

    vector<BasisPtsSf> bp_pts_only;
    vector<BasisDerivsSf> bp_pts_with_derivs;
    sf.computeBasisGrid(pars[0], pars[1],
			 bp_pts_only);
    sf.computeBasisGrid(pars[0], pars[1],
			 bp_pts_with_derivs);

    if (bp_pts_only.size() != bp_pts_with_derivs.size())
      {
	cout << "Test 2: Wrong base size" << endl;
	exit(-1);
      }

    for (size_t i = 0; i < bp_pts_only.size(); ++i)
      {
	BasisPtsSf pt1 = bp_pts_only[i];
	BasisDerivsSf pt2 = bp_pts_with_derivs[i];

	for (int j = 0; j < 2; ++j)
	  if (pt1.param[j] != pt2.param[j])
	    {
	      cout << "Test 2: Wrong param value at position " << i << "," << j << endl;
	      exit(-1);
	    }

	for (int j = 0; j < 2; ++j)
	  if (pt1.left_idx[j] != pt2.left_idx[j])
	    {
	      cout << "Test 2: Wrong left_idx value at position " << i << "," << j << endl;
	      exit(-1);
	    }

	if (pt1.basisValues.size() != pt2.basisValues.size())
	  {
	    cout << "Test 2: Wrong base size at position " << i << endl;
	    exit(-1);
	  }

	for (size_t j = 0; j < pt1.basisValues.size(); ++j)
	  if (pt1.basisValues[j] != pt2.basisValues[j])
	    {
	      cout << "Test 2: Different entries at (" << i << "," << j << "): " << pt1.basisValues[j] << " and " << pt2.basisValues[j] << endl;
	      exit(-1);
	    }
      }

    cout << "Test 2 complete" << endl;

}
