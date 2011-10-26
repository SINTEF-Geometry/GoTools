//===========================================================================
//
// File : test_grid.C
//
// Created: Thu Nov 20 11:04:15 2008
//
// Author: Kjell Fredrik Pettersen
//
// Revision: $Id: test_grid.C,v 1.1 2008-11-21 07:21:23 kfp Exp $
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


void dumpIfFardist(string msg,int i,int j,int k,Point p0,Point p1)
{
  if (p0.dist2(p1) > 1e-18)
    cout << "Error: " << msg << " Gridpos=(" << i << "," << j << "," << k << ")  grid = " << p0 << "  isolated = " << p1 << endl;
}
  

int main(int argc, char* argv[] )
{
    ALWAYS_ERROR_IF(argc != 6, "Usage: " << argv[0]
		    << " volumeinfile numpts_u numpts_v numpts_w derivs?" << endl);

    // Open input surface file
    ifstream is(argv[1]);
    ALWAYS_ERROR_IF(is.bad(), "Bad or no input filename");

    ObjectHeader head;
    is >> head;

    // Read volume from file
    SplineVolume vol;
    is >> vol;

    int nmb_u = atoi(argv[2]);
    int nmb_v = atoi(argv[3]);
    int nmb_w = atoi(argv[4]);

    bool derivs = atoi(argv[5]) > 0;

    vector<double> pts(nmb_u * nmb_v * nmb_w * 3);
    vector<double> du(nmb_u * nmb_v * nmb_w * 3);
    vector<double> dv(nmb_u * nmb_v * nmb_w * 3);
    vector<double> dw(nmb_u * nmb_v * nmb_w * 3);
    vector<double> paru(nmb_u);
    vector<double> parv(nmb_u);
    vector<double> parw(nmb_u);

    if (derivs)
      vol.gridEvaluator(nmb_u, nmb_v, nmb_w, pts, du, dv, dw, paru, parv, parw);
    else
      vol.gridEvaluator(nmb_u, nmb_v, nmb_w, pts, paru, parv, parw);

    vector<Point> pts2(4);

    int ptpos = 0;
    for (int k = 0; k < nmb_w; ++k)
      for (int j = 0; j < nmb_v; ++j)
	for (int i = 0; i < nmb_u; ++i)
	  {
	    vol.point(pts2, paru[i], parv[j], parw[k], 1);
	    Point p, px, py, pz;
	    p = Point(pts[ptpos],pts[ptpos+1],pts[ptpos+2]);
	    dumpIfFardist("Point",i,j,k,p,pts2[0]);
	    if (derivs)
	      {
		px = Point(du[ptpos],du[ptpos+1],du[ptpos+2]);
		dumpIfFardist("1st derivative",i,j,k,px,pts2[1]);
		py = Point(dv[ptpos],dv[ptpos+1],dv[ptpos+2]);
		dumpIfFardist("2nd derivative",i,j,k,py,pts2[2]);
		pz = Point(dw[ptpos],dw[ptpos+1],dw[ptpos+2]);
		dumpIfFardist("3rd derivative",i,j,k,pz,pts2[3]);
	      }
	    ptpos += 3;
	  }

}
