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

    int dim = vol.dimension();
    vector<double> pts(nmb_u * nmb_v * nmb_w * dim);
    vector<double> du(nmb_u * nmb_v * nmb_w * dim);
    vector<double> dv(nmb_u * nmb_v * nmb_w * dim);
    vector<double> dw(nmb_u * nmb_v * nmb_w * dim);
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
	    //p = Point(pts[ptpos],pts[ptpos+1],pts[ptpos+2]);
	    p = Point(pts.begin()+ptpos, pts.begin()+ptpos+dim);
	    dumpIfFardist("Point",i,j,k,p,pts2[0]);
	    if (derivs)
	      {
		//px = Point(du[ptpos],du[ptpos+1],du[ptpos+2]);
		px = Point(du.begin()+ptpos, du.begin()+ptpos+dim);
		dumpIfFardist("1st derivative",i,j,k,px,pts2[1]);
		//py = Point(dv[ptpos],dv[ptpos+1],dv[ptpos+2]);
		py = Point(dv.begin()+ptpos, dv.begin()+ptpos+dim);
		dumpIfFardist("2nd derivative",i,j,k,py,pts2[2]);
		//pz = Point(dw[ptpos],dw[ptpos+1],dw[ptpos+2]);
		pz = Point(dw.begin()+ptpos, dw.begin()+ptpos+dim);
		dumpIfFardist("3rd derivative",i,j,k,pz,pts2[3]);
	      }
	    ptpos += dim;
	  }

}
