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

#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/lrsplines2D/LRSplineSurface.h"
#include "GoTools/geometry/GoTools.h"

#include <iostream>
#include <fstream>
#include <vector>

using namespace Go;
using std::vector;

int main(int argc, char *argv[])
{
  if (argc != 4)
    {
      std::cout << "Usage: input surface, number of samples, number of derivatives" << std::endl;
      exit(-1);
    }

  std::ifstream input(argv[1]);
  int nmbpts = atoi(argv[2]);
  int ider = atoi(argv[3]);
  
  // Read volume from file
   // Read header specifying the type of geometry entity
  GoTools::init();
 ObjectHeader header;
  try {
    header.read(input);
  }
  catch (...)
    {
      std::cerr << "File empty or corrupt. Exiting" << std::endl;
      exit(-1);
    }
  
  shared_ptr<LRSplineSurface> surf(new LRSplineSurface());
  surf->read(input);

  double umin = surf->startparam_u();
  double umax = surf->endparam_u();
  double vmin = surf->startparam_v();
  double vmax = surf->endparam_v();
  double udel = (umax - umin)/(double)(nmbpts-1);
  double vdel = (vmax - vmin)/(double)(nmbpts-1);
  double upar = umin, vpar = vmin;

  vector<double> uvals(nmbpts), vvals(nmbpts);
  for (int ka=0; ka<nmbpts; ++ka, upar+=udel, vpar+=vdel)
    {
      uvals[ka] = upar;
      vvals[ka] = vpar;
    }

  std::ofstream of("der.txt");
  if (ider == 0)
    {
      vector<vector<double> > bvals;
      LRSplineSurface::Bspline2DIndexMap bmap;
      vector<LRBasisPtsSf> lrresultpt;
      surf->computeBasisGrid(uvals, vvals, lrresultpt);
      surf->computeBasisGrid(uvals, vvals, bvals, bmap);
       
      for (size_t ki=0; ki<lrresultpt.size(); ++ki)
	{
	  Point pnt1;
	  double u = lrresultpt[ki].param[0];
	  double v = lrresultpt[ki].param[1];
	  surf->point(pnt1, u, v);

	  double diff = 0.0;
	  Point pnt2(0.0, 0.0, 0.0);
	  for (size_t kj=0; kj<lrresultpt[ki].basisValues.size(); ++kj)
	    {
	      Point coef = lrresultpt[ki].bsplines[kj]->Coef();
	      pnt2 += coef*lrresultpt[ki].basisValues[kj];
	      size_t ix = bmap.at(lrresultpt[ki].bsplines[kj]);
	      diff += fabs(lrresultpt[ki].basisValues[kj] - bvals[ki][ix]);
	    }


	  of << u << " " << v << std::endl;
	  of << pnt1 << ", " << pnt2 << ", " << pnt1.dist(pnt2) << std::endl;
	  of << "bdiff: " << diff << std::endl << std::endl;
	}
    }
  
  else if (ider == 1)
    {
      vector<vector<double> > bvals, bder1, bder2;
      LRSplineSurface::Bspline2DIndexMap bmap;
      vector<LRBasisDerivsSf> lrresultder;
      surf->computeBasisGrid(uvals, vvals, lrresultder);
      surf->computeBasisGrid(uvals, vvals, bvals, bder1, bder2, bmap);
      
      for (size_t ki=0; ki<lrresultder.size(); ++ki)
	{
	  vector<Point> der1(4);
	  double u = lrresultder[ki].param[0];
	  double v = lrresultder[ki].param[1];
	  surf->point(der1, u, v, 1);

	  double diff = 0.0, diff1 = 0.0, diff2 = 0.0;
	  vector<Point> der2(3);
	  for (int ka=0; ka<3; ++ka)
	    der2[ka].setValue(0.0, 0.0, 0.0);
	  for (size_t kj=0; kj<lrresultder[ki].basisValues.size(); ++kj)
	    {
	      Point coef = lrresultder[ki].bsplines[kj]->Coef();
	      der2[0] += coef*lrresultder[ki].basisValues[kj];
	      der2[1] += coef*lrresultder[ki].basisDerivs_u[kj];
	      der2[2] += coef*lrresultder[ki].basisDerivs_v[kj];
	      size_t ix = bmap.at(lrresultder[ki].bsplines[kj]);
	      diff += fabs(lrresultder[ki].basisValues[kj] - bvals[ki][ix]);
	      diff1 += fabs(lrresultder[ki].basisDerivs_u[kj] - bder1[ki][ix]);
	      diff2 += fabs(lrresultder[ki].basisDerivs_v[kj] - bder2[ki][ix]);
	    }

	  of << u << " " << v << std::endl;
	  of << der1[0] << ", " << der2[0] << ", " << der1[0].dist(der2[0]) << std::endl;
	  of << der1[1] << ", " << der2[1] << ", " << der1[1].dist(der2[1]) << std::endl;
	  of << der1[2] << ", " << der2[2] << ", " << der1[2].dist(der2[2]) << std::endl;
	  of << "bdiff: " << diff << " " << diff1 << " " << diff2 << std::endl << std::endl;
	}
    }
  
  else if (ider == 2)
    {
      vector<LRBasisDerivsSf2> lrresultder;
      surf->computeBasisGrid(uvals, vvals, lrresultder);
      
      for (size_t ki=0; ki<lrresultder.size(); ++ki)
	{
	  vector<Point> der1(6);
	  double u = lrresultder[ki].param[0];
	  double v = lrresultder[ki].param[1];
	  surf->point(der1, u, v, 2);

	  vector<Point> der2(6);
	  for (int ka=0; ka<6; ++ka)
	    der2[ka].setValue(0.0, 0.0, 0.0);
	  for (size_t kj=0; kj<lrresultder[ki].basisValues.size(); ++kj)
	    {
	      Point coef = lrresultder[ki].bsplines[kj]->Coef();
	      der2[0] += coef*lrresultder[ki].basisValues[kj];
	      der2[1] += coef*lrresultder[ki].basisDerivs_u[kj];
	      der2[2] += coef*lrresultder[ki].basisDerivs_v[kj];
	      der2[3] += coef*lrresultder[ki].basisDerivs_uu[kj];
	      der2[4] += coef*lrresultder[ki].basisDerivs_uv[kj];
	      der2[5] += coef*lrresultder[ki].basisDerivs_vv[kj];
	    }

	  of << u << " " << v << std::endl;
	  of << der1[0] << ", " << der2[0] << ", " << der1[0].dist(der2[0]) << std::endl;
	  of << der1[1] << ", " << der2[1] << ", " << der1[1].dist(der2[1]) << std::endl;
	  of << der1[2] << ", " << der2[2] << ", " << der1[2].dist(der2[2]) << std::endl;
	  of << der1[3] << ", " << der2[3] << ", " << der1[3].dist(der2[3]) << std::endl;
	  of << der1[4] << ", " << der2[4] << ", " << der1[4].dist(der2[4]) << std::endl;
	  of << der1[5] << ", " << der2[5] << ", " << der1[5].dist(der2[5]) << std::endl;
	  of << std::endl;
	}
    }
  
  else if (ider == 3)
    {      
      vector<LRBasisDerivsSf3> lrresultder;
      surf->computeBasisGrid(uvals, vvals, lrresultder);
      
      for (size_t ki=0; ki<lrresultder.size(); ++ki)
	{
	  vector<Point> der1(10);
	  double u = lrresultder[ki].param[0];
	  double v = lrresultder[ki].param[1];
	  surf->point(der1, u, v, 3);

	  vector<Point> der2(10);
	  for (int ka=0; ka<10; ++ka)
	    der2[ka].setValue(0.0, 0.0, 0.0);
	  for (size_t kj=0; kj<lrresultder[ki].basisValues.size(); ++kj)
	    {
	      Point coef = lrresultder[ki].bsplines[kj]->Coef();
	      der2[0] += coef*lrresultder[ki].basisValues[kj];
	      der2[1] += coef*lrresultder[ki].basisDerivs_u[kj];
	      der2[2] += coef*lrresultder[ki].basisDerivs_v[kj];
	      der2[3] += coef*lrresultder[ki].basisDerivs_uu[kj];
	      der2[4] += coef*lrresultder[ki].basisDerivs_uv[kj];
	      der2[5] += coef*lrresultder[ki].basisDerivs_vv[kj];
	      der2[6] += coef*lrresultder[ki].basisDerivs_uuu[kj];
	      der2[7] += coef*lrresultder[ki].basisDerivs_uuv[kj];
	      der2[8] += coef*lrresultder[ki].basisDerivs_uvv[kj];
	      der2[9] += coef*lrresultder[ki].basisDerivs_vvv[kj];
	    }

	  of << u << " " << v << std::endl;
	  of << der1[0] << ", " << der2[0] << ", " << der1[0].dist(der2[0]) << std::endl;
	  of << der1[1] << ", " << der2[1] << ", " << der1[1].dist(der2[1]) << std::endl;
	  of << der1[2] << ", " << der2[2] << ", " << der1[2].dist(der2[2]) << std::endl;
	  of << der1[3] << ", " << der2[3] << ", " << der1[3].dist(der2[3]) << std::endl;
	  of << der1[4] << ", " << der2[4] << ", " << der1[4].dist(der2[4]) << std::endl;
	  of << der1[5] << ", " << der2[5] << ", " << der1[5].dist(der2[5]) << std::endl;
	  of << der1[6] << ", " << der2[6] << ", " << der1[6].dist(der2[6]) << std::endl;
	  of << der1[7] << ", " << der2[7] << ", " << der1[7].dist(der2[7]) << std::endl;
	  of << der1[8] << ", " << der2[8] << ", " << der1[8].dist(der2[8]) << std::endl;
	  of << der1[9] << ", " << der2[9] << ", " << der1[9].dist(der2[9]) << std::endl;
	  of << std::endl;
	}
    }
}
