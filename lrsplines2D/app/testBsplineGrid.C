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

#include "GoTools/lrsplines2D/LRSplineSurface.h"
#include "GoTools/lrsplines2D/DefineRefs2D.h"
#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/geometry/ObjectHeader.h"
#include <fstream>

using namespace Go;
using std::vector;


int main(int argc, char *argv[])
{
  if (argc != 4)
    {
      std::cout << "Usage: input spline surface, number of samples, number of derivatives" << std::endl;
      exit(-1);
    }

  std::ifstream input(argv[1]);
  int nmbpts = atoi(argv[2]);
  int ider = atoi(argv[3]);
  
  // Read spline surface from file
   // Read header specifying the type of geometry entity
  ObjectHeader header;
  try {
    header.read(input);
  }
  catch (...)
    {
      std::cerr << "File empty or corrupt. Exiting" << std::endl;
      exit(-1);
    }
  
  shared_ptr<SplineSurface> surf1(new SplineSurface());
  surf1->read(input);
  shared_ptr<LRSplineSurface> surf2(new LRSplineSurface(surf1.get(),
							1.0e-6));

  std::ofstream of("lrsurf.g2");
  surf2->writeStandardHeader(of);
  surf2->write(of);

  int uorder = surf1->order_u();
  int vorder = surf1->order_v();
  int nmbu = surf1->numCoefs_u();
  double umin = surf2->startparam_u();
  double umax = surf2->endparam_u();
  double vmin = surf2->startparam_v();
  double vmax = surf2->endparam_v();
  double udel = (umax - umin)/(double)(nmbpts-1);
  double vdel = (vmax - vmin)/(double)(nmbpts-1);
  double upar = umin, vpar = vmin;

  vector<double> uvals(nmbpts), vvals(nmbpts);
  for (int ka=0; ka<nmbpts; ++ka, upar+=udel, vpar+=vdel)
    {
      uvals[ka] = upar;
      vvals[ka] = vpar;
    }

  if (ider == 0)
    {
      vector<vector<double> > bvals1, bvals2;
      LRSplineSurface::Bspline2DIndexMap bmap;
      vector<BasisPtsSf> tpresultpt;
      vector<LRBasisPtsSf> lrresultpt;
      surf1->computeBasisGrid(uvals, vvals, bvals1);
      surf1->computeBasisGrid(uvals, vvals, tpresultpt);
      surf2->computeBasisGrid(uvals, vvals, bvals2, bmap);
      surf2->computeBasisGrid(uvals, vvals, lrresultpt);

      std::ofstream of1("tpgrid1.txt");
      std::ofstream of2("tpgrid2.txt");
      std::ofstream of3("lrgrid1.txt");
      std::ofstream of4("lrgrid2.txt");

      for (size_t ki=0; ki<bvals1.size(); ++ki)
	{
	  for (size_t kj=0; kj<bvals1[ki].size(); ++kj)
	    of1 << bvals1[ki][kj] << " ";
	  of1 << std::endl;
	}

      for (size_t ki=0; ki<bvals2.size(); ++ki)
	{
	  for (size_t kj=0; kj<bvals2[ki].size(); ++kj)
	    of3 << bvals2[ki][kj] << " ";
	  of3 << std::endl;
	}

      for (size_t ki=0; ki<tpresultpt.size(); ++ki)
	{
	  int ix = (tpresultpt[ki].left_idx[1]-vorder+1)*nmbu +
	    tpresultpt[ki].left_idx[0]-uorder+1;
	  int ix2 = 0;
	  for (size_t kj=0; kj<tpresultpt[ki].basisValues.size(); ++kj, ++ix2)
	    {
	      if (ix2 == uorder)
		{
		  ix += nmbu;
		  ix2 = 0;
		}
	      of2 <<  tpresultpt[ki].basisValues[kj] << " ";
	    }
	  of2 << std::endl;
	}

      for (size_t ki=0; ki<lrresultpt.size(); ++ki)
	{
	  for (size_t kj=0; kj<lrresultpt[ki].basisValues.size(); ++kj)
	    of4 << lrresultpt[ki].basisValues[kj] << " ";
	  of4 << std::endl;
	}
    }
  else if (ider == 1)
    {
      vector<vector<double> > bvals1, bder1_1, bder1_2, bvals2, bder2_1, bder2_2;
      LRSplineSurface::Bspline2DIndexMap bmap2;
      vector<BasisDerivsSf> tpresultder1;
      vector<LRBasisDerivsSf> lrresultder1;
      surf1->computeBasisGrid(uvals, vvals, bvals1, bder1_1, bder1_2);
      surf1->computeBasisGrid(uvals, vvals, tpresultder1);
      surf2->computeBasisGrid(uvals, vvals, bvals2, bder2_1, bder2_2, bmap2);
      surf2->computeBasisGrid(uvals, vvals, lrresultder1);

      std::ofstream of5("tpgrid1_1.txt");
      std::ofstream of6("tpgrid2_1.txt");
      std::ofstream of7("lrgrid1_1.txt");
      std::ofstream of8("lrgrid2_1.txt");

      for (size_t ki=0; ki<bvals1.size(); ++ki)
	{
	  for (size_t kj=0; kj<bvals1[ki].size(); ++kj)
	    of5 << bvals1[ki][kj] << " " << bder1_1[ki][kj] << " " << bder1_2[ki][kj] << " ";
	  of5 << std::endl;
	}

      for (size_t ki=0; ki<bvals2.size(); ++ki)
	{
	  for (size_t kj=0; kj<bvals2[ki].size(); ++kj)
	    of7 << bvals2[ki][kj] << " " << bder2_1[ki][kj] << " " << bder2_2[ki][kj] << " ";
	  of7 << std::endl;
	}

     for (size_t ki=0; ki<tpresultder1.size(); ++ki)
	{
	  int ix = (tpresultder1[ki].left_idx[1]-vorder+1)*nmbu +
	    tpresultder1[ki].left_idx[0]-uorder+1;
	  int ix2 = 0;
	  for (size_t kj=0; kj<tpresultder1[ki].basisValues.size(); ++kj, ++ix2)
	    {
	      if (ix2 == uorder)
		{
		  ix += nmbu;
		  ix2 = 0;
		}
	      of6 <<  tpresultder1[ki].basisValues[kj] << " " << tpresultder1[ki].basisDerivs_u[kj] << " " << tpresultder1[ki].basisDerivs_v[kj] << " ";
	    }
	  of6 << std::endl;
	}

      for (size_t ki=0; ki<lrresultder1.size(); ++ki)
	{
	  for (size_t kj=0; kj<lrresultder1[ki].basisValues.size(); ++kj)
	    of8 << lrresultder1[ki].basisValues[kj] << " "<< lrresultder1[ki].basisDerivs_u[kj] << " " << lrresultder1[ki].basisDerivs_v[kj] << " ";
	  of8 << std::endl;
	}
    }

  else if (ider == 2)
    {
      LRSplineSurface::Bspline2DIndexMap bmap2;
      vector<BasisDerivsSf2> tpresultder2;
      vector<LRBasisDerivsSf2> lrresultder2;
      surf1->computeBasisGrid(uvals, vvals, tpresultder2);
      surf2->computeBasisGrid(uvals, vvals, lrresultder2);

      std::ofstream of9("tpgrid2_2.txt");
      std::ofstream of10("lrgrid2_2.txt");
	
      for (size_t ki=0; ki<tpresultder2.size(); ++ki)
	{
	  int ix = (tpresultder2[ki].left_idx[1]-vorder+1)*nmbu +
	    tpresultder2[ki].left_idx[0]-uorder+1;
	  int ix2 = 0;
	  for (size_t kj=0; kj<tpresultder2[ki].basisValues.size(); ++kj, ++ix2)
	    {
	      if (ix2 == uorder)
		{
		  ix += nmbu;
		  ix2 = 0;
		}
	      of9 << tpresultder2[ki].basisValues[kj] << " ";
	      of9 << tpresultder2[ki].basisDerivs_u[kj] << " " << tpresultder2[ki].basisDerivs_v[kj] << " ";
	      of9 << tpresultder2[ki].basisDerivs_uu[kj] << " " << tpresultder2[ki].basisDerivs_uv[kj] << " "<< tpresultder2[ki].basisDerivs_vv[kj] << std::endl;
	    }
	  of9 << std::endl;
	}

      for (size_t ki=0; ki<lrresultder2.size(); ++ki)
	{
	  for (size_t kj=0; kj<lrresultder2[ki].basisValues.size(); ++kj)
	    {
	      of10 << " " << lrresultder2[ki].basisValues[kj] << " ";
	      of10 << lrresultder2[ki].basisDerivs_u[kj] << " " << lrresultder2[ki].basisDerivs_v[kj] << " ";
	      of10 << lrresultder2[ki].basisDerivs_uu[kj] << " " << lrresultder2[ki].basisDerivs_uv[kj] << " " << lrresultder2[ki].basisDerivs_vv[kj] << std::endl;
	    }
	  of10 << std::endl;
	}
    }
  else if (ider == 3)
    {
      LRSplineSurface::Bspline2DIndexMap bmap2;
      vector<BasisDerivsSf3> tpresultder3;
      vector<LRBasisDerivsSf3> lrresultder3;
      surf1->computeBasisGrid(uvals, vvals, tpresultder3);
      surf2->computeBasisGrid(uvals, vvals, lrresultder3);

      std::ofstream of11("tpgrid2_3.txt");
      std::ofstream of12("lrgrid2_3.txt");
      std::ofstream of13("der3.txt");
	
      for (size_t ki=0; ki<tpresultder3.size(); ++ki)
	{
	  int ix = (tpresultder3[ki].left_idx[1]-vorder+1)*nmbu +
	    tpresultder3[ki].left_idx[0]-uorder+1;
	  int ix2 = 0;
	  for (size_t kj=0; kj<tpresultder3[ki].basisValues.size(); ++kj, ++ix2)
	    {
	      if (ix2 == uorder)
		{
		  ix += nmbu;
		  ix2 = 0;
		}
	      of11 << tpresultder3[ki].basisValues[kj] << " ";
	      of11 << tpresultder3[ki].basisDerivs_u[kj] << " " << tpresultder3[ki].basisDerivs_v[kj] << " ";
	      of11 << tpresultder3[ki].basisDerivs_uu[kj] << " " << tpresultder3[ki].basisDerivs_uv[kj] << " "<< tpresultder3[ki].basisDerivs_vv[kj] << " ";
	      of11 << tpresultder3[ki].basisDerivs_uuu[kj] << " " << tpresultder3[ki].basisDerivs_uuv[kj] << " "<< tpresultder3[ki].basisDerivs_uvv[kj] << " "<< tpresultder3[ki].basisDerivs_vvv[kj] << std::endl;
	    }
	  of11 << std::endl;
	}

      for (size_t ki=0; ki<lrresultder3.size(); ++ki)
	{
	  vector<Point> deriv3(4);
	  for (int ka=0; ka<4; ++ka)
	    deriv3[ka].setValue(0.0, 0.0, 0.0);
	  
	  for (size_t kj=0; kj<lrresultder3[ki].basisValues.size(); ++kj)
	    {
	      of12 << lrresultder3[ki].basisValues[kj] << " ";
	      of12 << lrresultder3[ki].basisDerivs_u[kj] << " " << lrresultder3[ki].basisDerivs_v[kj] << " ";
	      of12 << lrresultder3[ki].basisDerivs_uu[kj] << " " << lrresultder3[ki].basisDerivs_uv[kj] << " " << lrresultder3[ki].basisDerivs_vv[kj] << " ";
	      of12 << lrresultder3[ki].basisDerivs_uuu[kj] << " " << lrresultder3[ki].basisDerivs_uuv[kj] << " " << lrresultder3[ki].basisDerivs_uvv[kj] << " " << lrresultder3[ki].basisDerivs_vvv[kj]<< std::endl;

	      Point coef = lrresultder3[ki].bsplines[kj]->Coef();
	      deriv3[0] += coef*lrresultder3[ki].basisDerivs_uuu[kj];
	      deriv3[1] += coef*lrresultder3[ki].basisDerivs_uuv[kj];
	      deriv3[2] += coef*lrresultder3[ki].basisDerivs_uvv[kj];
	      deriv3[3] += coef*lrresultder3[ki].basisDerivs_vvv[kj];
	    }

	  vector<Point> pts(10);
	  double u = lrresultder3[ki].param[0];
	  double v = lrresultder3[ki].param[1];
	  surf2->point(pts, u, v, 3);

	  of13 << u << " " << v << std::endl;
	  of13 << deriv3[0] << ", " << pts[6] << " " << deriv3[0].dist(pts[6]) << std::endl;
	  of13 << deriv3[1] <<  ", " << pts[7] << " " << deriv3[1].dist(pts[7]) <<  std::endl;
	  of13 << deriv3[2] <<  ", " << pts[8] << " " << deriv3[2].dist(pts[8]) <<  std::endl;
	  of13 << deriv3[3] <<  ", " << pts[9] << " " << deriv3[3].dist(pts[9]) <<  std::endl << std::endl;

	  of12 << std::endl;
	}
    }
}
    



  
 
