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


#include "GoTools/lrsplines2D/LRFeatureUtils.h"
#include "GoTools/geometry/RectDomain.h"

#include <fstream>

using namespace Go;
using std::vector;
using std::cout;
using std::endl;
using std::ifstream;

//#define DEBUG


//===========================================================================
void  LRFeatureUtils::writeCellInfo(const LRSplineSurface& srf, 
				    double tol, int ncell,
				    std::ostream &out)
//===========================================================================
{
  RectDomain dom = srf.containingDomain();
  double u1 = dom.umin();
  double u2 = dom.umax();
  double v1 = dom.vmin();
  double v2 = dom.vmax();
  double udel = (u2 - u1)/(double)ncell;
  double vdel = (v2 - v1)/(double)ncell;
  int ix = (srf.dimension() == 3) ? 4 : 2;

  int nc2 = ncell*ncell;
  vector<int> nmb_pts(3*nc2, 0);
  vector<double> cellinfo(18*nc2, 0.0);
  int *npt = &nmb_pts[0];
  int *nout_over = npt+nc2;
  int *nout_under = nout_over+nc2;
  double *avdist = &cellinfo[0];
  double *max_over = avdist+nc2;
  double *max_under = max_over+nc2;
  double *minheight = max_under+nc2;
  double *maxheight = minheight+nc2;
  double *nel = maxheight+nc2;
  double *stdd = nel+nc2;
  double *slope = stdd+nc2;
  double *avheight = slope+nc2;
  double *stddheight = avheight+nc2;
  double *av2 = stddheight+nc2;
  double *av_over = av2+nc2;
  double *av_under = av_over+nc2;
  double *av_height_sf = av_under+nc2;
  double *min_height_sf = av_height_sf+nc2;
  double *max_height_sf = min_height_sf+nc2;
  double *laplace = max_height_sf+nc2;

  // Adjust default
  for (int kr=0; kr<nc2; ++kr)
    {
      max_over[kr] = std::numeric_limits<double>::lowest();
      max_under[kr] = std::numeric_limits<double>::max();
      maxheight[kr] = std::numeric_limits<double>::lowest();
      minheight[kr] = std::numeric_limits<double>::max();
      max_height_sf[kr] = std::numeric_limits<double>::lowest();
      min_height_sf[kr] = std::numeric_limits<double>::max();
    }

  // Traverse elements and relate point and element information to
  // cell
  for (LRSplineSurface::ElementMap::const_iterator it=srf.elementsBegin();
       it != srf.elementsEnd(); ++it)
    {
      // Identify overlapping cells
      double uel1 = it->second->umin();
      double uel2 = it->second->umax();
      double vel1 = it->second->vmin();
      double vel2 = it->second->vmax();
      int i1 = (int)((uel1 - u1)/udel);
      int i2 = (int)((uel2 - u1)/udel);
      i1 = std::max(0, std::min(ncell-1, i1));
      i2 = std::max(0, std::min(ncell-1, i2));
      int j1 = (int)((vel1 - v1)/vdel);
      int j2 = (int)((vel2 - v1)/vdel);
      j1 = std::max(0, std::min(ncell-1, j1));
      j2 = std::max(0, std::min(ncell-1, j2));

      // Compute element fraction in cell
      for (int ki=i1; ki<=i2; ++ki)
	{
	  double ufrac = (std::min(uel2, u1+(ki+1)*udel) -
			  std::max(uel1, u1+ki*udel))/(uel2 - uel1);
	  for (int kj=j1; kj<=j2; ++kj)
	    {
	      double vfrac = (std::min(vel2, v1+(kj+1)*vdel) -
			      std::max(vel1, v1+kj*vdel))/(vel2 - vel1);
	      nel[kj*ncell+ki] += (ufrac*vfrac);
	    }
	}

     // Fetch associated data points
      vector<double>& points = it->second->getDataPoints();
      int nmb = it->second->nmbDataPoints();
      int del = it->second->getNmbValPrPoint();
      for (int kr=0; kr<nmb; ++kr)
	{
	  // Identify cell
	  i1 = (int)((points[kr*del] - u1)/udel);
	  i1 = std::max(0, std::min(ncell-1, i1));
	  j1 = (int)((points[kr*del+1] - v1)/vdel);
	  j1 = std::max(0, std::min(ncell-1, j1));
	  int ixcell = j1*ncell + i1;
	  av2[ixcell] += points[kr*del+ix+1];
	  avdist[ixcell] += fabs(points[kr*del+ix+1]);
	  avheight[ixcell] += points[kr*del+ix];
	  max_over[ixcell] = std::max(max_over[ixcell], points[kr*del+ix+1]);
	  max_under[ixcell] = std::min(max_under[ixcell], points[kr*del+ix+1]);
	  maxheight[ixcell] = std::max(maxheight[ixcell], points[kr*del+ix]);
	  minheight[ixcell] = std::min(minheight[ixcell], points[kr*del+ix]);
	  npt[ixcell]++;
	  if (points[kr*del+ix+1] > tol)
	    nout_over[ixcell]++;
	  if (points[kr*del+ix+1] > 0.0)
	    av_over[ixcell] += points[kr*del+ix+1];
	  if (points[kr*del+ix+1] < -tol)
	  if (points[kr*del+ix+1] < -tol)
	    nout_under[ixcell]++;
	  if (points[kr*del+ix+1] < 0.0)
	    av_under[ixcell] -= points[kr*del+ix+1];
	}
    }
  
  for (int kr=0; kr<nc2; ++kr)
    {
      if (npt[kr] > 0)
	{
	  av2[kr] /= (double)npt[kr];
	  avdist[kr] /= (double)npt[kr];
	  avheight[kr] /= (double)npt[kr];
	  av_over[kr] /= (double)npt[kr];
	  av_under[kr] /= (double)npt[kr];
	}
    }

  // Standard deviation
  for (LRSplineSurface::ElementMap::const_iterator it=srf.elementsBegin();
       it != srf.elementsEnd(); ++it)
    {
      // Identify overlapping cells
      double uel1 = it->second->umin();
      double uel2 = it->second->umax();
      double vel1 = it->second->vmin();
      double vel2 = it->second->vmax();

      // Fetch associated data points
      vector<double>& points = it->second->getDataPoints();
      int nmb = it->second->nmbDataPoints();
      int del = it->second->getNmbValPrPoint();
      for (int kr=0; kr<nmb; ++kr)
	{
	  // Identify cell
	  int i1 = (int)((points[kr*del] - u1)/udel);
	  i1 = std::max(0, std::min(ncell-1, i1));
	  int j1 = (int)((points[kr*del+1] - v1)/vdel);
	  j1 = std::max(0, std::min(ncell-1, j1));
	  int ixcell = j1*ncell + i1;
	  double tmp = points[kr*del+ix+1] - av2[ixcell];
	  stdd[ixcell] += tmp*tmp;
	  double tmp2 = points[kr*del+ix] - avheight[ixcell];
	  stddheight[ixcell] += tmp2*tmp2;
	}
    }

  for (int kr=0; kr<nc2; ++kr)
    {
      if (npt[kr] > 0)
	{
	  stdd[kr] /= (double)npt[kr];
	  stddheight[kr] /= (double)npt[kr];
	}
    }

  double upar, vpar;
  int ki, kj;
  for (ki=0, upar=u1; ki<ncell; ++ki, upar+=udel)
    for ( kj=0, vpar=v1; kj<ncell; ++kj, vpar+=vdel)
      {
	// Compute average slope in cell (9 points)
	double slope2 = 0.0;
	double laplace2 = 0.0;
	int ka, kb;
	double upar2, vpar2;
	double udel2 = udel/3.0;
	double vdel2 = vdel/3.0;
	vector<Point> der(6);
	for (ka=0, upar2=upar+0.5*udel2; ka<3; ++ka, upar2+=udel2)
	  for (kb=0, vpar2=vpar+0.5*vdel2; kb<3; ++kb, vpar2+=vdel2)
	    {
	      srf.point(der, upar2, vpar2, 2);
	      slope2 += sqrt(der[1][0]*der[1][0]+der[2][0]*der[2][0]);
	      laplace2 += der[3][0]*der[3][0]+der[5][0]*der[5][0];
	      min_height_sf[kj*ncell+ki] = std::min(min_height_sf[kj*ncell+ki],
						    der[0][0]);
	      max_height_sf[kj*ncell+ki] = std::max(max_height_sf[kj*ncell+ki],
						    der[0][0]);
	      av_height_sf[kj*ncell+ki] += der[0][0];
	    }
	slope[kj*ncell+ki] = slope2/9.0;
	laplace[kj*ncell+ki] = laplace2/9.0;
	av_height_sf[kj*ncell+ki] /= 9.0;
      }

  // Normalize to harmonize the different entries
  double minval = 0.0;
  double maxval = 10.0;
  vector<vector<double> > outval(17);
  for (int kh=0; kh<17; ++kh)
    outval[kh].resize(nc2);
  double currmin = std::numeric_limits<double>::max();
  double currmax = std::numeric_limits<double>::lowest();
  for (int kr=0; kr<nc2; ++kr)
    {
      currmin = std::min(slope[kr], currmin);
      currmax = std::max(slope[kr], currmax);
      outval[0][kr] = slope[kr];
    }
  for (int kr=0; kr<nc2; ++kr)
    outval[0][kr] = (currmax < 1.0e-6) ? 0 : outval[0][kr]*maxval/currmax;

  currmin = std::numeric_limits<double>::max();
  currmax = std::numeric_limits<double>::lowest();
  for (int kr=0; kr<nc2; ++kr)
    {
      currmin = std::min(av_height_sf[kr], currmin);
      currmax = std::max(av_height_sf[kr], currmax);
      outval[1][kr] = av_height_sf[kr];
    }
  for (int kr=0; kr<nc2; ++kr)
    outval[1][kr] = (currmax < 1.0e-6) ? 0 : outval[1][kr]*maxval/currmax;

  currmin = std::numeric_limits<double>::max();
  currmax = std::numeric_limits<double>::lowest();
  for (int kr=0; kr<nc2; ++kr)
    {
      double diff = max_height_sf[kr] - min_height_sf[kr];
      currmin = std::min(diff, currmin);
      currmax = std::max(diff, currmax);
      outval[2][kr] = diff;
    }
  for (int kr=0; kr<nc2; ++kr)
    outval[2][kr] = (currmax < 1.0e-6) ? 0 : outval[2][kr]*maxval/currmax;

  currmin = std::numeric_limits<double>::max();
  currmax = std::numeric_limits<double>::lowest();
  for (int kr=0; kr<nc2; ++kr)
    {
      currmin = std::min(avdist[kr], currmin);
      currmax = std::max(avdist[kr], currmax);
      outval[3][kr] = avdist[kr];
    }
  for (int kr=0; kr<nc2; ++kr)
    outval[3][kr] = (currmax < 1.0e-6) ? 0 : outval[3][kr]*maxval/currmax;

  currmin = std::numeric_limits<double>::max();
  currmax = std::numeric_limits<double>::lowest();
  for (int kr=0; kr<nc2; ++kr)
    {
      double tmp = (npt[kr] == 0) ? 0.0 :
	std::max(fabs(max_over[kr]), fabs(max_under[kr]));
      currmin = std::min(tmp, currmin);
      currmax = std::max(tmp, currmax);
      outval[4][kr] = tmp;
    }
  for (int kr=0; kr<nc2; ++kr)
    outval[4][kr] = (currmax < 1.0e-6) ? 0 : outval[4][kr]*maxval/currmax;


  currmin = std::numeric_limits<double>::max();
  currmax = std::numeric_limits<double>::lowest();
  for (int kr=0; kr<nc2; ++kr)
    {
      currmin = std::min(avheight[kr], currmin);
      currmax = std::max(avheight[kr], currmax);
      outval[5][kr] = avheight[kr];
    }
  for (int kr=0; kr<nc2; ++kr)
    {
      if (currmin < 0)
	outval[5][kr] -= currmin;
      outval[5][kr] =(currmax < 1.0e-6) ? 0 :  outval[5][kr]*maxval/currmax;
    }

  currmin = std::numeric_limits<double>::max();
  currmax = std::numeric_limits<double>::lowest();
  for (int kr=0; kr<nc2; ++kr)
    {
      double tmp = (npt[kr] == 0) ? 0.0 : maxheight[kr] - minheight[kr];
      currmin = std::min(tmp, currmin);
      currmax = std::max(tmp, currmax);
      outval[6][kr] = tmp;
    }
  for (int kr=0; kr<nc2; ++kr)
    outval[6][kr] =(currmax < 1.0e-6) ? 0 :  outval[6][kr]*maxval/currmax;

  currmin = std::numeric_limits<double>::max();
  currmax = std::numeric_limits<double>::lowest();
  for (int kr=0; kr<nc2; ++kr)
    {
      currmin = std::min(stdd[kr], currmin);
      currmax = std::max(stdd[kr], currmax);
      outval[7][kr] = stdd[kr];
    }
  for (int kr=0; kr<nc2; ++kr)
    outval[7][kr] = (currmax < 1.0e-6) ? 0 : outval[7][kr]*maxval/currmax;
 
  currmin = std::numeric_limits<double>::max();
  currmax = std::numeric_limits<double>::lowest();
  for (int kr=0; kr<nc2; ++kr)
    {
      currmin = std::min(stddheight[kr], currmin);
      currmax = std::max(stddheight[kr], currmax);
      outval[8][kr] = stddheight[kr];
    }
  for (int kr=0; kr<nc2; ++kr)
    outval[8][kr] = (currmax < 1.0e-6) ? 0 : outval[8][kr]*maxval/currmax;

  currmin = std::numeric_limits<double>::max();
  currmax = std::numeric_limits<double>::lowest();
  for (int kr=0; kr<nc2; ++kr)
    {
      double div = std::max(fabs(max_under[kr]), fabs(max_over[kr]));
      double tmp = (npt[kr] == 0 || div < 1.0e-6) ? 0.0 : 
	avdist[kr]/div;
      currmin = std::min(tmp, currmin);
      currmax = std::max(tmp, currmax);
      outval[9][kr] = tmp;
    }
  for (int kr=0; kr<nc2; ++kr)
    outval[9][kr] = (currmax < 1.0e-6) ? 0 : outval[9][kr]*maxval/currmax;

  currmin = std::numeric_limits<double>::max();
  currmax = std::numeric_limits<double>::lowest();
  for (int kr=0; kr<nc2; ++kr)
    {
      double tmp = (npt[kr] == 0) ? 0.0 : 
	std::max(max_over[kr], 0.0) - std::min(max_under[kr], 0.0);
      currmin = std::min(tmp, currmin);
      currmax = std::max(tmp, currmax);
      outval[10][kr] = tmp;
    }
  for (int kr=0; kr<nc2; ++kr)
    outval[10][kr] = (currmax < 1.0e-6) ? 0 : outval[10][kr]*maxval/currmax;
  
  currmin = std::numeric_limits<double>::max();
  currmax = std::numeric_limits<double>::lowest();
  for (int kr=0; kr<nc2; ++kr)
    {
      currmin = std::min(av_over[kr], currmin);
      currmax = std::max(av_over[kr], currmax);
      outval[11][kr] = av_over[kr];
    }
  for (int kr=0; kr<nc2; ++kr)
    outval[11][kr] = (currmax < 1.0e-6) ? 0 : outval[11][kr]*maxval/currmax;

  currmin = std::numeric_limits<double>::max();
  currmax = std::numeric_limits<double>::lowest();
  for (int kr=0; kr<nc2; ++kr)
    {
      currmin = std::min(av_under[kr], currmin);
      currmax = std::max(av_under[kr], currmax);
      outval[12][kr] = av_under[kr];
    }
  for (int kr=0; kr<nc2; ++kr)
    outval[12][kr] = (currmax < 1.0e-6) ? 0 : outval[12][kr]*maxval/currmax;
      
   currmin = std::numeric_limits<double>::max();
  currmax = std::numeric_limits<double>::lowest();
  for (int kr=0; kr<nc2; ++kr)
    {
      double tmp = (npt[kr] == 0) ? 0.0 : (double)(nout_under[kr])/(double)npt[kr];
      currmin = std::min(tmp, currmin);
      currmax = std::max(tmp, currmax);
      outval[13][kr] = tmp;
    }
  for (int kr=0; kr<nc2; ++kr)
    outval[13][kr] = (currmax < 1.0e-6) ? 0 : outval[13][kr]*maxval/currmax;
 
  currmin = std::numeric_limits<double>::max();
  currmax = std::numeric_limits<double>::lowest();
  for (int kr=0; kr<nc2; ++kr)
    {
      double tmp = (npt[kr] == 0) ? 0.0 : (double)(nout_over[kr])/(double)npt[kr];
      currmin = std::min(tmp, currmin);
      currmax = std::max(tmp, currmax);
      outval[14][kr] = tmp;
    }
  for (int kr=0; kr<nc2; ++kr)
    outval[14][kr] = (currmax < 1.0e-6) ? 0 : outval[14][kr]*maxval/currmax;


  currmin = std::numeric_limits<double>::max();
  currmax = std::numeric_limits<double>::lowest();
  for (int kr=0; kr<nc2; ++kr)
    {
      currmin = std::min(nel[kr], currmin);
      currmax = std::max(nel[kr], currmax);
      outval[15][kr] = nel[kr];
    }
  for (int kr=0; kr<nc2; ++kr)
    outval[15][kr] = outval[15][kr]*maxval/currmax;

   currmin = std::numeric_limits<double>::max();
  currmax = std::numeric_limits<double>::lowest();
  for (int kr=0; kr<nc2; ++kr)
    {
      currmin = std::min(laplace[kr], currmin);
      currmax = std::max(laplace[kr], currmax);
      outval[16][kr] = laplace[kr];
    }
  for (int kr=0; kr<nc2; ++kr)
    outval[16][kr] = (currmax < 1.0e-9) ? 0 : outval[16][kr]*maxval/currmax;


  // Write to file
  out << ncell << "  " << ncell << "  " << "17" << std::endl;
  for (int kr=0; kr<nc2; ++kr)
    {
      for (int kh=0; kh<17; ++kh)
	out << outval[kh][kr] << " ";
      out << std::endl;
    }
}

