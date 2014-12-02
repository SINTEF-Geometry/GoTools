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

#include "GoTools/utils/config.h"
#include "GoTools/geometry/PointCloud.h"
#include "GoTools/utils/Array.h"
#include "GoTools/geometry/Utils.h"
#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/lrsplines2D/LRSplineSurface.h"
#include "GoTools/lrsplines2D/LRSurfApprox.h"
#include <iostream>
#include <fstream>
#include <string.h>

#define DEBUG

using namespace Go;
using std::vector;
using std::pair;

int main(int argc, char *argv[])
{
  if (argc != 7) {
    std::cout << "Usage: input rainfall data (.mat/.csv), nmb use, output info (.txt), output surface(.g2), tolerance, maximum number iterations " << std::endl;
    return -1;
  }

 // Read input arguments
  std::ifstream input_rain(argv[1]);
  int nmb_in = atoi(argv[2]);
  std::ofstream out_info(argv[3]);
  std::ofstream out_sf(argv[4]);
  double aepsge = atof(argv[5]);
  int max_iter = atoi(argv[6]);
  (void)out_info.precision(15);
  
  // Read rain data
  char location[40];
  vector<double> xyz;
  vector<double> rain;
  double val;
  int nmb_loc = 0;
  char curr[40];

  int ki, kj;
  input_rain >> location;
  while (!input_rain.eof())
    {
      nmb_loc++;
      for (ki=0; ki<3; ++ki)
	{
	  input_rain >> val;
	  xyz.push_back(val);
	}
      while (true)
	{
	  if (input_rain.eof())
	    break;
	  input_rain >> curr;
	  if (isdigit(curr[0]))
	    {
	      val = atof(curr);
	      rain.push_back(val);
	    }
	  else if (curr[0] == '-')
	    {
	      val = atof(curr);
	      rain.push_back(val);
	    }
	  else 
	    {
	      break;
	    }
	}
    }

  // Identify identical xy-coordinates
  vector<pair<int, int> > corr_sites;
  double site_tol = 1.0e-3;
  double site_tol2 = site_tol*site_tol;
  for (ki=0; ki<nmb_loc; ++ki)
    for (kj=ki+1; kj<nmb_loc; ++kj)
      {
	double tdist2 = (xyz[3*ki]-xyz[3*kj])*(xyz[3*ki]-xyz[3*kj])+
	  (xyz[3*ki+1]-xyz[3*kj+1])*(xyz[3*ki+1]-xyz[3*kj+1]);
	if (tdist2 < site_tol2)
	  corr_sites.push_back(std::make_pair(ki,kj));
      }
  std::cout << "Identical measurment stations: " << std::endl;
  for (ki=0; ki<(int)corr_sites.size(); ++ki)
    {
      std::cout << corr_sites[ki].first << ": (" << xyz[3*corr_sites[ki].first] << ",";
      std::cout << xyz[3*corr_sites[ki].first+1] << ",";
      std::cout << xyz[3*corr_sites[ki].first+2] << ") and ";
      std::cout << corr_sites[ki].second <<": (";
      std::cout << xyz[3*corr_sites[ki].second] << ",";
      std::cout << xyz[3*corr_sites[ki].second+1] << ",";
      std::cout << xyz[3*corr_sites[ki].second+2] << ")" << std::endl;
    }
  
  // Compute xy bounding box for rainfall data
  Point low(xyz[0], xyz[1]);
  Point high(xyz[0], xyz[1]);
  for (ki=1; ki<nmb_loc; ++ki)
    {
      for (kj=0; kj<2; ++kj)
	{
	  low[kj] = std::min(low[kj], xyz[3*ki+kj]);
	  high[kj] = std::max(high[kj], xyz[3*ki+kj]);
	}
    }

  // Make extended parameter domain to allow for some extrapolation
  Point mid = 0.5*(low + high);
  Point del = high - low;
  double domain[4];
  double fac = 0.01;
  domain[0] = low[0] - 0.1*del[0];
  domain[1] = high[0] + 0.1*del[0];
  domain[2] = low[1] - 0.1*del[1];
  domain[3] = high[1] + 0.1*del[1];
  domain[0] -= mid[0];
  domain[1] -= mid[0];
  domain[2] -= mid[1];
  domain[3] -= mid[1];
  for (ki=0; ki<4; ++ki)
    domain[ki] *= fac;
 
  int nmb_rain = (int)rain.size()/nmb_loc;
  std::cout << std::endl;
  std::cout << "Input points read and pre processed. Ready for ";
  std::cout << nmb_rain << " surface creations." << std::endl << std::endl;

  // Approximate rainfall
  for (kj=0; kj<nmb_rain; ++kj)
    {
      // Collect data
      vector<double> data;
      data.reserve(3*nmb_loc);
      int nmb_pts = 0;
      for (ki=0; ki<nmb_loc; ++ki)
	{
	  if (rain[ki*nmb_rain+kj] >= 0)
	    {
	      data.push_back(xyz[3*ki]);
	      data.push_back(xyz[3*ki+1]);
	      
	      size_t kr;
	      for (kr=0; kr<corr_sites.size(); ++kr)
		if (corr_sites[kr].first == ki || corr_sites[kr].second == ki)
		  break;
	      if (kr < corr_sites.size())
		data.push_back(std::max(rain[corr_sites[kr].first*nmb_rain+kj],
					rain[corr_sites[kr].second*nmb_rain+kj]));
	      else
		data.push_back(rain[ki*nmb_rain+kj]);
	      if (ki<nmb_in)
		nmb_pts++;
	    }
	}

      // Translate domain to origo to improve numerical accuracy
      int nmb_pts2 = (int)data.size()/3;   // nmb_pts <= nmb_loc
      for (ki=0; ki<nmb_pts2; ++ki)
      	for (int kr=0; kr<2; ++kr)
      	  {
      	    data[3*ki+kr] -= mid[kr];
	    data[3*ki+kr] *= fac;
      	  }

#ifdef DEBUG
      // Write translated surface and points to g2 format
      PointCloud3D cloud(data.begin(), nmb_pts);
      std::ofstream of2("translated_points.g2");
      cloud.writeStandardHeader(of2);
      cloud.write(of2);   
      PointCloud3D cloud2(data.begin()+nmb_pts*3, nmb_pts2-nmb_pts);
      cloud2.writeStandardHeader(of2);
      cloud2.write(of2);   
#endif
  
      std::cout << std::endl << "Starting observation nr. " << kj+1 << std::endl;

       // Initiate approximation
      int nmb_coef = 10;
      int order = 3; 
      bool init_tp = false;
      vector<double> data2(data.begin(), data.begin()+nmb_pts*3);
      LRSurfApprox approx(nmb_coef, order, nmb_coef, order, data2, 1, domain,  
			  aepsge,
			  // init_tp,
			  true, true);
      //approx.setVerbose(true);
      approx.setUseMBA(true);
      approx.addLowerConstraint(0.0);

      // Perform approximation
      double maxdist, avdist, avdist_total; // will be set below
      int nmb_out_eps;        // will be set below
      shared_ptr<LRSplineSurface> surf = approx.getApproxSurf(maxdist, avdist_total,
							      avdist, nmb_out_eps, 
							      max_iter);
#ifdef DEBUG
      std::cout << std::endl;
      std::cout << "Approximation nr " << kj+1 << " completed. " << std::endl;
      
      std::cout << "Total number of points: " << nmb_pts << std::endl;
      std::cout << "Number of elements: " << surf->numElements() << std::endl;
      std::cout << "Maximum distance: " << maxdist << std::endl;
      std::cout << "Average distance: " << avdist_total << std::endl;
      std::cout << "Average distance for points outside of the tolerance: " << avdist << std::endl;
      std::cout << "Number of points outside the tolerance: " << nmb_out_eps << std::endl;
#endif

#ifdef DEBUG
      std::ofstream of1("translated_sf_3d.g2");
      shared_ptr<LRSplineSurface> surf2(surf->clone());
      surf2->to3D();
      surf2->writeStandardHeader(of1);
      surf2->write(of1);
      LineCloud lines2 = surf2->getElementBds();
      lines2.writeStandardHeader(of1);
      lines2.write(of1);
#endif



      // Translate surface back to the initial domain
      double fac2 = 1.0/fac;
      surf->setParameterDomain(fac2*domain[0] + mid[0], fac2*domain[1] + mid[0],
      			       fac2*domain[2] + mid[1], fac2*domain[3] + mid[1]);
      // surf->setParameterDomain(domain[0] + mid[0], domain[1] + mid[0],
      // 			       domain[2] + mid[1], domain[3] + mid[1]);

      // Write to output files
      std::cout << std::endl << "Accuracy with regard to input rain values" << std::endl;
      
      // Write accuracy information
      double minval = 1.0e8;
      double maxval = -1.0e8;
      double maxerr = -1.0;
      double xerr, yerr;

      out_info << std::endl << std::endl;
      out_info << "Observation nr " << kj << ".Format: x, y, z, computed rain, difference from observed " << std::endl;
      out_info << "=========================================================================" << std::endl;
      for (ki=0; ki<nmb_loc; ++ki)
	{
	  if (rain[ki*nmb_rain+kj] >= 0)
	    {
	      Point pos;
	      surf->point(pos, xyz[3*ki], xyz[3*ki+1]);
	      out_info << xyz[3*ki] << " " << xyz[3*ki+1] << " ";
	      out_info << xyz[3*ki+2] << " " << pos[0] << " ";
	      out_info << rain[ki*nmb_rain+kj]-pos[0] <<std::endl;
	      minval = std::min(minval, pos[0]);
	      maxval = std::max(maxval, pos[0]);
	      if (fabs(rain[ki*nmb_rain+kj]-pos[0]) > maxerr)
		{
		  maxerr = fabs(rain[ki*nmb_rain+kj]-pos[0]);
		  xerr = xyz[3*ki];
		  yerr = xyz[3*ki+1];
		}
	    }
	}
      out_info << std::endl << "Largest rain value: " << maxval << std::endl;
      out_info << "Largest difference: " << maxerr << " at (" << xerr << ", ";
      out_info << yerr << ")" << std::endl;

      std::cout << "Largest rain value: " << maxval << std::endl;
      //std::cout << "Smallest rain value: " << minval << std::endl;
      std::cout << "Largest difference: " << maxerr << " at (" << xerr << ", ";
      std::cout << yerr << ")" << std::endl;
	
      // Write surface
      surf->writeStandardHeader(out_sf);
      surf->write(out_sf);
     }
}
