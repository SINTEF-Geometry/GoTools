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

//#define DEBUG

using namespace Go;
using std::vector;
using std::pair;

int main(int argc, char *argv[])
{

  if (argc != 3)
    {
      std::cout << "Input parameters: input file (.csv), output file " << std::endl;
      return -1;
    }

 // Define parameters
  std::ifstream input_rain(argv[1]);
  std::ofstream out_rain(argv[2]);
  double aepsge = 0.05;
  int max_iter = 16;
  
  // Read rain data
  char location[40];
  vector<double> xyz;
  vector<double> rain;
  double val;
  int nmb_loc = 0;
  char curr[40];

  int ki, kj, kh;
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
  domain[0] = low[0] - 0.1*del[0];
  domain[1] = high[0] + 0.1*del[0];
  domain[2] = low[1] - 0.1*del[1];
  domain[3] = high[1] + 0.1*del[1];
  domain[0] -= mid[0];
  domain[1] -= mid[0];
  domain[2] -= mid[1];
  domain[3] -= mid[1];

  // Approximate rainfall
  int nmb_rain = (int)rain.size()/nmb_loc;
  vector<double> computed_rain(nmb_rain*nmb_loc);
  for (kj=0; kj<nmb_rain; ++kj)
    {
      std::cout << "Rainfall observation nr " << kj+1 << std::endl;

      for (kh=0; kh<nmb_loc; ++kh)
	{
	  // Collect data
	  vector<double> data;
	  data.reserve(3*(nmb_loc-1));
	  for (ki=0; ki<nmb_loc; ++ki)
	    {
	      if (ki == kh)
		continue;  // This station is turned off

	      if (rain[ki*nmb_rain+kj] >= 0)
		{
	      
		  size_t kr;
		  for (kr=0; kr<corr_sites.size(); ++kr)
		    {
		      if (corr_sites[kr].first == kh || corr_sites[kr].second == kh)
			break;  // One of the double sites is ignored. Use none
		      if (corr_sites[kr].first == ki || corr_sites[kr].second == ki)
			break;
		    }
		  if (kr < corr_sites.size())
		    {
		      if (corr_sites[kr].first == kh || corr_sites[kr].second == kh);
		      else if (corr_sites[kr].first == kh)
			{
			  // Use the first occurance of the location and the maximum
			  // rain value
			  data.push_back(xyz[3*ki]);
			  data.push_back(xyz[3*ki+1]);
			  data.push_back(std::max(rain[corr_sites[kr].first*nmb_rain+kj],
						  rain[corr_sites[kr].second*nmb_rain+kj]));
			}
		    }
		  else
		    {
		      data.push_back(xyz[3*ki]);
		      data.push_back(xyz[3*ki+1]);
		      data.push_back(rain[ki*nmb_rain+kj]);
		    }
		}
	    }

	  // Translate domain to origo to improve numerical accuracy
	  int nmb_pts = (int)data.size()/3;   // nmb_pts <= nmb_loc
	  for (ki=0; ki<nmb_pts; ++ki)
	    for (int kr=0; kr<2; ++kr)
	      {
		data[3*ki+kr] -= mid[kr];
	      }

	  // Initiate approximation
	  int nmb_coef = 10;
	  int order = 3; 
	  bool init_tp = false;
	  LRSurfApprox approx(nmb_coef, order, nmb_coef, order, data, 1, domain,  
			      aepsge, init_tp, true, true);
	  //approx.setVerbose(true);
	  approx.setUseMBA(true);
	  approx.addLowerConstraint(0.0);

	  // Perform approximation
	  double maxdist, avdist, avdist_total; // will be set below
	  int nmb_out_eps;        // will be set below
	  shared_ptr<LRSplineSurface> surf = approx.getApproxSurf(maxdist, avdist_total,
								  avdist, nmb_out_eps, 
								  max_iter);
      
	  surf->setParameterDomain(domain[0] + mid[0], domain[1] + mid[0],
				   domain[2] + mid[1], domain[3] + mid[1]);

	  // Collect output
	  if (rain[kh*nmb_rain+kj] >= 0)
	    {
	      Point pos;
	      surf->point(pos, xyz[3*kh], xyz[3*kh+1]);

	      size_t kr;
	      double rain_val = rain[kh*nmb_rain+kj];
	      for (kr=0; kr<corr_sites.size(); ++kr)
		{
		  if (corr_sites[kr].first == kh || corr_sites[kr].second == kh)
		    {
		      // One of the double sites is ignored. Compare with the
		      // maximum number
		      rain_val = std::max(rain[corr_sites[kr].first*nmb_rain+kj],
					  rain[corr_sites[kr].second*nmb_rain+kj]);
		      break;  
		    }
		}
	      computed_rain[kh*nmb_rain+kj] = pos[0] - rain_val;
	    }
	  else
	    computed_rain[kh*nmb_rain+kj] = 0.0;
	}
    }


  // Write output to file
  for (ki=0; ki<nmb_loc; ++ki)
    {
      for (kj=0; kj<nmb_rain; ++kj)
	out_rain << computed_rain[ki*nmb_rain+kj] << "  ";
      out_rain << std::endl;
    }
}

