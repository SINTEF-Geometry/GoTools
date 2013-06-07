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

#include <fstream>
#include "GoTools/trivariate/SplineVolume.h"
#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/utils/errormacros.h"
#include "GoTools/geometry/BsplineBasis.h"




using namespace std;
using namespace Go;

void write_basis(ofstream &os, double start, double end, int ord, int nmb)
{
  os << nmb << " " << ord << endl;

  double step = (end-start)/(double)(nmb-ord+1);
  double val = start;
  for (int i = 0; i < nmb + ord - 1; ++i)
    {
      os << val << " ";
      if (i >= ord-1 && i < nmb) val += step;
    }
  os << val << endl;
}


int main(int argc, char* argv[] )
{
  ALWAYS_ERROR_IF(argc != 14, "Usage: " << argv[0]
		  << " theta_start theta_end theta_order theta_pts"
		  << " phi_start phi_end phi_order phi_pts"
		  << " r_start r_end r_order r_pts"
		  << " outfile");
  ofstream os(argv[13]);
  ALWAYS_ERROR_IF(os.bad(), "Bad or no output filename");

  double t_start = atof(argv[1]);
  double t_end = atof(argv[2]);
  int t_ord = atoi(argv[3]);
  int t_nmb = atoi(argv[4]);
  double p_start = atof(argv[5]);
  double p_end = atof(argv[6]);
  int p_ord = atoi(argv[7]);
  int p_nmb = atoi(argv[8]);
  double r_start = atof(argv[9]);
  double r_end = atof(argv[10]);
  int r_ord = atoi(argv[11]);
  int r_nmb = atoi(argv[12]);

  os << "700 1 0 0" << endl;
  os << "3 0" << endl;

  write_basis(os,t_start, t_end, t_ord, t_nmb);
  write_basis(os,p_start, p_end, p_ord, p_nmb);
  write_basis(os,r_start, r_end, r_ord, r_nmb);

  double t_step = (t_end-t_start)/(double)(t_nmb-1);
  double p_step = (p_end-p_start)/(double)(p_nmb-1);
  double r_step = (r_end-r_start)/(double)(r_nmb-1);

  for (int i = 0; i < r_nmb; ++i)
    {
      double r = r_start + r_step * (double)(i);
      for (int j = 0; j < p_nmb; ++j)
	{
	  double p = p_start + p_step * (double)(j);
	  double cp = cos(p);
	  double sp = sin(p);
	  for (int k = 0; k < t_nmb; ++k)
	    {
	      double t = t_start + t_step * (double)(k);
	      os << r*cp*cos(t) << " " << r*cp*sin(t) << " " << r*sp << endl;
	    }
	}
    }
}
