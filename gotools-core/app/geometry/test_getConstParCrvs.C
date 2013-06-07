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

#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/geometry/BoundedSurface.h"
#include "GoTools/geometry/SplineCurve.h"
#include "GoTools/geometry/BoundedUtils.h"
#include "GoTools/geometry/ObjectHeader.h"
#include <memory>
#include <fstream>

using namespace Go;
using std::ifstream;
using std::ofstream;
using std::vector;
using std::cout;
using std::endl;

int main(int argc, char** argv)
{

  if (argc != 5)
    {
      cout << "Usage: " << argv[0] << " surfaceinfile curvesoutfile nmb_crvs_u nmb_crvs_v" << endl;
      exit(-1);
    }

  ifstream filein(argv[1]);
  ALWAYS_ERROR_IF(filein.bad(), "Bad or no surface input filename");
  ObjectHeader head;
  filein >> head;
  if (head.classType() != SplineSurface::classType()) {
    THROW("Not a spline surface");
  }
  SplineSurface ss;
  filein >> ss;

  ofstream fileout(argv[2]);
  ALWAYS_ERROR_IF(fileout.bad(), "Bad curves output filename");

  int nmb_u = atoi(argv[3]);
  int nmb_v = atoi(argv[4]);
  vector<double> par_u, par_v;

  double step;
  double par;

  par = ss.startparam_u();
  if (nmb_u == 1)
    step = 0.0;
  else
    step = (ss.endparam_u()-par) / (double)(nmb_u-1);
  for (int i = 0; i < nmb_u; ++i, par += step)
    par_u.push_back(par);

  par = ss.startparam_v();
  if (nmb_v == 1)
    step = 0.0;
  else
    step = (ss.endparam_v()-par) / (double)(nmb_v-1);
  for (int i = 0; i < nmb_v; ++i, par += step)
    par_v.push_back(par);

  vector<shared_ptr<SplineCurve> > curves_u, curves_v;

  ss.getConstParamCurves(par_u, par_v, curves_u, curves_v);

  for (int i = 0; i < (int)curves_u.size(); ++i)
    {
      curves_u[i]->writeStandardHeader(fileout);
      curves_u[i]->write(fileout);
    }
  for (int i = 0; i < (int)curves_v.size(); ++i)
    {
      curves_v[i]->writeStandardHeader(fileout);
      curves_v[i]->write(fileout);
    }

  double eps = 1.0e-8;
  shared_ptr<BoundedSurface> bd_surf = 
    shared_ptr<BoundedSurface>(BoundedUtils::convertToBoundedSurface(ss, eps));

  size_t ki;
  for (ki=0; ki<par_u.size(); ++ki)
    {
      vector<shared_ptr<ParamCurve> > cvs = 
	bd_surf->constParamCurves(par_u[ki], false);
      for (size_t kj=0; kj<cvs.size(); ++kj)
	{
	  shared_ptr<SplineCurve> cc = 
	    shared_ptr<SplineCurve>(cvs[kj]->geometryCurve());
	  cc->writeStandardHeader(fileout);
	  cc->write(fileout);
	}
    }

  for (ki=0; ki<par_v.size(); ++ki)
    {
      vector<shared_ptr<ParamCurve> > cvs = 
	bd_surf->constParamCurves(par_v[ki], true);
      for (size_t kj=0; kj<cvs.size(); ++kj)
	{
	  shared_ptr<SplineCurve> cc = 
	    shared_ptr<SplineCurve>(cvs[kj]->geometryCurve());
	  cc->writeStandardHeader(fileout);
	  cc->write(fileout);
	}
    }
								      
}
