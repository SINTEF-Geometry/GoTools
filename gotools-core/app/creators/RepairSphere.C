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
#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/geometry/BoundedSurface.h"
#include "GoTools/geometry/GoTools.h"


using namespace std;
using namespace Go;


int main(int argc, char** argv)
{
  GoTools::init();

  if (argc != 3) {
    cout << "Usage:  " << argv[0] << " infile outfile" << endl;
    return 1;
  }

  ifstream ins(argv[1]);
  ObjectHeader header;
  BoundedSurface* in_surf = new BoundedSurface();
  ins >> header;
  ins >> *in_surf;
  ins.close();

  /*
  shared_ptr<Sphere> sphere = dynamic_pointer_cast<Sphere>(in_surf->underlyingSurface());

  CurveLoop loop = in_surf->outerBoundaryLoop();
  for (int i = 0; i < loop.size(); ++i)
    {
      shared_ptr<ParamCurve> curve = loop[i];
      shared_ptr<CurveOnSurface> cos = dynamic_pointer_cast<CurveOnSurface>(curve);
      shared_ptr<ParamCurve> geo_curve = cos->spaceCurve();
      shared_ptr<Circle> curve_circ = dynamic_pointer_cast<Circle>(geo_curve);
      if (curve_circ.get())
	{
	  cout << "Got circle at i = " << i << endl;
	  shared_ptr<ParamCurve> par_curve = sphere->getElementaryParamCurve(curve_circ.get(), 1.0e-3);
	  cos->setParameterCurve(par_curve);
	}
    }
  */

  ofstream outs(argv[2]);
  in_surf->writeStandardHeader(outs);
  in_surf->write(outs);
  outs.close();
}
