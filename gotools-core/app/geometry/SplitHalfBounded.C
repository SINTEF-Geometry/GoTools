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

  ofstream outs(argv[2]);

  shared_ptr<Cylinder> under_surf = dynamic_pointer_cast<Cylinder>(in_surf->underlyingSurface());
  RectDomain rd = under_surf->parameterDomain();

  for (int i = 0; i < 2; ++i)
    {
      double from_upar = (i == 0) ? rd.umin() : M_PI;
      double to_upar = (i == 0) ? M_PI : rd.umax();
      double from_vpar = rd.vmin();
      double to_vpar = rd.vmax();
      cout << "Splitting to [" << from_upar << ", " << to_upar << "]x[" << from_vpar << ", " << to_vpar << "]" << endl;
      vector<shared_ptr<ParamSurface> > sub_surfs = in_surf->subSurfaces(from_upar, from_vpar, to_upar, to_vpar);
      cout << "Splitting completed, number of surfaces is " << sub_surfs.size() << endl;
      sub_surfs[0]->writeStandardHeader(outs);
      sub_surfs[0]->write(outs);
    }
  outs.close();

}
