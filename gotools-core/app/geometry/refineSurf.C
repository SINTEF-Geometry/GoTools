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
#include "GoTools/geometry/ObjectHeader.h"
#include <fstream>
#include <vector>

using namespace Go;
using namespace std;

int main(int argc, char** argv)
{
  if (argc < 3) {
      cerr << "Usage: " << argv[0]
	   << " inputfile outputfile [max_coefs_u max_coefs_v]" << endl;
      return 1;
  }

  ifstream in(argv[1]);
  ofstream out(argv[2]);

  if (!in || !out) {
    cout << "Bad file(s) or filename(s)." << endl;
    return 1;
  }

  ObjectHeader oh;
  SplineSurface sf;

  in >> oh >> sf;


  int m = sf.numCoefs_v() - sf.order_v() + 1;
  int n = sf.numCoefs_u() - sf.order_u() + 1;
  if (argc >= 5) {
      // Note the weird order (v then u)
      m = min(atoi(argv[4])-sf.numCoefs_v(), m);
      n = min(atoi(argv[3])-sf.numCoefs_u(), n);
  }
  int i;
  vector<double> newknots_v;
  vector<double> newknots_u;
  for (i = 0; i < m; ++i) {
    vector<double>::const_iterator it = sf.basis_v().begin();
    double newknot = 0.5*it[sf.order_v()+i-1] + 0.5*it[sf.order_v()+i];
    newknots_v.push_back(newknot);
  }
  for (i = 0; i < n; ++i) {
    vector<double>::const_iterator it = sf.basis_u().begin();
    double newknot = 0.5*it[sf.order_u()+i-1] + 0.5*it[sf.order_u()+i];
    newknots_u.push_back(newknot);
  }

  sf.insertKnot_v(newknots_v);
  sf.insertKnot_u(newknots_u);

  out << oh << sf;
  return 0;
}


