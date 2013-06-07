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
#include "GoTools/geometry/GeometryTools.h"

using namespace Go;
using namespace std;

int main()
{
     ObjectHeader h;
     SplineSurface s;
     cin >> h >> s;

     int peru = GeometryTools::analyzePeriodicityDerivs(s, 0, s.basis_u().order()-2);
     cout << "U-periodicity is (derivative-based): " << peru << endl;
     int perku = GeometryTools::analyzePeriodicity(s, 0);
     cout << "U-periodicity is (knot and control point-based): " << perku << endl;
     int perv = GeometryTools::analyzePeriodicityDerivs(s, 1, s.basis_v().order()-2);
     cout << "V-periodicity is (derivative-based): " << perv << endl;
     int perkv = GeometryTools::analyzePeriodicity(s, 1);
     cout << "V-periodicity is (knot and control point-based): " << perkv << endl;
}
