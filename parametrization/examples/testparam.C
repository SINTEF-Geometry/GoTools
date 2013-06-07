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

#include "GoTools/parametrization/PrTriangulation_OP.h"
#include "GoTools/parametrization/PrPlanarGraph_OP.h"
#include "GoTools/parametrization/PrParametrizeBdy.h"
#include "GoTools/parametrization/PrPrmShpPres.h"
#include <fstream>
#include <memory>
using std::cout;
using std::endl;

int main()
{
  cout << "Here is an explicit triangulation:\n";
  cout << "\n";
  cout << "   2-------3 \n";
  cout << "   |\\      | \n";
  cout << "   |\\\\     | \n";
  cout << "   | |\\    | \n";
  cout << "   | | \\   | \n";
  cout << "   | |  \\  | \n";
  cout << "   | 4-  \\ | \n";
  cout << "   |/  \\--\\| \n";
  cout << "   0-------1 \n";

  int numpnts=5, numtrs=4;
  vector<double> points(3*numpnts);
  vector<int> triangles(3*numtrs);
  points[0] = 0.0; points[1] = 0.0; points[2] = 0.0; 
  points[3] = 1.0; points[4] = 0.0; points[5] = 0.0; 
  points[6] = 0.0; points[7] = 1.0; points[8] = 0.0; 
  points[9] = 1.0; points[10] = 1.0; points[11] = 1.0; 
  points[12] = 0.25; points[13] = 0.25; points[14] = 0.5; 
  triangles[0] = 0; triangles[1] = 1; triangles[2] = 4;
  triangles[3] = 4; triangles[4] = 2; triangles[5] = 0;
  triangles[6] = 1; triangles[7] = 3; triangles[8] = 2;
  triangles[9] = 1; triangles[10] = 2; triangles[11] = 4;

  
  shared_ptr<PrTriangulation_OP>
      pr_triang(new PrTriangulation_OP(&points[0],
				       numpnts,
				       &triangles[0],
				       numtrs));

  cout << "The data of pr_triang is \n";
  pr_triang->print(cout);
  cout << "The information concerning pr_triang is \n";
  pr_triang->printInfo(cout);

  cout << "The nodes of pr_triang are \n";
  pr_triang->printXYZNodes(cout);
  cout << "The edges of pr_triang are \n";
  pr_triang->printXYZEdges(cout);
  cout << "The triangles of pr_triang are \n";
  pr_triang->printXYZFaces(cout);



  //----------------------- Added testing code from Atgeirr ----------------


  PrParametrizeBdy bdy;
  bdy.attach(pr_triang);
  bdy.setParamKind(PrCENTRIPETAL);

  cout << "Parametrizing boundary..." << endl;
  bdy.parametrize();

  if (pr_triang->getNumNodes() - pr_triang->findNumBdyNodes() > 0) {
      PrPrmShpPres interior;
      interior.attach(pr_triang);
      interior.setStartVectorKind(PrBARYCENTRE);
      interior.setBiCGTolerance(1.0e-6);
      cout << "Parametrizing interior..." << endl;
      interior.parametrize();
  }

  cout << "Saving parametrization to file..." << endl;
  std::ofstream pout("edges_out");
  pr_triang->printUVEdges(pout);

  cout << "Quitting..." << endl;



  return 0;
}
