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
#include "GoTools/parametrization/PrPrmSymMeanValue.h"
#include "GoTools/parametrization/PrPrmEDDHLS.h"
#include "GoTools/parametrization/PrPrmMeanValue.h"
#include "GoTools/parametrization/PrPrmWachspress.h"
#include "GoTools/parametrization/PrPrmExperimental.h"
#include <fstream>
#include <memory>
using std::vector;
using std::cin;
using std::cout;
using std::endl;
using std::max;

int main()
{
	int numpnts;
    int numtrs;
    cin >> numpnts >> numtrs;

    vector<double> points(3*numpnts);
    vector<int> triangles(3*numtrs);
    for (int i = 0; i < numpnts; ++i) {
	cin >> points[3*i];
	cin >> points[3*i+1];
	cin >> points[3*i+2];
	points[3*i+2] = 0;
    }
    for (int i = 0; i < numtrs; ++i) {
	cin >> triangles[3*i];
	cin >> triangles[3*i+1];
	cin >> triangles[3*i+2];
    }
  
    shared_ptr<PrTriangulation_OP>
	pr_triang(new PrTriangulation_OP(&points[0],
					 numpnts,
					 &triangles[0],
					 numtrs));

//     cout << "The data of pr_triang is \n";
//     pr_triang->print(cout);
//     cout << "The information concerning pr_triang is \n";
//     pr_triang->printInfo(cout);

//     cout << "The nodes of pr_triang are \n";
//     pr_triang->printXYZNodes(cout);
//     cout << "The edges of pr_triang are \n";
//     pr_triang->printXYZEdges(cout);
//     cout << "The triangles of pr_triang are \n";
//     pr_triang->printXYZFaces(cout);



    //----------------------- Added testing code from Atgeirr ----------------

    for (int i = 0; i < numpnts; ++i) {
	if (pr_triang->isBoundary(i)) {
	    Vector3D p = pr_triang->get3dNode(i);
	    pr_triang->setU(i, p[0]);
	    pr_triang->setV(i, p[1]);
	} else {
	    pr_triang->setU(i, 0.0);
	    pr_triang->setV(i, 0.0);
	}
    }

    if (pr_triang->getNumNodes() - pr_triang->findNumBdyNodes() > 0) {
	PrPrmExperimental interior;
	//	PrPrmWachspress interior;
	//	PrPrmEDDHLS interior;
	// 	PrPrmSymMeanValue interior;
	//	PrPrmMeanValue interior;
	interior.attach(pr_triang);
	interior.setStartVectorKind(PrBARYCENTRE);
	interior.setBiCGTolerance(1.0e-6);
	cout << "Parametrizing interior..." << endl;
	interior.parametrize();
    }

    cout << "Saving parametrization to file..." << endl;
    std::ofstream pout("edges_out");
    pr_triang->printUVEdges(pout);

    double maxd = 0.0;
    for (int i = 0; i < numpnts; ++i) {
	Vector3D p = pr_triang->get3dNode(i);
	double du = fabs(p[0] - pr_triang->getU(i));
	double dv = fabs(p[1] - pr_triang->getV(i));
	maxd = max(maxd, du);
	maxd = max(maxd, dv);
//  	cout << "du = " << du << "    dv = " << dv << endl;
    }
    cout << "Max diff = " << maxd << endl;


    cout << "Quitting..." << endl;
}
