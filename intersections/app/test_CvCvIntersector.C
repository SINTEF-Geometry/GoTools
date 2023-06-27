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

#include "GoTools/intersections/CvCvIntersector.h"
#include "GoTools/intersections/ParamCurveInt.h"
#include "GoTools/intersections/SplineCurveInt.h"
#include "GoTools/intersections/IntersectionPoint.h"
#include "GoTools/intersections/IntersectionCurve.h"
#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/geometry/SplineCurve.h"
#include <fstream>
#include <iomanip>


using namespace std;
using namespace Go;


int main(int argc, char** argv)
{
 
    if (argc != 4) {
	cout << "Usage: test_CvCvIntersector FileCv1 FileCv2 aepsge"
	     << endl;
	return 0;
    }


    ObjectHeader header;

    // Read the first curve from file
    ifstream input1(argv[1]);
    if (input1.bad()) {
	cerr << "File #1 error (no file or corrupt file specified)."
	     << std::endl;
	return 1;
    }
    header.read(input1);
    shared_ptr<ParamCurve> curve1(new SplineCurve());
    curve1->read(input1);
    input1.close();
    
    // Read the second curve from file
    ifstream input2(argv[2]);
    if (input2.bad()) {
	cerr << "File #2 error (no file or corrupt file specified)."
	     << std::endl;
	return 1;
    }
    header.read(input2);
    shared_ptr<ParamCurve> curve2(new SplineCurve());
    curve2->read(input2);
    input2.close();

    double aepsge;
    aepsge = atof(argv[3]);

//     cout << "\nFile : " << argv[2] << " Parameter range u: "
// 	 << curve1->startparam_u() <<" " << curve1->endparam_u()
// 	 << " Parameter range v: " << curve1->startparam_v() <<" " 
// 	 << curve1->endparam_v();

//     cout << "\nFile : " << argv[2] << " Parameter range u: "
// 	 << curve2->startparam_u() <<" " << curve2->endparam_u()
// 	 << " Parameter range v: " << curve2->startparam_v() <<" " 
// 	 << curve2->endparam_v() << endl;

    shared_ptr<ParamGeomInt> scurveint1 =
	shared_ptr<ParamGeomInt>(new SplineCurveInt (curve1));
    shared_ptr<ParamGeomInt> scurveint2 =
	shared_ptr<ParamGeomInt>(new SplineCurveInt (curve2));

    CvCvIntersector cvcvintersect (scurveint1, scurveint2, aepsge);
    cvcvintersect.compute();

    std::vector<shared_ptr<IntersectionPoint> > intpts;
    std::vector<shared_ptr<IntersectionCurve> > intcrv;
    cvcvintersect.getResult(intpts, intcrv);
    printf("Number of points: %d \n", int(intpts.size()));
    printf("Number of curves: %d \n", int(intcrv.size()));

    if (intpts.size() > 0)
      {
	std::ofstream of("intpts.g2");
	of << "400 1 0 4 255 0 0 255" << std::endl;
	of << intpts.size() << std::endl;
	int ki, kj;
	for (ki=0; ki < int(intpts.size()); ki++) {
	  std::vector<double> par = intpts[ki]->getPar();
	  for (kj=0; kj<int(par.size()); kj++)
	    std::cout << par[kj] << " ";
	  std::cout << std::endl;

	  Point pos = intpts[ki]->getPoint();
	  of << pos;
	  if (pos.dimension() == 2)
	    of << " " << 0.0;
	  of << std::endl;
	}
      }

    if (intcrv.size() > 0)
      {
	std::ofstream of("intcvs.g2");
	for (size_t ki=0; ki<intcrv.size(); ++ki)
	  {
	    shared_ptr<ParamCurve> splcrv = intcrv[ki]->getCurve();
	    if (splcrv.get()) {
	    of << "100 1 0 4 50 205 50 255" << endl;
	    splcrv->write(of);

	    double start, end;
	    intcrv[ki]->getParamSpan(start, end);
	    Point pos, tan;
	    intcrv[ki]->evaluateAt(start, pos, tan);
	    of << "400 1 0 4 255 0 0 255" << endl;
	    of << "1" << std::endl;
	    of << pos << std::endl;
	    intcrv[ki]->evaluateAt(0.5*(start+end), pos, tan);
	    of << "400 1 0 4 255 0 0 255" << endl;
	    of << "1" << std::endl;
	    of << pos << std::endl;
	    intcrv[ki]->evaluateAt(end, pos, tan);
	    of << "400 1 0 4 255 0 0 255" << endl;
	    of << "1" << std::endl;
	    of << pos << std::endl;
	    }
	  }
	}

	
    return 0;
}
