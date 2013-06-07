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

#include <iostream>
#include <fstream>
#include "GoTools/geometry/SplineCurve.h"
#include "GoTools/geometry/ObjectHeader.h"

using std::cin;
using std::cout;
using std::endl;
using std::string;
using std::ifstream;
using std::ofstream;

using namespace Go;

int main(int argc, char* argv[])
{
    ifstream infile(argv[1]);
    ALWAYS_ERROR_IF(infile.bad(), "Wrong filename or corrupted file.");
    ALWAYS_ERROR_IF(argc != 6,
		"Arguments: Infile1 Infile2 order continuity Outfile");

    SplineCurve curve, other_curve;
    ObjectHeader header;

    infile >> header >> curve;
    other_curve = curve;
    ifstream infile2(argv[2]);
    ALWAYS_ERROR_IF(infile2.bad(), "Wrong filename or corrupted file.");
    infile2 >> header >> other_curve;

    int order = atoi(argv[3]);
    int order1 = curve.order();
    int other_order = other_curve.order();
    order = std::max(order, order1);
    order = std::max(order, other_order);
    if (order1 < order) {
      curve.raiseOrder(order - order1);
    } 
    if (other_order < order) {
      other_curve.raiseOrder(order - other_order);
    }

    // double par = 0.5*(curve.startparam()+curve.endparam());
    // if (/*curve.rational() && */curve.order() == curve.numCoefs())
    //   curve.insertKnot(par);
    // par = 0.5*(other_curve.startparam()+other_curve.endparam());
    // if (/*other_curve.rational() && */other_curve.order() == other_curve.numCoefs())
    //   other_curve.insertKnot(par);
#ifdef GEOMETRY_DEBUG
   // We insert knot into second cv. Rather case specific...
   double new_knot = other_curve.startparam() +
     0.00001*(other_curve.endparam()-other_curve.startparam());
   other_curve.insertKnot(new_knot);
   other_curve.insertKnot(new_knot);
#endif // GEOMETRY_DEBUG

    double dist = 0;
    curve.appendCurve(&other_curve, atoi(argv[4]), dist, true);
    cout << "Estimated difference between original and smooth curve: " << 
	dist << endl;

    ofstream outfile(argv[5]);
    outfile << header << curve;	

    return 0;

}
