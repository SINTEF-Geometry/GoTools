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

#include "GoTools/intersections/SfCvIntersector.h"
#include "GoTools/intersections/ParamCurveInt.h"
#include "GoTools/intersections/ParamSurfaceInt.h"
#include "GoTools/intersections/SplineCurveInt.h"
#include "GoTools/intersections/SplineSurfaceInt.h"
#include "GoTools/intersections/IntersectionPoint.h"
#include "GoTools/intersections/IntersectionCurve.h"
#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/geometry/SplineCurve.h"
#include "GoTools/geometry/SplineSurface.h"
#include <memory>
#include <fstream>
#include <iomanip>


using std::cout;
using std::endl;
using std::ifstream;
using std::ofstream;
using std::ostream;
using std::cerr;
using std::vector;
using std::setprecision;
using namespace Go;


int main(int argc, char** argv)
{
#ifndef _MSC_VER
    setenv("DEBUG", "1", 1);
    setenv("DEBUG_PAR", "1", 1);
    setenv("DEBUG_FINISH", "1", 1);
    setenv("DEBUG_IMPL", "1", 1);
    setenv("SUBDIV_SFCV", "1", 1);

    if (getenv("DEBUG") && (*getenv("DEBUG"))=='1')
	cout << "DEBUG=1" << endl;
    if (getenv("DEBUG_PAR") && (*getenv("DEBUG_PAR"))=='1')
	cout << "DEBUG_PAR=1" << endl;
    if (getenv("DEBUG_FINISH") && (*getenv("DEBUG_FINISH"))=='1')
	cout << "DEBUG_FINISH=1" << endl;
    if (getenv("DEBUG_IMPL") && (*getenv("DEBUG_IMPL"))=='1')
	cout << "DEBUG_IMPL=1" << endl;
    if (getenv("SUBDIV_SFCV") && (*getenv("SUBDIV_SFCV"))=='1')
	cout << "SUBDIV_SFCV=1" << endl;
#endif // _MSC_VER
 
    if (argc != 5 && argc != 6) {
	cout << "Usage: test_CvCvIntersector FileSf FileCv"
	     << " aepsge switch-order (OutputFile)" << endl;
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
    shared_ptr<ParamSurface> surf(new SplineSurface());
    surf->read(input1);
    input1.close();
    
    // Read the second curve from file
    ifstream input2(argv[2]);
    if (input2.bad()) {
	cerr << "File #2 error (no file or corrupt file specified)."
	     << std::endl;
	return 1;
    }
    header.read(input2);
    shared_ptr<ParamCurve> curve(new SplineCurve());
    curve->read(input2);
    input2.close();

    double aepsge;
    aepsge = atof(argv[3]);

    int seq;
    seq = atoi(argv[4]);

    ofstream tmpout;
    if (argc == 6)
	tmpout.open(argv[5]);
    ostream& out = (argc == 6) ? tmpout : cout;
    out << setprecision(8);

    // Setup the SplineSurfaceInt and SplineCurveInt objects
    shared_ptr<ParamGeomInt> ssurfint =
	shared_ptr<ParamGeomInt>(new SplineSurfaceInt (surf));
    shared_ptr<ParamGeomInt> scurveint =
	shared_ptr<ParamGeomInt>(new SplineCurveInt (curve));

    // Set up and run the intersection algorithm
    if (seq)
    {
	SfCvIntersector sfcvintersect(scurveint, ssurfint, aepsge);
	sfcvintersect.compute();

	// Get the results
	vector<shared_ptr<IntersectionPoint> > intpts;
	vector<shared_ptr<IntersectionCurve> > intcrv;
	sfcvintersect.getResult(intpts, intcrv);

	cout << "Number of points: " << intpts.size() << endl;
	cout << "Number of curves: " << intcrv.size() << endl;

	for (int ki = 0; ki < int(intpts.size()); ki++) {
	    double dist = intpts[ki]->getDist();
	    vector<double> par = intpts[ki]->getPar();
	    for (int kj = 0; kj < int(par.size()); kj++)
		cout << par[kj] << " ";
	    cout << ": dist = " << dist << endl;
	}

	// Write regression test data
	int npoints = (int)intpts.size();
	int ncurves = (int)intcrv.size();
	vector<double> par;
	out << "-1 ______________Case sf-cv intersector______________"
	    << endl;
	out << "0 Surf-Curve:    npoints = " << npoints
	    << "    ncurves = " << ncurves << endl;
	for (int i = 0; i < npoints; ++i) {
	    out << "1 Point: ";
	    par = intpts[i]->getPar();
	    for (int j=0; j < int(par.size()); j++) {
		out << par[j] << " ";
	    }
	    out << endl;
	}
	for (int i = 0; i < ncurves; ++i) {
	    out << "2 Curve:    Start = ";
	    int nguide = intcrv[i]->numGuidePoints(); 
	    shared_ptr<IntersectionPoint> startpt
		= intcrv[i]->getGuidePoint(0);
	    par = startpt->getPar();
	    for (int j=0; j < int(par.size()); j++) {
		out << par[j] << " ";
	    }
	    out << "   End = ";
	    shared_ptr<IntersectionPoint> endpt
		= intcrv[i]->getGuidePoint(nguide-1);
	    par = endpt->getPar();
	    for (int j=0; j < int(par.size()); j++) {
		out << par[j] << " ";
	    }
	    out << endl;
	}

    }
    else
    {
	SfCvIntersector sfcvintersect(ssurfint, scurveint, aepsge);
	sfcvintersect.compute();

	// Get the results
	vector<shared_ptr<IntersectionPoint> > intpts;
	vector<shared_ptr<IntersectionCurve> > intcrv;
	sfcvintersect.getResult(intpts, intcrv);
	cout << "Number of points: " << intpts.size() << endl;
	cout << "Number of curves: " << intcrv.size() << endl;

	for (int ki = 0; ki < int(intpts.size()); ki++) {
	    double dist = intpts[ki]->getDist();
	    vector<double> par = intpts[ki]->getPar();
	    for (int kj = 0; kj < int(par.size()); kj++)
		cout << par[kj] << " ";
	    cout << ": dist = " << dist << endl;
	}

	// Write regression test data
	int npoints = (int)intpts.size();
	int ncurves = (int)intcrv.size();
	vector<double> par;
	out << "-1 ______________Case sf-cv intersector______________"
	    << endl;
	out << "0 Surf-Curve:    npoints = " << npoints
	    << "    ncurves = " << ncurves << endl;
	for (int i = 0; i < npoints; ++i) {
	    out << "1 Point: ";
	    par = intpts[i]->getPar();
	    for (int j=0; j < int(par.size()); j++) {
		out << par[j] << " ";
	    }
	    out << endl;
	}
	for (int i = 0; i < ncurves; ++i) {
	    out << "2 Curve:    Start = ";
	    int nguide = intcrv[i]->numGuidePoints(); 
	    shared_ptr<IntersectionPoint> startpt
		= intcrv[i]->getGuidePoint(0);
	    par = startpt->getPar();
	    for (int j=0; j < int(par.size()); j++) {
		out << par[j] << " ";
	    }
	    out << "   End = ";
	    shared_ptr<IntersectionPoint> endpt
		= intcrv[i]->getGuidePoint(nguide-1);
	    par = endpt->getPar();
	    for (int j=0; j < int(par.size()); j++) {
		out << par[j] << " ";
	    }
	    out << endl;
	}

    }

    return 0;

}
