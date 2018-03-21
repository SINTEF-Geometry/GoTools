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
#include "GoTools/geometry/Line.h"
#include "GoTools/geometry/Circle.h"
#include "GoTools/geometry/Ellipse.h"
#include "GoTools/geometry/Hyperbola.h"
#include "GoTools/geometry/Parabola.h"
#include "GoTools/geometry/SplineCurve.h"
#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/utils/errormacros.h"
#include <memory>

using namespace Go;
using std::ifstream;
using std::ofstream;
using std::endl;
using std::vector;

int main(int argc, char* argv[] )
{
    if (argc != 3)
    {
	std::cout << "Usage: " << argv[0]
		  << " input_file output_file" << endl;
	return -1;
    }

    // Open input surface file
    ifstream is(argv[1]);
    ofstream os(argv[2]);
    if (is.bad())
    {
	std::cout << "Bad or no input filename" << std::endl;
	return -1;
    }

    // Read surface from file
    ObjectHeader head;
    shared_ptr<ElementaryCurve> crv;
    is >> head;
    if (head.classType() == Class_Line)
      crv = shared_ptr<ElementaryCurve>(new Line());
    else if (head.classType() == Class_Circle)
      crv = shared_ptr<ElementaryCurve>(new Circle());
    else if (head.classType() == Class_Ellipse)
      crv = shared_ptr<ElementaryCurve>(new Ellipse());
    else if (head.classType() == Class_Hyperbola)
      crv = shared_ptr<ElementaryCurve>(new Hyperbola());
    else if (head.classType() == Class_Parabola)
      crv = shared_ptr<ElementaryCurve>(new Parabola());
    if (!crv.get())
      {
	std::cout << "Input file did not contain an elementary curve" << std::endl;
	return -1;
      }
    crv->read(is);

    crv->writeStandardHeader(os);
    crv->write(os);

    double t1 = crv->startparam();
    double t2 = crv->endparam();
    std::cout << "Parameter bounds: " << t1 << ", " << t2 << std::endl;

    Point p1, p2;
    vector<Point> der1(3);
    vector<Point> der2(3);
    crv->point(p1, 0.1*t1+0.9*t2);
    crv->point(der1, 0.1*t1+0.9*t2, 2);
    std::cout << "Der1: " << der1[0] << std::endl
	      << der1[1] << std::endl
	      << der1[2] << std::endl;

    double tpar0, dist0;
    Point clo_pt0;
    crv->closestPoint(p1, t1, t2, tpar0, clo_pt0, dist0);
    std::cout << "Par: " << tpar0 << ", dist: " << dist0 
	      << ", closest point: " << clo_pt0 << std::endl;

    crv->reverseParameterDirection();
    crv->point(p2, 0.1*t1+0.9*t2);
    crv->point(der2, 0.1*t1+0.9*t2, 2);
    std::cout << std::endl << "Der2: " << der2[0] << std::endl
	      << der2[1] << std::endl
	      << der2[2] << std::endl;
    
    os << "400 1 0 4 255 0 0 255" << std::endl;
    os << "1" << std::endl;
    os << p1 << std::endl;

    os << "400 1 0 4 0 255 0 255" << std::endl;
    os << "1" << std::endl;
    os << p2 << std::endl;

    double tpar2, dist2;
    Point clo_pt2;
    crv->closestPoint(p1, t1, t2, tpar2, clo_pt2, dist2);
    std::cout << "Par: " << tpar2 << ", dist: " << dist2 
	      << ", closest point: " << clo_pt2 << std::endl;


    shared_ptr<SplineCurve> spline_cv(crv->geometryCurve());
 
    spline_cv->writeStandardHeader(os);
    spline_cv->write(os);

     crv->setParameterInterval(t1, 2*t2);
    double t3, t4;
    t3 = crv->startparam();
    t4 = crv->endparam();
    std::cout << std::endl << "Parameter interval: " << t3 << ", " << t4 << std::endl;

    vector<Point> der3(3);
    crv->point(der3, 0.1*t3+0.9*t4, 2);
    std::cout << "Der3: " << der3[0] << std::endl
	      << der3[1] << std::endl
	      << der3[2] << std::endl;

    crv->setParameterInterval(t1+t2, t1+4*fabs(t2-1));
    t3 = crv->startparam();
    t4 = crv->endparam();
    std::cout << std::endl << "Parameter interval: " << t3 << ", " << t4 << std::endl;

    crv->point(der3, 0.1*t3+0.9*t4, 2);
    std::cout << "Der3(2): " << der3[0] << std::endl
	      << der3[1] << std::endl
	      << der3[2] << std::endl;

    shared_ptr<ElementaryCurve> crv2(crv->subCurve(t3, 0.1*t3+0.9*t4));
    double t5, t6;
    t5 = crv2->startparam();
    t6 = crv2->endparam();
    std::cout << std::endl << "Parameter interval: " << t5 << ", " << t6 << std::endl;
    vector<Point> der5(3);
    crv2->point(der5, t6, 2);
    std::cout << "Der5: " << der5[0] << std::endl
	      << der5[1] << std::endl
	      << der5[2] << std::endl;
    crv2->writeStandardHeader(os);
    crv2->write(os);

    shared_ptr<SplineCurve> spline_cv2(crv2->geometryCurve());
    vector<Point> der4(3);
    spline_cv2->point(der4, 0.1*t3+0.9*t4, 2);
    std::cout << std::endl << "Der4: " << der4[0] << std::endl
	      << der4[1] << std::endl
	      << der4[2] << std::endl;
    spline_cv2->writeStandardHeader(os);
    spline_cv2->write(os);
 
    double mid = 0.5*(spline_cv2->startparam()+spline_cv2->endparam());
    Point p3 = spline_cv2->ParamCurve::point(mid);
    std::cout << std::endl << "Par: " << mid << ", pos: " << p3 << std::endl;
    double tpar, dist;
    Point clo_pt;
    crv->closestPoint(p3, t3, t4, tpar, clo_pt, dist);
    std::cout << "Par: " << tpar << ", dist: " << dist 
	      << ", closest point: " << clo_pt << std::endl;

    shared_ptr<ElementaryCurve> crv3(crv->subCurve(0.1*t3+0.9*t4, t4));
    crv3->writeStandardHeader(os);
    crv3->write(os);

    crv2->appendCurve(crv3.get(), false);
    crv2->writeStandardHeader(os);
    crv2->write(os);
    
     return 0;
}
