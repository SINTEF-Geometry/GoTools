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
#include "GoTools/geometry/Plane.h"
#include "GoTools/geometry/Cylinder.h"
#include "GoTools/geometry/Cone.h"
#include "GoTools/geometry/Sphere.h"
#include "GoTools/geometry/Torus.h"
#include "GoTools/geometry/Disc.h"
#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/geometry/ElementaryCurve.h"
#include "GoTools/geometry/SplineCurve.h"
#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/geometry/RectDomain.h"
#include "GoTools/creators/CurveCreators.h"
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

    double eps = 1.0e-6;

    // Read surface from file
    ObjectHeader head;
    shared_ptr<ElementarySurface> surf;
    is >> head;
    if (head.classType() == Class_Plane)
      surf = shared_ptr<ElementarySurface>(new Plane());
    else if (head.classType() == Class_Cylinder)
      surf = shared_ptr<ElementarySurface>(new Cylinder());
    else if (head.classType() == Class_Cone)
      surf = shared_ptr<ElementarySurface>(new Cone());
    else if (head.classType() == Class_Sphere)
      surf = shared_ptr<ElementarySurface>(new Sphere());
    else if (head.classType() == Class_Torus)
      surf = shared_ptr<ElementarySurface>(new Torus());
    else if (head.classType() == Class_Disc)
      surf = shared_ptr<ElementarySurface>(new Disc());
    if (!surf.get())
      {
	std::cout << "Input file did not contain an elementary surface" << std::endl;
	return -1;
      }
    surf->read(is);

    if (surf->isBounded())
      {
	surf->writeStandardHeader(os);
	surf->write(os);
      }

    RectDomain dom1 = surf->containingDomain();
    std::cout << "Parameter interval: " << dom1.umin() << ", " << dom1.umax() 
	      << ", " << dom1.vmin() << ", " << dom1.vmax() << std::endl;

    DirectionCone cone = surf->normalCone();
    std::cout << "Cone, dir: " << cone.centre() << ", angle: "
	      << cone.angle() << ", gtpi: " << cone.greaterThanPi() << std::endl;
    if (!surf->isBounded())
      {
	if (surf->instanceType() == Class_Plane)
	  surf->setParameterBounds(0.0, 0.0, 0.5, 0.5);
	else
	  surf->setParameterBounds(dom1.umin(), 0.0, dom1.umax(), 0.05);
	
	dom1 = surf->containingDomain();
	std::cout << "Parameter interval: " << dom1.umin() << ", " << dom1.umax() 
		  << ", " << dom1.vmin() << ", " << dom1.vmax() << std::endl;
      }

    Point p1, p2, n1, n2;
    vector<Point> der1(3);
    vector<Point> der2(3);
    surf->point(p1, 
		0.1*dom1.umin()+0.9*dom1.umax(), 0.5*(dom1.vmin()+dom1.vmax()));
    surf->normal(n1, 
		0.1*dom1.umin()+0.9*dom1.umax(), 0.5*(dom1.vmin()+dom1.vmax()));
     surf->point(der1, 
		0.1*dom1.umin()+0.9*dom1.umax(), 0.5*(dom1.vmin()+dom1.vmax()),
		1);
    std::cout << "Der1: " << der1[0] << std::endl
	      << der1[1] << std::endl
	      << der1[2] << std::endl;
    std::cout << "Normal: " << n1 << std::endl;

    double upar0, vpar0, dist0;
    Point clo_pt0;
    surf->closestPoint(p1, upar0, vpar0, clo_pt0, dist0, eps);
    std::cout << "Par: (" << upar0 << ", "<< vpar0 << "), dist: " << dist0 
	      << ", closest point: " << clo_pt0 << std::endl;

    surf->swapParameterDirection();
    RectDomain dom2 = surf->containingDomain();
    surf->point(p2, 
		0.1*dom2.umin()+0.9*dom2.umax(), 0.5*(dom2.vmin()+dom2.vmax()));
    surf->normal(n2, 
		0.1*dom2.umin()+0.9*dom2.umax(), 0.5*(dom2.vmin()+dom2.vmax()));
    surf->point(der2, 
		0.1*dom2.umin()+0.9*dom2.umax(), 0.5*(dom2.vmin()+dom2.vmax()),
		1);
    std::cout << std::endl << "Der2: " << der2[0] << std::endl
	      << der2[1] << std::endl
	      << der2[2] << std::endl;
    std::cout << "Normal: " << n2 << std::endl;
    
    os << "400 1 0 4 255 0 0 255" << std::endl;
    os << "1" << std::endl;
    os << p1 << std::endl;

    os << "400 1 0 4 0 255 0 255" << std::endl;
    os << "1" << std::endl;
    os << p2 << std::endl;

    double upar2, vpar2, dist2;
    Point clo_pt2;
    surf->closestPoint(p1, upar2, vpar2, clo_pt2, dist2, eps);
    std::cout << "Par: (" << upar2 << ", "<< vpar2 << "), dist: " << dist2 
	      << ", closest point: " << clo_pt2 << std::endl;

    shared_ptr<SplineSurface> spline_sf(surf->geometrySurface());
 
    if (spline_sf.get())
      {
	spline_sf->writeStandardHeader(os);
	spline_sf->write(os);
      }
    
    vector<Point> der3_0(3);
    surf->point(der3_0, 0.5*(dom2.umin()+dom2.umax()),
		0.5*(dom2.vmin()+dom2.vmax()), 1);
    std::cout << "Der3: " << der3_0[0] << std::endl
    	      << der3_0[1] << std::endl
    	      << der3_0[2] << std::endl;

    double del1 = dom2.umax()-dom2.umin();
    double del2 = dom2.vmax()-dom2.vmin();
    surf->setParameterDomain(dom2.umin(), dom2.umax()+del1, 
			       dom2.vmin()-0.5*del2, dom2.vmax()+0.5*del2);
    RectDomain dom3 = surf->containingDomain();
    std::cout << "Parameter interval: " << dom3.umin() << ", " << dom3.umax() 
	      << ", " << dom3.vmin() << ", " << dom3.vmax() << std::endl;

    vector<Point> der3(3);
    surf->point(der3, 0.5*(dom3.umin()+dom3.umax()),
		0.5*(dom3.vmin()+dom3.vmax()), 1);
    std::cout << "Der3: " << der3[0] << std::endl
    	      << der3[1] << std::endl
    	      << der3[2] << std::endl;
							
    surf->writeStandardHeader(os);
    surf->write(os);

    vector<CurveLoop> bd_loops = surf->allBoundaryLoops();
    for (size_t ki=0; ki<bd_loops.size(); ++ki)
      {
	for (int kj=0; kj<bd_loops[ki].size(); ++kj)
	  {
	    shared_ptr<ParamCurve> curr = bd_loops[ki][kj];
	    curr->writeStandardHeader(os);
	    curr->write(os);
	    os << "400 1 0 4 100 100 55 255" << std::endl;
	    os << "1" << std::endl;
	    os << curr->point(curr->startparam()) << std::endl;
	  }
      }

    del1 = dom3.umax() - dom3.umin();
    del2 = dom3.vmax() - dom3.vmin();
    vector<shared_ptr<ParamSurface> > sfs2 = 
      surf->subSurfaces(dom3.umin(), dom3.vmin(), dom3.umax(), 
			dom3.vmin()+0.4*del2);
    vector<shared_ptr<ParamSurface> > sfs3 = 
      surf->subSurfaces(dom3.umin()+0.2*del1, dom3.vmin(), 
			dom3.umax()-0.1*del1, dom3.vmin()+0.4*del2);
    vector<shared_ptr<ParamSurface> > sfs4 = 
      surf->subSurfaces(dom3.umin(), dom3.vmin()+0.4*del2, dom3.umax(), 
			dom3.vmax());
    RectDomain dom_3 = sfs3[0]->containingDomain();

    vector<Point> der5(6);
    sfs3[0]->point(der5, 0.5*(dom_3.umin()+dom_3.umax()),
		 0.5*(dom_3.vmin()+dom_3.vmax()), 2);
    std::cout << "Der5: " << der5[0] << std::endl
    	      << der5[1] << std::endl
    	      << der5[2] << std::endl;
    std::cout << der5[3] << std::endl
    	      << der5[4] << std::endl
    	      << der5[5] << std::endl;
    surf->writeStandardHeader(os);
    surf->write(os);
    sfs2[0]->writeStandardHeader(os);
    sfs2[0]->write(os);
    sfs3[0]->writeStandardHeader(os);
    sfs3[0]->write(os);
    sfs4[0]->writeStandardHeader(os);
    sfs4[0]->write(os);

    Point mid1, mid2;
    sfs3[0]->point(mid1, dom_3.umin(), 0.5*(dom_3.vmin()+dom_3.vmax()));
    sfs3[0]->point(mid2, dom_3.umax(), 0.5*(dom_3.vmin()+dom_3.vmax()));
    os << "400 1 0 4 255 0 0 255" << std::endl;
    os << "1" << std::endl;
    os << mid1 << std::endl;

    os << "400 1 0 4 255 0 0 255" << std::endl;
    os << "1" << std::endl;
    os << mid2 << std::endl;

    shared_ptr<SplineSurface> spline_sf2(sfs3[0]->asSplineSurface());
    if (spline_sf.get())
      {
	vector<Point> der4(6);
	spline_sf2->point(der4, 0.5*(dom_3.umin()+dom_3.umax()),
			  0.5*(dom_3.vmin()+dom_3.vmax()), 2);
	std::cout << "Der4: " << der4[0] << std::endl
		  << der4[1] << std::endl
		  << der4[2] << std::endl;
	std::cout << der4[3] << std::endl
		  << der4[4] << std::endl
		  << der4[5] << std::endl;
	spline_sf2->writeStandardHeader(os);
	spline_sf2->write(os);
      }
 
    vector<shared_ptr<ParamCurve> > cvs1 = 
      sfs3[0]->constParamCurves(dom_3.umin(), false);
    vector<shared_ptr<ParamCurve> > cvs2 = 
      sfs3[0]->constParamCurves(dom_3.umax(), false);
    vector<shared_ptr<ParamCurve> > cvs3 = 
      sfs3[0]->constParamCurves(dom_3.vmin(), true);
    vector<shared_ptr<ParamCurve> > cvs4 = 
      sfs3[0]->constParamCurves(dom_3.vmax(), true);
    vector<shared_ptr<ParamCurve> > cvs5 = 
      sfs3[0]->constParamCurves(0.2*dom_3.umin()+0.8*dom_3.umax(), false);
   vector<shared_ptr<ParamCurve> > cvs6 = 
      sfs3[0]->constParamCurves(0.2*dom_3.vmin()+0.8*dom_3.vmax(), true);

   if (cvs1.size() > 0)
     {
       cvs1[0]->writeStandardHeader(os);
       cvs1[0]->write(os);
     }
   if (cvs2.size() > 0)
     {
       cvs2[0]->writeStandardHeader(os);
       cvs2[0]->write(os);
     }
   if (cvs3.size() > 0)
     {
       cvs3[0]->writeStandardHeader(os);
       cvs3[0]->write(os);
     }
   if (cvs4.size() > 0)
     {
       cvs4[0]->writeStandardHeader(os);
       cvs4[0]->write(os);
     }
   if (spline_sf2.get())
     {
       double m1 = 0.5*(spline_sf2->startparam_u()+spline_sf2->endparam_u());
       double m2 = 0.5*(spline_sf2->startparam_v()+spline_sf2->endparam_v());
       Point p3 = spline_sf2->ParamSurface::point(m1, m2);
       std::cout << std::endl << "Par: (" << m1 << ", " << m2 << "), pos: " 
		 << p3 << std::endl;
       double upar, vpar, dist;
       Point clo_pt;
       surf->closestPoint(p3, upar, vpar, clo_pt, dist, eps);
       std::cout << "Par: (" << upar << ", "<< vpar << "), dist: " << dist 
		 << ", closest point: " << clo_pt << std::endl;
     }

    if (cvs5.size() > 0 && cvs6.size() > 0)
      {
	shared_ptr<ElementarySurface> elem =
	  dynamic_pointer_cast<ElementarySurface,ParamSurface>(sfs3[0]);
	shared_ptr<ElementaryCurve> elem_cv1 =
	  dynamic_pointer_cast<ElementaryCurve,ParamCurve>(cvs5[0]);
	shared_ptr<ElementaryCurve> elem_cv2 =
	  dynamic_pointer_cast<ElementaryCurve,ParamCurve>(cvs6[0]);
	shared_ptr<ParamCurve> par1 = 
	  elem->getElementaryParamCurve(elem_cv1.get(), eps);
	shared_ptr<ParamCurve> par2 = 
	  elem->getElementaryParamCurve(elem_cv2.get(), eps);
	if (par1.get())
	  {
	    shared_ptr<SplineCurve> scv1(CurveCreators::liftParameterCurve(par1, 
									   sfs3[0],
									   10.0*eps));
	    elem_cv1->writeStandardHeader(os);
	    elem_cv1->write(os);
	    scv1->writeStandardHeader(os);
	    scv1->write(os);
	  }
	if (par2.get())
	  {
	    shared_ptr<SplineCurve> scv2(CurveCreators::liftParameterCurve(par2, 
									   sfs3[0],
									   10.0*eps));
	    elem_cv2->writeStandardHeader(os);
	    elem_cv2->write(os);
	    scv2->writeStandardHeader(os);
	    scv2->write(os);
	  }
      }
    
    
     return 0;
}
