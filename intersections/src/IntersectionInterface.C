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

#include "GoTools/intersections/IntersectionInterface.h"
#include "GoTools/geometry/SplineCurve.h"
#include "GoTools/intersections/CvCvIntersector.h"
#include "GoTools/intersections/SplineCurveInt.h"
#include "GoTools/intersections/IntersectionPoint.h"
#include "GoTools/intersections/IntersectionCurve.h"
#include <fstream>

namespace Go
{

using namespace std;

//---------------------------------------------------------------------------
 void intersectCurves(shared_ptr<ParamCurve> crv1, shared_ptr<ParamCurve> crv2,
		      double tol, vector<pair<double, double> >& intersection_points)
//---------------------------------------------------------------------------
 {
     // Get spline curves
     SplineCurve *scurve1 = crv1->geometryCurve();
     SplineCurve *scurve2 = crv2->geometryCurve();
     if (!scurve1 || !scurve2)
     {
	 // Not spline curves. This can currently not be handled
	 return;
     }

     // To avoid shared pointer problems. We need a better solution
     shared_ptr<SplineCurve> curve1 = shared_ptr<SplineCurve>(scurve1->clone());
     shared_ptr<SplineCurve> curve2 = shared_ptr<SplineCurve>(scurve2->clone());

     // DEBUG. Draw curves
     std::ofstream out_file("int_crvs.g2");
      curve1->writeStandardHeader(out_file);
      curve1->write(out_file);
      curve2->writeStandardHeader(out_file);
      curve2->write(out_file);
     
    shared_ptr<ParamGeomInt> scurveint1 =
	shared_ptr<ParamGeomInt>(new SplineCurveInt (curve1));
    shared_ptr<ParamGeomInt> scurveint2 =
	shared_ptr<ParamGeomInt>(new SplineCurveInt (curve2));

    CvCvIntersector cvcvintersect (scurveint1, scurveint2, tol);
    cvcvintersect.compute();

    // Intersect
    vector<shared_ptr<IntersectionPoint> > intpts;
    vector<shared_ptr<IntersectionCurve> > intcrv;
    cvcvintersect.getResult(intpts, intcrv);

    // Report results
    size_t ki;
    for (ki=0; ki<intpts.size(); ki++)
	intersection_points.push_back(make_pair(intpts[ki]->getPar(0), intpts[ki]->getPar(1)));

    for (ki=0; ki<intcrv.size(); ki++)
    {
	// The middle point on the curve is reported
	int nmb_guide = intcrv[ki]->numGuidePoints();
	double par1 = 0.5*(intcrv[ki]->getGuidePoint(0)->getPar(0) +
			   intcrv[ki]->getGuidePoint(nmb_guide-1)->getPar(0));
	double par2 = 0.5*(intcrv[ki]->getGuidePoint(0)->getPar(1) +
			   intcrv[ki]->getGuidePoint(nmb_guide-1)->getPar(1));

	Point pnt1 = crv1->point(par1);
	Point pnt2;
	double dist, clo_par;
	crv2->closestPoint(pnt1, crv2->startparam(), crv2->endparam(), clo_par, pnt2, dist, &par2);
	
	intersection_points.push_back(make_pair(par1, clo_par));
    }

    
 }

} // namespace Go

