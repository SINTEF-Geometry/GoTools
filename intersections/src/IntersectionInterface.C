//===========================================================================
//                                                                           
// File: IntersectionInterface.C
//                                                                           
// Created: Dec. 08
//                                                                           
// Author: Vibeke Skytt
//                                                                           
// Revision: 
//                                                                           
//===========================================================================

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

