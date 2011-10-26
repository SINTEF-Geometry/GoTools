//===========================================================================
//                                                                           
// File: Singular.C
//                                                                           
// Created: April 2008
//                                                                           
// Author: Vibeke Skytt
//                                                                           
// Revision:
//                                                                           
// Description:
//                                                                           
//===========================================================================

#include "GoTools/intersections/Singular.h"
#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/geometry/BoundedSurface.h"
#include "GoTools/intersections/SfPtIntersector.h"
#include "GoTools/intersections/CvPtIntersector.h"
#include "GoTools/intersections/SplineSurfaceInt.h"
#include "GoTools/intersections/SplineCurveInt.h"
#include "GoTools/intersections/ParamPointInt.h"
#include "GoTools/intersections/GeoTol.h"
#include "GoTools/intersections/IntersectionPoint.h"
#include "GoTools/intersections/IntersectionCurve.h"
#include "GoTools/utils/Array.h"
#include "GoTools/geometry/CurveBoundedDomain.h"

using std::vector;
using std::shared_ptr;
using std::dynamic_pointer_cast;

namespace Go {

/// Compute vanishing curve tangents and surface normals


    void Singular::vanishingNormal(std::shared_ptr<ParamSurface> srf, double tol,
				   vector<Point>& singular_pts,  // Singular points in the parameter domain 
				   vector<vector<Point> >& singular_sequences)  // Sequences of parameter points
                                                                        	// making a singular curve
	{
	    double partol = 1.0e-6;  // For the time being

	    // Get the underlying surface of a bounded surface
	    shared_ptr<ParamSurface> curr_srf = srf;
	    while (curr_srf->instanceType() == Class_BoundedSurface)
	    {
		shared_ptr<BoundedSurface> bd_srf = dynamic_pointer_cast<BoundedSurface,ParamSurface>(curr_srf);
		curr_srf = bd_srf->underlyingSurface();
	    }

	    // For the time being this should not be any problem
	    if (curr_srf->instanceType() != Class_SplineSurface)
		return;

	    vector<shared_ptr<IntersectionPoint> > intpts;
	    vector<shared_ptr<IntersectionCurve> > intcrv;

	    shared_ptr<SplineSurface> splinesf = 
		dynamic_pointer_cast<SplineSurface,ParamSurface>(curr_srf);

	    // Check if the surface is degenerate. In that case it is known that there
	    // are vanishing normals, and it takes a long time to compute it
	    double umin = splinesf->startparam_u();
	    double umax = splinesf->endparam_u();
	    double vmin = splinesf->startparam_v();
	    double vmax = splinesf->endparam_v();
	    bool b, r, t, l;
	    bool is_degenerate = splinesf->isDegenerate(b, r, t, l, tol);
	    if (is_degenerate)
	      {
		// Pick a part of the surface to avoid the degeneracy
		double fac = 1.0e-4;
		double del_u = umax - umin;
		double del_v = vmax - vmin;
		double umin2 = l ? umin + fac*del_u : umin;
		double umax2 = r ? umax - fac*del_u : umax;
		double vmin2 = b ? vmin + fac*del_v : vmin;
		double vmax2 = t ? vmax - fac*del_v : vmax;
		if (b)
		  vmin += fac*del_v;
		if (t)
		  vmax -= fac*del_v;
		if (l)
		  umin += fac*del_u;
		if (r)
		  umax += fac*del_u;

		shared_ptr<SplineSurface> sub_sf = 
		  shared_ptr<SplineSurface>(splinesf->subSurface(umin2, vmin2, 
								 umax2, vmax2));
		splinesf = sub_sf;
	      }

	    // Make normal surface
	    shared_ptr<SplineSurface> normalsf = shared_ptr<SplineSurface>(splinesf->normalSurface());
	    shared_ptr<ParamSurfaceInt> nsurf = shared_ptr<ParamSurfaceInt>(new SplineSurfaceInt(normalsf));

	    // Define origo
	    int dim = srf->dimension();
	    shared_ptr<Point> origo(new Point(dim));
	    origo->setValue(0.0);
	    shared_ptr<ParamPointInt> origo_int(new ParamPointInt(origo));
		
	    // Intersect
	    shared_ptr<GeoTol> eps = shared_ptr<GeoTol>(new GeoTol(tol));
	    shared_ptr<SfPtIntersector>
		sfptint(new SfPtIntersector(nsurf, origo_int, eps));
	    sfptint->setSelfintCase(1);

	    sfptint->compute();
							
	    sfptint->getResult(intpts, intcrv);

	    size_t ki;
	    int kj;
	    if (srf->instanceType() == Class_BoundedSurface)
	    {
		// Trim self intersection results with the surface domain
		shared_ptr<BoundedSurface> bd_srf = dynamic_pointer_cast<BoundedSurface,ParamSurface>(srf);
		const CurveBoundedDomain& domain = bd_srf->parameterDomain();
		for (ki=0; ki<intpts.size(); ki++)
		{
		    Array<double,2> point(intpts[ki]->getPar(0),intpts[ki]->getPar(1));
		    if (domain.isInDomain(point, partol))
			singular_pts.push_back(intpts[ki]->getPar1Point());
		}

		for (ki=0; ki<intcrv.size(); ki++)
		{
		    vector<Point> current_seq;
		    int num_pts = intcrv[ki]->numGuidePoints();
		    for (kj=1; kj<num_pts; kj++)
		    {
			Point p1 = intcrv[ki]->getGuidePoint(kj-1)->getPar1Point();
			Point p2 = intcrv[ki]->getGuidePoint(kj)->getPar1Point();
			SplineCurve linear_seg(p1, 0.0, p2, 1.0);
			vector<double> params;
			domain.findPcurveInsideSegments(linear_seg, partol, params);
			for (size_t kr=1; kr<params.size(); kr+=2)
			{
			    Point curr = (1.0 - params[kr-1])*p1 + params[kr-1]*p2;
			    current_seq.push_back(curr);
			    curr = (1.0 - params[kr])*p1 + params[kr]*p2;
			    current_seq.push_back(curr);
			    if (1.0 - params[kr] > partol)
			    {
				singular_sequences.push_back(current_seq);
				current_seq.clear();
			    }
			}
				
		    }
		    
		}
		// Add information due to degeneracy
		vector<Point> deg_seq;
		if (b)
		  {
		    deg_seq.push_back(Point(umin,vmin));
		    deg_seq.push_back(Point(umax,vmin));
		  }
		if (t)
		  {
		    deg_seq.push_back(Point(umin,vmax));
		    deg_seq.push_back(Point(umax,vmax));
		  }
		if (l)
		  {
		    deg_seq.push_back(Point(umin,vmin));
		    deg_seq.push_back(Point(umin,vmax));
		  }
		if (r)
		  {
		    deg_seq.push_back(Point(umax,vmin));
		    deg_seq.push_back(Point(umax,vmax));
		  }
		for (size_t kh=0; kh<deg_seq.size(); kh+=2)
		  {
		    SplineCurve linear_seg(deg_seq[kh], 0.0, deg_seq[kh+1], 1.0);
		    vector<double> params;
		    domain.findPcurveInsideSegments(linear_seg, partol, params);
		    vector<Point> current_seq;
		    for (size_t kr=1; kr<params.size(); kr+=2)
		      {
			Point curr = (1.0 - params[kr-1])*deg_seq[kh] + 
			  params[kr-1]*deg_seq[kh+1];
			current_seq.push_back(curr);
			curr = (1.0 - params[kr])*deg_seq[kh] + params[kr]*deg_seq[kh+1];
			current_seq.push_back(curr);
			if (1.0 - params[kr] > partol)
			  {
			    singular_sequences.push_back(current_seq);
			    current_seq.clear();
			  }
		      }
		  }
	    }
	    else if (srf->instanceType() == Class_SplineSurface)
	    {
		for (ki=0; ki<intpts.size(); ki++)
		    singular_pts.push_back(intpts[ki]->getPar1Point());

		for (ki=0; ki<intcrv.size(); ki++)
		{
		    vector<Point> current_seq;
		    int num_pts = intcrv[ki]->numGuidePoints();
		    for (kj=0; kj<num_pts; kj++)
			current_seq.push_back(intcrv[ki]->getGuidePoint(kj)->getPar1Point());
		    singular_sequences.push_back(current_seq);
		}

		// Add information due to degeneracy
		if (b)
		  {
		    vector<Point> current_seq;
		    current_seq.push_back(Point(umin,vmin));
		    current_seq.push_back(Point(umax,vmin));
		    singular_sequences.push_back(current_seq);
		  }
		if (t)
		  {
		    vector<Point> current_seq;
		    current_seq.push_back(Point(umin,vmax));
		    current_seq.push_back(Point(umax,vmax));
		    singular_sequences.push_back(current_seq);
		  }
		if (l)
		  {
		    vector<Point> current_seq;
		    current_seq.push_back(Point(umin,vmin));
		    current_seq.push_back(Point(umin,vmax));
		    singular_sequences.push_back(current_seq);
		  }
		if (r)
		  {
		    vector<Point> current_seq;
		    current_seq.push_back(Point(umax,vmin));
		    current_seq.push_back(Point(umax,vmax));
		    singular_sequences.push_back(current_seq);
		  }
	    }
	}


    void Singular::vanishingTangent(shared_ptr<ParamCurve> crv, 
				    double start, double end, double tol,
				    vector<double>& singular_pts,  // Singular points in the parameter domain 
				    vector<vector<double> >& singular_sequences)  // Sequences of parameter points
                                                                        	// making a singular curve
	{
	    double partol = 1.0e-6;  // For the time being

	    // Pick the relevant part of the curve
	    shared_ptr<ParamCurve> curr_crv;
	    if (start > crv->startparam()+partol || end < crv->endparam()-partol)
		curr_crv = shared_ptr<ParamCurve>(crv->subCurve(start, end));
	    else
		curr_crv = crv;

	    // Get spline curve
	    SplineCurve *spline = curr_crv->geometryCurve();
	    if (spline == 0)
		return;   // Elementary curve. No vanising tangent (I hope, nov 08)

	    // Compute derivative curve
	    shared_ptr<SplineCurve> deriv_crv = shared_ptr<SplineCurve>(spline->derivCurve(1));

	    // Intersect with origo
	    // Define origo
	    int dim = crv->dimension();
	    shared_ptr<Point> origo(new Point(dim));
	    origo->setValue(0.0);
	    shared_ptr<ParamPointInt> origo_int(new ParamPointInt(origo));
	    
	    // Intersect
	    vector<shared_ptr<IntersectionPoint> > intpts;
	    vector<shared_ptr<IntersectionCurve> > intcrv;
	    shared_ptr<ParamCurveInt> crv_int = shared_ptr<ParamCurveInt>(new SplineCurveInt(deriv_crv));	    
	    shared_ptr<GeoTol> eps = shared_ptr<GeoTol>(new GeoTol(tol));
	    shared_ptr<CvPtIntersector>
		cvptint(new CvPtIntersector(crv_int, origo_int, eps));
	    cvptint->setSelfintCase(1);

	    cvptint->compute();
							
	    cvptint->getResult(intpts, intcrv);

	    size_t ki, kj;
	    for (ki=0; ki<intpts.size(); ki++)
		singular_pts.push_back(intpts[ki]->getPar1Point()[0]);

	    for (ki=0; ki<intcrv.size(); ki++)
	    {
		vector<double> current_seq;
		int num_pts = intcrv[ki]->numGuidePoints();
		for (kj = 0; (int)kj < num_pts; kj++)
		    current_seq.push_back(intcrv[ki]->getGuidePoint((int)kj)->getPar1Point()[0]);
		singular_sequences.push_back(current_seq);
	    }
	}


} // end namespace Go

