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

//#include <iostream>
#include "GoTools/compositemodel/ftCurve.h"
#include "GoTools/compositemodel/ftPlane.h"
#include "GoTools/geometry/CurveOnSurface.h"
#include "GoTools/geometry/BoundedSurface.h"

#include "sislP.h"
#include "GoTools/geometry/SISLconversion.h"
#include <vector>

using std::vector;

namespace Go
{

namespace // anon namespace
{

//----------------------- Helper function for ChopOff -----------------------
vector<double> intersect_curve_plane(const SplineCurve& cv,
				   const ftPlane& plane)
//---------------------------------------------------------------------------
{
    SISLCurve* sislcv = Curve2SISL(cv, false);
    double* pt = const_cast<double*>(plane.point().begin());
    double* norm = const_cast<double*>(plane.normal().begin());
    int dim = 3;
    double epsco = 1e-15; // Not used
    double epsge = 1e-6;
    int numintpt;
    double* intpar = 0;
    int numintcu;
    SISLIntcurve** intcurves = 0;
    int stat;
    // Find the topology of the intersection
    s1850(sislcv, pt, norm, dim, epsco, epsge,
	  &numintpt, &intpar, &numintcu, &intcurves, &stat);
    vector<double> pars(intpar, intpar + numintpt);
    for (int i = 0; i < numintcu; ++i) {
	int n = intcurves[i]->ipoint;
	pars.push_back(intcurves[i]->epar1[0]);
	pars.push_back(intcurves[i]->epar1[n-1]);
    }
    free(intpar);
    freeIntcrvlist(intcurves, numintcu);
    freeCurve(sislcv);
    return pars;
}

} // anon namespace

//===========================================================================
vector<ftCurveSegment> ftCurveSegment::chopOff(const BoundingBox& box,
					       bool& eraseme)
//===========================================================================
{
    SplineCurve* spacespline
	= dynamic_cast<SplineCurve*>(space_curve_.get());
    if (spacespline == 0) {
	CurveOnSurface* sc
	    = dynamic_cast<CurveOnSurface*>(space_curve_.get());
	if (sc) {
	    spacespline = dynamic_cast<SplineCurve*>(sc->spaceCurve().get());
	} else {
	    THROW("ChopOff() can only be called on splinecurve segments.");
	}
    }
    ASSERT(spacespline != 0);

    int i, j;
    vector<ftCurveSegment> nseg;
    eraseme = false;

    // Check if the bounding box of space_curve is contained in or
    // does not overlap the box

    const BoundingBox& this_box = space_curve_->boundingBox();
    if (box.containsBox(this_box))
	return nseg;
    if (!box.overlaps(this_box)) {
	// We should destroy this segment!
	eraseme = true;
	return nseg;
    }

    // Intersect space_curve_ with every plane that bounds the box
    vector<double> intpars;
    vector<double> newpars;
    // We start off with the start parameter...
    intpars.push_back(space_curve_->startparam());
    for (i = 0; i < 2; ++i) {
	Point p = (i==0) ? box.low() : box.high();
	for (j = 0; j < 3; ++j) {
	    Point n(0.0, 0.0, 0.0);
	    n[j] = 1.0;
	    ftPlane pl(n, p);
	    if (pl.intersectsBox(space_curve_->boundingBox())) {
		newpars = intersect_curve_plane(*spacespline, pl);
		intpars.insert(intpars.end(), newpars.begin(), newpars.end());
	    }
	}
    }
    // ... and end with the end parameter
    intpars.push_back(space_curve_->endparam());

    // Sort intpars and remove any duplicate members
    std::sort(intpars.begin(), intpars.end());
    intpars.erase(std::unique(intpars.begin(), intpars.end()), intpars.end());


    // Now we go through each successive pair of intpars.
    // We check every subsegment to see if it is inside the box,
    // and split the main curve and the parameter curves into subcurves.

    int n = (int)intpars.size();
    bool any_inside = false;
    vector<bool> interval_inside(n - 1, false);
    Point halfway_point(space_curve_->dimension());
    for (i = 0; i < n - 1; ++i) {
	double halfway = 0.5*intpars[i] + 0.5*intpars[i+1];
	space_curve_->point(halfway_point, halfway);
	if (box.containsPoint(halfway_point)) {
	    interval_inside[i] = true;
	    any_inside = true;
	}
    }

    if (!any_inside) {
	eraseme = true;
	return nseg;
    }

    // We will need the SISLCurve version of the spline curves
    ParamCurve* cvs[3] = {0, 0, 0};
    ParamCurve* newcvs[3] = {0, 0, 0};
    const int numcv = 3;
    cvs[0] = parameter_curve_[0].get();
    cvs[1] = parameter_curve_[1].get();
    cvs[2] = space_curve_.get();

    bool first_pick = true;
    for (i = 0; i < n - 1; ++i) {
	if (interval_inside[i]) {
	    // Find the next interval that is outside, if any.
	    // The while loop skips over tangential cases
	    int kk = i+1;
	    while (kk < n - 1 && interval_inside[kk]) ++kk;
	    // We've got an interval inside with the param range
	    // [intpars[i], intpars[kk]]
	    if (i == 0 && kk == n - 1) {
		// The whole curve is inside! We skip out of this.
		break;
	    }
	    // We have to pick curve parts.
	    for (j = 0; j < numcv; ++j) {
		if (cvs[j] != 0)
		    newcvs[j] = cvs[j]->subCurve(intpars[i], intpars[kk]);
	    }

	    ftCurveSegment newseg(segment_type_, JOINT_DISC,
				  underlying_face_[0], underlying_face_[1],
				  shared_ptr<ParamCurve>(newcvs[0]),
				  shared_ptr<ParamCurve>(newcvs[1]),
				  shared_ptr<ParamCurve>(newcvs[2]));
	    if (first_pick) {
		// Modify this ftCurveSegment.
		first_pick = false;
		(*this) = newseg;
	    } else {
		// Add a new element to nseg.
		nseg.push_back(newseg);
	    }
	    // Let i jump to next possible interval inside
	    i = kk + 1;
	}
    }

    return nseg;
}

//===========================================================================
void ftCurveSegment::redefineSpaceCurve(double eps_go)
//===========================================================================
{
//     if (space_curve_.get() != 0 && 
// 	space_curve_->instanceType() != Class_SplineCurve) {
// 	ALWAYS_ERROR_IF(space_curve_->instanceType != Class_CurveOnSurface,
// 		    "Unrecognized curve type encountered in redefineSpaceCurve().");
// 	// if we got here, we have a spacecurve that is a CurveOnSurface
// 	ALWAYS_ERROR_IF(underlying_face_[0] != 0 || underlying_face_[1] != 0 ||
// 		    parameter_curve_[0] != 0 || parameter_curve_[1] != 0,
// 		    "Possibly conflicting information between this segment's space curve "
// 		    "(which is detected to be of type CurveOnSurface) and this segment's "
// 		    "underlying face and parameter curve. Aborting execution of "
// 		    "redefineSpaceCurve().", 
// 		    InputError());
// 	// We know that the CurveOnSurface object is the only information present.

// 	if (space_curve_->parPref() == false && 
// 	    space_curve_->spaceCurve()->instanceType() == Class_SplineCurve) {
// 	    // The geometric information present in the CurveOnSurface object is kept, the
// 	    // rest is thrown away
// 	    space_curve_ = space_curve_->spaceCurve();
// 	} else {
// 	    // We will 'dissect' the CurveOnSurface, and put its information into
// 	    // underlying_face_[0] and parameter_curve_[0], after which we will reconstruct
// 	    // a spacecurve of type SplineCurve
	
// 	    underlying_face_[0] = new ftSurface(space_curve_->underlyingSurface(), -1);
// 	    parameter_curve_[0] = space_curve_->parameterCurve();
// 	    space_curve_.reset(); // don't worry!  It will be "reborn" a little later...
// 	}
//     }

    if (space_curve_.get() == 0) {
	// This object does not have defined any spacecurve.  
	// Let us define it!
	shared_ptr<ParamSurface> ul_face;
	shared_ptr<ParamCurve> ul_curve;

	if (underlying_face_[0] && parameter_curve_[0]) {
	    ul_face = underlying_face_[0]->surface();
	    ul_curve = parameter_curve_[0];
	} else if (underlying_face_[1] && parameter_curve_[1]) {
	    ul_face = underlying_face_[1]->surface();
	    ul_curve = parameter_curve_[1];
	} 
	ALWAYS_ERROR_IF(ul_face.get() == 0 || ul_curve.get() == 0,
		    "redefineSpaceCurve() was not able to find any curve information in "
			"this object.");

	SISLSurf* ul_spline_surf = 0;
	SISLCurve* param_spline_curve = 0;
	
	shared_ptr<SplineSurface> tmp_spline;
	SplineSurface* splinesf = ul_face->getSplineSurface();
	if (!splinesf)
	  {
	    // Convert to spline surface
	    tmp_spline = shared_ptr<SplineSurface>(ul_face->asSplineSurface());
	    splinesf = tmp_spline.get();
	  }
	ul_spline_surf = GoSurf2SISL(*splinesf, false);
	
	// if (ul_face->instanceType() == Class_SplineSurface) {
	//     ul_spline_surf = 
	// 	GoSurf2SISL(*(dynamic_cast<SplineSurface*>(ul_face.get())), false);
	// } else if (ul_face->instanceType() == Class_BoundedSurface) {
	//     BoundedSurface* bsurf = dynamic_cast<BoundedSurface*>(ul_face.get());
	//     shared_ptr<ParamSurface> psurf = bsurf->underlyingSurface();
	//     if (psurf->instanceType() == Class_SplineSurface) {
	// 	ul_spline_surf = 
	// 	    GoSurf2SISL(*(dynamic_cast<SplineSurface*>(psurf.get())), false);
	//     }
	// }
	if (ul_curve->instanceType() == Class_SplineCurve) {
	    param_spline_curve = 
		Curve2SISL(*(dynamic_cast<SplineCurve*>(ul_curve.get())), false); 
	}
	
	if (ul_spline_surf == 0 || param_spline_curve == 0)
	  THROW("redefineSpaceCurve() was not able to convert the curve information into any object it could use.");


	// making a spacecurve in SISL format
	SISLCurve* dummycurve1 = 0;
	SISLCurve* dummycurve2 = 0;
	SISLCurve* resultcurve = 0;
	int stat;
	s1383(ul_spline_surf, // underlying surface
	      param_spline_curve, // curve in parameter domain of underlying surface
	      eps_go, // tolerance
	      0, // maximum step length (0 -> calculated from bounding box)
	      0, // calculate no derivatives
	      &resultcurve,
	      &dummycurve1,
	      &dummycurve2,
	      &stat);
	ALWAYS_ERROR_IF(stat < 0, 
			"Error occured during execution of SISL routine s1383");

	// making a space curve in Go format
	space_curve_ = shared_ptr<ParamCurve>(SISLCurve2Go(resultcurve));

	// freeing memory used for temporary SISL objects
	if (resultcurve) freeCurve(resultcurve);
	if (dummycurve1) freeCurve(dummycurve1);
	if (dummycurve2) freeCurve(dummycurve2);
	if (param_spline_curve) freeCurve(param_spline_curve);
	if (ul_spline_surf) freeSurf(ul_spline_surf);
    }
}

//===========================================================================
void ftCurveSegment::reparametrize(double eps_go)
//===========================================================================
{
    // We will reparametrize this segment from 0 to its length.
    // Our algorithm is as follows:
    //    o Convert the segment to a series of Bezier curves,
    //      with only order-tuple knots.
    //    o Compute the arclengths of these Bezier-curves
    //    o Make the difference between the order-tuple sets of
    //      knots equal to the arclengths

    // checking whether spacecurve exist, generate it if it does not.

    SplineCurve* spacespline = dynamic_cast<SplineCurve*>(space_curve_.get());
    if (!spacespline) {
	// the space_curve_ was either 0 or not a SplineCurve.  We must transform it
	redefineSpaceCurve(eps_go);
	spacespline = dynamic_cast<SplineCurve*>(space_curve_.get());
	ALWAYS_ERROR_IF(!spacespline,
			"reparametrize() was not able to find or generate a proper "
			"spatial representation of the curveSegment in question.");

    }

    int i;

    SISLCurve* sislcv = Curve2SISL(*spacespline);
    SISLCurve* newcv;
    int stat;
    s1730(sislcv, &newcv, &stat);
    int ord = newcv->ik;
    int numbez = newcv->in / ord;
    double startpar, endpar;
    vector<double> coef(ord * 4);
    vector<double> arclen(numbez);
    vector<double> knots(ord*2);
    std::fill(knots.begin(), knots.begin() + ord, 0.0);
    std::fill(knots.begin() + ord, knots.end(), 1.0);

    for (i = 0; i < numbez; ++i) {
	s1732(newcv, i, &startpar, &endpar, &coef[0], &stat);
	SISLCurve* bez = newCurve(ord, ord, &knots[0], &coef[0], 1, 3, 0);
	double eps = 0.0001;
	s1240(bez, eps, &arclen[i], &stat);
	freeCurve(bez);
    }
    double par = 0.0;
    std::fill(newcv->et, newcv->et + ord, par);
    for (i = 0; i < numbez; ++i) {
	par += arclen[i];
	std::fill(newcv->et + ord*(i+1), newcv->et + ord*(i+2), par);
    }
    
    space_curve_ = shared_ptr<SplineCurve>(new SplineCurve(newcv->in,
							       newcv->ik,
							       newcv->et,
							       newcv->ecoef,
							       3));
    freeCurve(newcv);
    freeCurve(sislcv);

    parameter_curve_[0].reset();
    parameter_curve_[1].reset();
}

} // namespace Go
