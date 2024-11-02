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

#include "GoTools/geometry/GeometryTools.h"
#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/geometry/BoundedSurface.h"
#include "GoTools/geometry/SplineCurve.h"
#include "GoTools/geometry/CurveOnSurface.h"
#include "GoTools/geometry/SplineUtils.h"
#include "GoTools/geometry/Utils.h"
#include "GoTools/geometry/ElementaryUtils.h"
#include "GoTools/geometry/ElementarySurface.h"
// #include "values.h"

#include <fstream>


namespace Go
{

using namespace std;

//---------------------------------------------------------------------------
int GeometryTools::analyzePeriodicityDerivs(const ParamCurve& cv, int max_derivs, double tol)
//---------------------------------------------------------------------------
{
    // Compute derivatives at start and end.
    std::vector<Point> ptsbeg(max_derivs + 1);
    std::vector<Point> ptsend(max_derivs + 1);
    cv.point(ptsbeg, cv.startparam(), max_derivs);
    cv.point(ptsend, cv.endparam(), max_derivs);

    // Compare computed points
    int i = 0;
    for (; i <= max_derivs; ++i) {
	double diff = ptsbeg[i].distInf(ptsend[i]);
	if (diff > tol) break;
    }
    return i - 1;
}

//---------------------------------------------------------------------------
int GeometryTools::analyzePeriodicity(const SplineCurve& cv, double knot_tol)
//---------------------------------------------------------------------------
{
    // Analyze the knot vectors for periodicity.
    // That is, check if the first knots are repeated
    // (translated) at the end.
    int basper = analyzePeriodicity(cv.basis(), knot_tol);
    if (basper == -1) return -1;

    // Analyze the coefficient vector for periodicity.
    // The (basper+1) first must equal the (basper+1) last.
    int n = cv.numCoefs() - basper - 1;
    int dim = cv.dimension();
    std::vector<double>::const_iterator c = cv.coefs_begin();
    if (cv.rational()) {
	++dim;
	c = cv.rcoefs_begin();
    }
    for (int i = 0; i < basper + 1; ++i) {
	for (int dd = 0; dd < dim; ++dd) {
	    if (c[dim*i + dd] != c[dim*(i+n) + dd]) {
		return -1;
	    }
	}
    }
    return basper;
}


//---------------------------------------------------------------------------
int GeometryTools::analyzePeriodicityDerivs(const SplineSurface& sf,
			     int direction,
			     int max_derivs,
			     double tol)
//---------------------------------------------------------------------------
{
    if (direction < 0 || direction > 1) {
	THROW("Error in direction parameter. Must be 0 or 1.");
    }
    // Make a hypercurve from this surface.
    // Direction parameter is one more for representSurface...() :-P
    shared_ptr<SplineCurve> cv
	= GeometryTools::representSurfaceAsCurve(sf, direction + 1);
    return GeometryTools::analyzePeriodicityDerivs(*cv, max_derivs, tol);
}



//---------------------------------------------------------------------------
int GeometryTools::analyzePeriodicity(const SplineSurface& sf, int direction, double knot_tol)
//---------------------------------------------------------------------------
{
    if (direction < 0 || direction > 1) {
	THROW("Error in direction parameter. Must be 0 or 1.");
    }
    // Make a hypercurve from this surface.
    // Direction parameter is one more for representSurface...() :-P
    shared_ptr<SplineCurve> cv
	= representSurfaceAsCurve(sf, direction + 1);
    return analyzePeriodicity(*cv, knot_tol);
}

//---------------------------------------------------------------------------
int GeometryTools::analyzePeriodicity(const BsplineBasis& basis, double knot_tol)
//---------------------------------------------------------------------------
{
    // Check if the first knots are repeated
    // (translated) at the end.
    // We must first find our 'target' number of repeats,
    // by checking the multiplicity of the startparam() knot.
    int emul = basis.endMultiplicity(true);
    if (basis.endMultiplicity(false) != emul) {
	return -1;
    }
    int k = basis.order();
    int n = basis.numCoefs();
    double delta = basis.endparam() -  basis.startparam();
    std::vector<double>::const_iterator kn = basis.begin();
    int i = 1;
    for (; i < k-emul+1; ++i) {
	double d1 = kn[n+emul-1+i] - kn[k-1+i];
	double d2 = kn[n-i] - kn[k-emul-i];
	// We only allow minimal floating point errors
	if (std::max(std::abs(d1-delta), std::abs(d2-delta)) > knot_tol) {
	    break;
	}
    }
    return i - 2;
}

//---------------------------------------------------------------------------
  bool GeometryTools::commonSeam(shared_ptr<CurveOnSurface> bd_cv1,
				 shared_ptr<CurveOnSurface> bd_cv2,
				 double tol, double angtol,
				 int& pardir, double& parval1, 
				 double& parval2)
//---------------------------------------------------------------------------
{
  pardir = -1;
  parval1 = parval2 = 0.0;

  shared_ptr<ParamSurface> surf1 = bd_cv1->underlyingSurface();
  shared_ptr<ParamSurface> surf2 = bd_cv2->underlyingSurface();
  ElementarySurface *elem1 = surf1->elementarySurface();
  ElementarySurface *elem2 = surf2->elementarySurface();
  if (elem1 == NULL || elem2 == NULL)
    return false;  // Not elementary surfaces

  if (!ElementaryUtils::sameElementarySurface(surf1.get(), surf2.get(),
					      tol, angtol))
    return false;

  int dir1, dir2;
  double val1, val2;
  bool constpar1 = bd_cv1->isConstantCurve(tol, dir1, val1);
  bool constpar2 = bd_cv2->isConstantCurve(tol, dir2, val2);

  if (constpar1 == false || constpar2 == false)
    return false;

  if (dir1 != dir2)
    return false;

  pardir = dir1;
  parval1 = val1;
  parval2 = val2;

  bool seam1 = elem1->atSeam(dir1, val1);
  bool seam2 = elem2->atSeam(dir2, val2);
  bool period = elem1->fullPeriod(dir1, val1, val2);

  if ((seam1 && seam2) || period)
    return true;

  return false;
}

//---------------------------------------------------------------------------
// Note that the surface size is computed with respect to the containing
  // domain
void  GeometryTools::estimateSurfaceSize(const ParamSurface& srf, double& length_u, 
			  double& length_v, double* area)
//---------------------------------------------------------------------------
{
  const int ncvsample = 3;
  const int nptsample = 5;
  length_u = length_v = 0.0;

  int ki, kj;
  double umin, vmin, umax, vmax;
  if (area)
    {
      umin = area[0];
      umax = area[1];
      vmin = area[2];
      vmax = area[3];
    }
  else
    {
      RectDomain domain = srf.containingDomain();
      umin = domain.umin();
      vmin = domain.vmin();
      umax = domain.umax();
      vmax = domain.vmax();
    }
  double uint = (umax - umin)/(double)(ncvsample-1);
  double vint = (vmax - vmin)/(double)(nptsample-1);

  double par_u, par_v;
  Point pt1, pt2;
  for (ki=0, par_u=umin; ki<ncvsample; ki++, par_u+=uint)
    {
      par_v = vmin;
      pt1 = srf.point(par_u, par_v);
      for (kj=1; kj<=nptsample; kj++, par_v+=vint)
	{
	  pt2 = srf.point(par_u, par_v);
	  length_v += pt1.dist(pt2);
	  pt1 = pt2;
	}
    }
  length_v /= (double)ncvsample;
  
  uint = (umax - umin)/(double)(nptsample-1);
  vint = (vmax - vmin)/(double)(ncvsample-1);
  for (ki=0, par_v=vmin; ki<ncvsample; ki++, par_v+=vint)
    {
      par_u = umin;
      pt1 = srf.point(par_u, par_v);
      for (kj=1; kj<=nptsample; kj++, par_u+=uint)
	{
	  pt2 = srf.point(par_u, par_v);
	  length_u += pt1.dist(pt2);
	  pt1 = pt2;
	}
    }
  length_u /= (double)ncvsample;
}

//-------------------------------------------------------------------
void GeometryTools::estimateIsoCurveLength(const ParamSurface& srf, bool dir_u, 
			    double par, double& length)
//-------------------------------------------------------------------
{
  int nptsample = 5;
  length = 0;
  int ki;
  double sfpar[2];
  RectDomain dom = srf.containingDomain();
  double ta = (dir_u) ? dom.umin() : dom.vmin();
  double tb = (dir_u) ? dom.umax() : dom.vmax();
  double tint = (tb - ta)/(double)(nptsample-1);
  int idxc = (dir_u) ? 1 : 0; // Index of constant parameter, i.e. par
  int idxv = 1 - idxc;
  sfpar[idxc] = par;
  sfpar[idxv] = ta;
  Point pt1, pt2;
  pt1 = srf.ParamSurface::point(sfpar[0],sfpar[1]);
  for (ki=1, sfpar[idxv]+=tint; ki<nptsample; ki++, sfpar[idxv]+=tint)
    {
      pt2 = srf.ParamSurface::point(sfpar[0],sfpar[1]);
      length += pt1.dist(pt2);
      pt1 = pt2;
    }
}

//-------------------------------------------------------------------
bool GeometryTools::degenerateToCurve(const SplineSurface& srf, bool dir_u, double tol)
//-------------------------------------------------------------------
{
    // For each row of coefficients in the appropriate parameter
    // direction, make a rough over estimate of the curve length by
    // computing the length of the control polygon. Check if this length
    // exceeds the given tolerance.
    int nmb1 = (dir_u) ? srf.numCoefs_u() : srf.numCoefs_v();
    int nmb2 = (dir_u) ? srf.numCoefs_v() : srf.numCoefs_u();
    int dim = srf.dimension();
    vector<double>::const_iterator coefs = srf.coefs_begin();
    int ki, kj;
    int kr1 = (dir_u) ? dim : nmb2*dim;
    int kr2 = (dir_u) ? nmb1*dim : dim;
    for (kj=0; kj<nmb2; kj++)
    {
	vector<double>::const_iterator c1 = coefs + kr2;
	double tlen = 0.0;
	for (ki=1; ki<nmb1; ki++, c1+=kr1)
	{
	    double d1 = Utils::distance_squared(&c1[0], &c1[0]+dim, &c1[0]+dim);
	    tlen += d1;
	}
	if (tlen > tol)
	    break;
    }

    return (kj < nmb2) ? false : true;
}

//-------------------------------------------------------------------
void GeometryTools::makeBdDegenerate(SplineSurface& srf, int bd_idx)  // left, right, bottom, top
//-------------------------------------------------------------------
{
    bool dir_u = (bd_idx == 2 || bd_idx == 3);
    int nmb1 = (dir_u) ? srf.numCoefs_u() : srf.numCoefs_v();
    int nmb2 = (dir_u) ? srf.numCoefs_v() : srf.numCoefs_u();
    int dim = srf.dimension();
    vector<double>::iterator coefs = srf.coefs_begin();
    int kr1 = (dir_u) ? dim : nmb1*dim;

    // First compute mean coefficient
    vector<double>::iterator c1 = coefs;
    if (bd_idx == 3)
	c1 += (nmb2-1)*nmb1*dim;
    if (bd_idx == 1)
	c1 += (nmb2-1)*dim;

    vector<double> deg_pt(dim, 0.0);
    for (int ki=0; ki<nmb1; ki++, c1+=kr1)
    {
	vector<double>::iterator d1 = c1;
	for (int kj=0; kj<dim; kj++, d1++)
	    deg_pt[kj] += c1[kj];
    }
    for (int kj=0; kj<dim; kj++)
	deg_pt[kj] /= (double)nmb1;

    c1 = coefs;
    if (bd_idx == 3)
	c1 += (nmb2-1)*nmb1*dim;
    if (bd_idx == 1)
	c1 += (nmb2-1)*dim;

    for (int ki=0; ki<nmb1; ki++, c1+=kr1)
    {
	vector<double>::iterator d1 = c1;
	for (int kj=0; kj<dim; kj++, d1++)
	    c1[kj] = deg_pt[kj];
    }
}

//-------------------------------------------------------------------
 bool GeometryTools::checkConstantCoef(SplineCurve& cv, int idx, double val, 
			double max_dist, double tol)

// Check if a curve coefficient is equal to a constant in a specified dimension
// provided it already lies close
//-------------------------------------------------------------------
{
    int ncoef = cv.numCoefs();
    int dim = cv.dimension();
    vector<double>::const_iterator coef = cv.coefs_begin();

    if (idx >= dim)
	return false;   // Question does not make sense

    bool identical = true;
    for (int ki=0; ki<ncoef; ki++, coef+=dim)
    {
	double d1 = fabs(coef[idx]-val);
	if (d1 > max_dist)
	    return false;   // Curve does not lie close the given value
	if (d1 > tol)
	    identical = false;
    }
    return (!identical);
}

//-------------------------------------------------------------------
 void GeometryTools::setSfBdCoefToConst(SplineSurface& srf, int bd_idx, int idx_d, double val,
			 double deg_tol)

// Modify surface along specified boundary to match a specific constant in 
// one direction. Currently only non-rational surfaces
//-------------------------------------------------------------------
{
    bool dir_u = (bd_idx == 2 || bd_idx == 3);
    int nmb1 = (dir_u) ? srf.numCoefs_u() : srf.numCoefs_v();
    int nmb2 = (dir_u) ? srf.numCoefs_v() : srf.numCoefs_u();
    int dim = srf.dimension();
    vector<double>::iterator coefs = srf.coefs_begin();
    int kr1 = (dir_u) ? dim : nmb1*dim;
    int kr2 = (dir_u) ? nmb2*dim : dim;

    if (idx_d >= dim)
	return;  // Does not make sense

    // Check degeneracy
    bool bd[] = { false, false, false, false }; // left, right, bottom, top
    // bool deg = srf.isDegenerate(bd[2], bd[1], bd[3], bd[0], deg_tol);
    bool start_deg = ((bd_idx <=1 && bd[2]) || (bd_idx > 1 && bd[0]));
    bool end_deg = ((bd_idx <=1 && bd[3]) || (bd_idx > 1 && bd[1]));

    // Modify coefficients
    vector<double>::iterator c1 = coefs;
    if (bd_idx == 3)
	c1 += (nmb2-1)*nmb1*dim;
    if (bd_idx == 1)
	c1 += (nmb2-1)*dim;

    for (int ki=0; ki<nmb1; ki++, c1+=kr1)
    {
	c1[idx_d] = val;
    }

    if (start_deg)
    {
	// Make sure that the degeneracy is maintained
	vector<double>::iterator d1 = coefs;
	for (int ki=0; ki<nmb2; ki++, d1+=kr2)
	{
	    d1[idx_d] = val;
	}
    }

    if (end_deg)
    {
	// Make sure that the degeneracy is maintained
	vector<double>::iterator d1 = coefs;
	if (bd_idx <= 1)
	    d1 += (nmb1-1)*nmb2*dim;
	else
	    d1 += (nmb1-1)*dim;
	for (int ki=0; ki<nmb2; ki++, d1+=kr2)
	{
	    d1[idx_d] = val;
	}
    }

}

//-------------------------------------------------------------------
void GeometryTools::getGnJoints(const ParamCurve& curve, const vector<double>& cont, vector<double>& gn_joints)
//-------------------------------------------------------------------
{
    int n = (int)cont.size() - 1;
    // If input should be of higher order we should let kink be a corr vector.
    ASSERT(n >= 0 && n < 3); // Currently only implemented for pos, tangent and 2nd order geom cont.

    int ki;
    vector<double> poss_disc; // We store all possible discontinuities.

    // We begin by extracting all parameter values (knots) which 
    const SplineCurve* spline_cv;
    if (curve.instanceType() == Class_SplineCurve) {
	spline_cv = dynamic_cast<const SplineCurve*>(&curve);
    } else if (curve.instanceType() == Class_CurveOnSurface) {
	const CurveOnSurface* cv_on_sf = dynamic_cast<const CurveOnSurface*>(&curve);
	if (cv_on_sf->parPref()) {
	    spline_cv = dynamic_cast<const SplineCurve*>((cv_on_sf->parameterCurve()).get());
	} else {
	    spline_cv = dynamic_cast<const SplineCurve*>((cv_on_sf->spaceCurve()).get());
	}
    } else {
	THROW("Unexpected curve type!");
    }
    ASSERT(spline_cv != 0);
    int order = spline_cv->order();
    vector<double>::const_iterator iter = spline_cv->basis().begin() + order;
    int mult = order - n; // If mult equal knots, we split.
    while ((iter < spline_cv->basis().end()) && (iter[0] != spline_cv->endparam())) {
	if (iter[mult-1] == iter[0]) {
	    poss_disc.push_back(iter[0]);
	    iter += mult;
	} else {
	    ++iter;
	}
    }

    // We loop through all possible discontinuities, checking the value.
    gn_joints.push_back(spline_cv->startparam());
    int pds = (int)poss_disc.size();
    for (ki = 0; ki < pds; ++ki) {
	vector<Point> from_left = spline_cv->ParamCurve::point(poss_disc[ki], n, false);
	vector<Point> from_right = spline_cv->ParamCurve::point(poss_disc[ki], n);
	ASSERT(int(from_left.size()) == n+1);
	double dist = from_left[0].dist(from_right[0]);
	if (dist > cont[0]) // Testing for G0 cont
	  {
	    gn_joints.push_back(poss_disc[ki]);
	    continue;
	  }

	if (n > 0) { // Testing for G1 cont (direction of tangents).
	    from_left[1].normalize();
	    from_right[1].normalize();
	    dist = from_left[1].dist(from_right[1]);
	    if (dist > cont[1])
	      {
		gn_joints.push_back(poss_disc[ki]);
		continue;
	      }

	    if (n > 1) { // Testing for G2 cont
		// (We may not need to split into segments, but as the algorithm expects 
		// it we do anyway.
		// Would be smarter to save segments as we proceed, but since when was
		// speed an issue ...)
		shared_ptr<SplineCurve> left_seg
		    (spline_cv->subCurve(spline_cv->startparam(), poss_disc[ki]));
		shared_ptr<SplineCurve> right_seg
		    (spline_cv->subCurve(poss_disc[ki], spline_cv->endparam()));

		// Let b1, b2 & b3 be the last 3 ctrl pts on left seg, and c0, c1, c2 the
		// first 3 ctrl pts
		// on the right seg. Since we are guaranteed G0 cont we know that b0 = c0.
		// Assuming they are coplanar, we define d as the intersection point
		// between the lines given by (b1, b2) & (c1, c2) (i.e. passing through the pts).
		// If all pts are colinear we let d = b0 (= c0).
		// Defining
		// r_- = ratio(b1, b2, d)
		// r   = ratio(b2, b3, c1) (= ratio(b2, c0, c1)
		// r_+ = ratio(d, c1, c2),
		// where ratio(a, b, c) = |b-a|/|c-b| for the colinear pts a, b & c.
		// By [Gerald Farin: Curves and Surfaces for CAGD, 12.2 (p.182)] G2 continuity is
		// given by the formula
		// r^2 = r_- * r_+.
		// Using the sinus sentence we calculate the unknown |d - b2| & |c1 - d|.
		int dim = left_seg->dimension();
		Point b1(left_seg->coefs_end() - 3*dim, left_seg->coefs_end() - 2*dim);
		Point b2(left_seg->coefs_end() - 2*dim, left_seg->coefs_end() - dim);
		Point b3(left_seg->coefs_end() - dim, left_seg->coefs_end());
		Point c0(right_seg->coefs_begin(), right_seg->coefs_begin() + dim);
		Point c1(right_seg->coefs_begin() + dim, right_seg->coefs_begin() + 2*dim);
		Point c2(right_seg->coefs_begin() + 2*dim, right_seg->coefs_begin() + 3*dim);

		Point left_dir = b2 - b1;
		Point mid_dir = c1 - b2;
		Point right_dir = c1 - c2;
		double ang_left, ang_right;
		try {
		  ang_left = left_dir.angle_smallest(mid_dir);
		  ang_right = right_dir.angle_smallest(mid_dir);
		}
		catch ( ... ) {
		  MESSAGE("Failed computing angle between joints, skipping.");
		  continue;
		}
		double ang_mid = M_PI - ang_left - ang_right;
		double sin_left = sin(ang_left);
		double sin_right = sin(ang_right);
		double sin_mid = sin(ang_mid);
		double ratio_left = left_dir.length()*sin_mid/(mid_dir.length()*sin_right);
		double ratio_mid = b2.dist(b3)/c0.dist(c1);
		double ratio_right = mid_dir.length()*sin_left/(right_dir.length()*sin_mid);

		double error = ratio_mid*ratio_mid - ratio_left*ratio_right;
		if (error > cont[2]) // @@sbr Not sure what estimate to apply.
		  {
		    gn_joints.push_back(poss_disc[ki]);
		    continue;
		  }
	    }
	}
	// Since we got here all tests passed! This is NOT a segment
	// gn_joints.push_back(poss_disc[ki]);
    }
    gn_joints.push_back(spline_cv->endparam());

    // We run through gn_joints checking length of segments.
    if (gn_joints.size() > 2) {
	for (ki = 0; ki < int(gn_joints.size()) - 1; ++ki) {
	    double segment_length =
		curve.ParamCurve::estimatedCurveLength(gn_joints[ki], gn_joints[ki+1], 1000);
	    if (segment_length < cont[0]) {
		if (ki == 0) {
		    gn_joints.erase(gn_joints.begin() + ki + 1);
		} else {
		    gn_joints.erase(gn_joints.begin() + ki);
		}
		--ki;
	    }
	    if (gn_joints.size() == 2) {
		break;
	    }
	}
    }
}

//-------------------------------------------------------------------
void GeometryTools::getGnJoints(const CurveLoop& loop, const vector<double>& cont,
		 vector<vector<double> >& gn_joints)
//-------------------------------------------------------------------
{
    int n = (int)cont.size() - 1;

    int ki, kj;
    vector<vector<double> > poss_disc;
    for (ki = 0; ki < loop.size(); ++ki) {
	vector<double> gn_joints;
	GeometryTools::getGnJoints(*(loop[ki]), cont, gn_joints);
	poss_disc.push_back(gn_joints);
    }

    // We then must check cont between consecutive segments.
    int pds = (int)poss_disc.size();
    for (ki = 0; ki < pds; ++ki) {
	int fi = ki;
	int si = (fi + 1)%pds;
	vector<Point> from_left = loop[fi]->ParamCurve::point(poss_disc[fi].back(), n, false);
	vector<Point> from_right = loop[si]->ParamCurve::point(poss_disc[si].front(), n);
	ASSERT(int(from_left.size()) == n+1);
	for (kj = 0; kj < n+1; ++kj) {
	    if (kj != 0) {
		from_left[kj].normalize();
		from_right[kj].normalize();
	    }
	    double dist = from_left[kj].dist(from_right[kj]);
	    if (dist > cont[kj])
		break;
	}
	poss_disc[fi].erase(poss_disc[fi].end() - 1);
	if (kj == n+1) {
	    poss_disc[si].erase(poss_disc[si].begin());
	}
    }

    gn_joints = poss_disc;
}

// Describe a surface as a curve in a given direction
//-------------------------------------------------------------------
shared_ptr<SplineCurve>
GeometryTools::representSurfaceAsCurve(const SplineSurface& surface, int cv_dir)
//-------------------------------------------------------------------
{
  int dim = surface.dimension();
  int kdim = dim + (surface.rational() ? 1 : 0);
  const std::vector<double>::const_iterator co = (surface.rational()) ?
    surface.rcoefs_begin() : surface.coefs_begin();
  std::vector<double> huge_curve_coefs;
  std::vector<double>::const_iterator coefstart; 
  if (cv_dir != 2) {
    // We must flip the surface coefficients
    int nu = surface.numCoefs_u();
    int nv = surface.numCoefs_v();
    huge_curve_coefs.reserve(nu*nv*kdim);
    for (int i = 0; i < nu; ++i)
      for (int j = 0; j < nv; ++j)
	for (int k = 0; k < kdim; ++k)
	  huge_curve_coefs.push_back(co[(j*nu+i)*kdim + k]);
    coefstart = huge_curve_coefs.begin();
  } else {
    coefstart = co;
  }
  
    const BsplineBasis& bas = surface.basis(2-cv_dir);
    const BsplineBasis& other_bas = surface.basis(cv_dir-1);

    int num = other_bas.numCoefs();
    int order = other_bas.order();
    std::vector<double>::const_iterator knotstart = other_bas.begin();

    int dim2 = bas.numCoefs() * kdim;
    shared_ptr<SplineCurve> curve(new SplineCurve(num, order,
						  knotstart, coefstart,
						  dim2, 
						  false));
    return curve;
}

// Describe a curve as a surface in a given direction
//-------------------------------------------------------------------
shared_ptr<SplineSurface>
GeometryTools::representCurveAsSurface(const SplineCurve& curve,
			int cv_dir,
			const BsplineBasis& other_bas,
			bool rational)
//-------------------------------------------------------------------
{
    if (curve.rational()) {
	THROW("It does not make sense to have a rational hypercurve.");
    }
  int kdim = curve.dimension()/other_bas.numCoefs();
  std::vector<double>::const_iterator co
      = curve.coefs_begin();
  std::vector<double> surf_coefs;
  std::vector<double>::const_iterator coefstart;
  const BsplineBasis& bas = curve.basis();
  if (cv_dir!=2) {
    // We must flip the curve coefficients
    int nu = curve.numCoefs();
    int nv = other_bas.numCoefs();
    surf_coefs.reserve(nu*nv*kdim);
    for (int i = 0; i < nv; ++i)
      for (int j = 0; j < nu; ++j)
	for (int k = 0; k < kdim; ++k)
	  surf_coefs.push_back(co[(j*nv+i)*kdim + k]);
    coefstart = surf_coefs.begin();
  } else {
    coefstart = co;
  }
  int dim = rational ? kdim - 1 : kdim;
  shared_ptr<SplineSurface> surface
    (new SplineSurface((cv_dir==1) ? bas : other_bas,
			 (cv_dir==1) ? other_bas : bas, 
			 coefstart, dim, rational));
  return surface;
}

//-------------------------------------------------------------------
shared_ptr<SplineSurface>
GeometryTools::joinPatches(const vector<shared_ptr<SplineSurface> >& patches,
	    const SplineSurface& spline_space)
//-------------------------------------------------------------------
{
    // As we want to know how many patches there exist in u- and v-direction we
    // insert knots into copy of input sf.
    SplineSurface cp_spline_sf = spline_space;

    cp_spline_sf.makeBernsteinKnotsU();
    cp_spline_sf.makeBernsteinKnotsV();

    const BsplineBasis& orig_basis_u = spline_space.basis_u();
    const BsplineBasis& orig_basis_v = spline_space.basis_v();
    const BsplineBasis& basis_u = cp_spline_sf.basis_u();
    const BsplineBasis& basis_v = cp_spline_sf.basis_v();

    // Assuming both input basises are represented with Bezier knots.
    // (i.e. all knots should have mult equal to degree).
    int num_u = basis_u.numCoefs();
    int num_v = basis_v.numCoefs();
    int order_u = basis_u.order();
    int order_v = basis_v.order();
    int numpat_u = num_u / order_u; //num_u + 1 - order_u
    int numpat_v = num_v / order_v;; //num_v + 1 - order_v;

    bool repar = false; // @@sbr Alternatively we could update input knots (spline_space) as we append.
    double fuzzy_tol = 1e-06;
    SplineSurface* total_sf;
    for (int ki = 0; ki < numpat_v; ++ki) {
	SplineSurface sf_u = *patches[ki*numpat_u];
	for (int kj = 1; kj <numpat_u; ++kj) {
	    // @@sbr Perhaps make sure that we are not fooled by numerical noise?
	    // suppose this could be done somewhat more efficient, but it'll do for now.
	    double knot_u = sf_u.endparam_u();
	    orig_basis_u.knotIntervalFuzzy(knot_u, fuzzy_tol);
	    int mult_u = orig_basis_u.knotMultiplicity(knot_u);
	    int cont = order_u - mult_u - 1;
	    double dist;
	    sf_u.appendSurface(patches[ki*numpat_u+kj].get(), 1, cont, dist, repar);
	}

	if (ki == 0) {
	    total_sf = sf_u.clone();
	} else {
	    double knot_v = total_sf->endparam_v();
	    orig_basis_v.knotIntervalFuzzy(knot_v, fuzzy_tol);
	    int mult_v = orig_basis_v.knotMultiplicity(knot_v);
	    int cont = order_v - mult_v - 1;
	    double dist;
	    total_sf->appendSurface(&sf_u, 2, cont, dist, repar);
	}
    }

    return shared_ptr<SplineSurface>(total_sf);
}

//-------------------------------------------------------------------
void GeometryTools::surfaceKinks(const SplineSurface& sf, double max_normal_angle,
		  vector<double>& g1_disc_u, vector<double>& g1_disc_v,
		  bool compute_g1_disc)
//-------------------------------------------------------------------
{
    int dim = sf.dimension();
    bool rational = sf.rational();
    int kdim = dim + (rational);
    ALWAYS_ERROR_IF(dim != 3,
		"Expecting surface to be of dimension 3.");
    g1_disc_u.clear();
    g1_disc_v.clear();

    bool dir_is_u = true; // We traverse dir given by dir_is_u, looking for kinks.
    for (int ki = 0; ki < 2; ++ki) { // We perform calculations in both parameter directions.
        const BsplineBasis& basis = sf.basis(!dir_is_u);
	const BsplineBasis& other_basis = sf.basis(dir_is_u);
	std::vector<double> disc_candidates;
	int order = basis.order();
	int other_order = other_basis.order();
	int num_coefs = basis.numCoefs();
	int other_num_coefs = other_basis.numCoefs();
	std::vector<double> coefs(rational ? sf.rcoefs_begin() : sf.coefs_begin(), 
				  rational ? sf.rcoefs_end() : sf.coefs_end());
	if (!dir_is_u) // We may now assume that direction is u.
	    SplineUtils::transpose_array(kdim, num_coefs, other_num_coefs, &coefs[0]);

	int nmb_int_samples = (rational) ?
	    9*other_order - 11 : 4*other_order - 5; // # of evaluations for each knot interval

	if (!compute_g1_disc)
	{
	    // Check for parametric continuity. Fewer evaluations are required
	    nmb_int_samples = (rational) ? 4*other_order - 1 : other_order-1;
	}

	int g1_mult = order - 2; // Multiplicity for which g1 is assured.
	std::vector<double>::const_iterator iter = basis.begin() + order;
	std::vector<double>::const_iterator end_iter = basis.end() - order;
	for (; iter != end_iter; ++iter)
	    if (iter[0] == iter[g1_mult]) {
		disc_candidates.push_back(*iter);
		while (iter[0] == iter[1])
		    ++iter;
	    }

	// For each possible discontinuity we calculate angles between normals.
	for (size_t kj = 0; kj < disc_candidates.size(); ++kj) {
	    double tpar = disc_candidates[kj];
	    int derivs = 1;
	    std::vector<double> basis_vals_left = // Computing from the left.
		basis.computeBasisValuesLeft(tpar, derivs);
	    std::vector<double> basis_vals_right = // Computing from the right.
		basis.computeBasisValues(tpar, derivs);
	    int coef_li = basis.knotIntervalFuzzy(tpar) - order + 1;
	    int mult = basis.knotMultiplicity(tpar);
	    int coef_left_li = coef_li - mult;

	    // To increase performance we avoid recalculating basis values when next
	    // knots is in the same interval.
	    std::vector<double>::const_iterator iter = other_basis.begin();
	    while (iter != other_basis.end()) {
		std::vector<double>::const_iterator next_iter = iter + 1;
		while ((next_iter < other_basis.end()) && (*next_iter == *iter))
		    ++next_iter;
		if (next_iter == other_basis.end()) {
		    break;
		}
		double tstep = (*next_iter - *iter)/(nmb_int_samples - 1);
		double tmin = other_basis.startparam();
		int kh;
		for (kh = 0; kh < nmb_int_samples; ++kh) {
		    double other_tpar = tmin + kh*tstep;
		    // We must calculate normals from both sides of the iso-curve.
		    std::vector<double> other_basis_vals =
			other_basis.computeBasisValues(other_tpar, derivs);
		    int other_coef_li = other_basis.knotIntervalFuzzy(other_tpar) - other_order + 1;

		    Point part_left(kdim);
		    part_left.setValue(0.0);
		    Point part_right = part_left;
		    Point other_part = part_left;
		    Point pt = part_left;
		    // We compute the partial derivatives in both dirs, assuming surface is cont.
		    for (int ii = 0; ii < order; ++ii)
			for (int jj = 0; jj < other_order; ++jj)
			    for (int kk = 0; kk < kdim; ++kk) {
				double coef =
				    coefs[kdim*((other_coef_li+jj)*num_coefs+coef_li+ii)+kk];
				double coef_l =
				    coefs[kdim*((other_coef_li+jj)*num_coefs+coef_left_li+ii)+kk];
				pt[kk] +=
				    coef*basis_vals_right[2*ii]*other_basis_vals[2*jj];
				part_left[kk] +=
				    coef_l*basis_vals_left[2*ii+1]*other_basis_vals[2*jj];
				part_right[kk] +=
				    coef*basis_vals_right[2*ii+1]*other_basis_vals[2*jj];
				other_part[kk] +=
				    coef*basis_vals_left[2*ii]*other_basis_vals[2*jj+1];
			    }
		    if (rational) {
			std::vector<double> restemp(3*kdim);
			std::copy(pt.begin(), pt.end(), restemp.begin());
			std::copy(part_left.begin(), part_left.end(), restemp.begin() + kdim);
			std::copy(other_part.begin(), other_part.end(), restemp.begin() + 2*kdim);
			std::vector<double> restemp2(3*dim);
			SplineUtils::surface_ratder(&restemp[0], dim, 1, &restemp2[0]);
			part_left.resize(dim);
			part_left.setValue(&restemp2[dim]);
			std::copy(part_right.begin(), part_right.end(), restemp.begin() + kdim);
			SplineUtils::surface_ratder(&restemp[0], dim, 1, &restemp2[0]);
			part_right.resize(dim);
			part_right.setValue(&restemp2[dim]);
			other_part.resize(dim);
			other_part.setValue(&restemp2[2*dim]);
		    }

		    double angle;
		    if (compute_g1_disc)
		    {
			// Compute angle between surface normals
			Point left_normal = part_left%other_part;
			Point right_normal = part_right%other_part;
			angle = left_normal.angle(right_normal);
		    }
		    else
		    {
			// Compute difference between partial derivatives 
			angle = part_left.dist(part_right);
		    }

		    if (angle > max_normal_angle) {
			if (dir_is_u)
			    g1_disc_u.push_back(tpar);
			else
			    g1_disc_v.push_back(tpar);
			break;
		    }
		}
		if (kh < nmb_int_samples)
		    break;
		else
		    iter = next_iter;
	    }
	}
	dir_is_u = false;
    }
}

//-------------------------------------------------------------------
std::vector<shared_ptr<SplineSurface> > 
GeometryTools::splitInKinks(const SplineSurface& sf,
	     const std::vector<double>& u_kinks,
	     const std::vector<double>& v_kinks)
//-------------------------------------------------------------------
{
    std::vector<shared_ptr<SplineSurface> > return_sfs;

    std::vector<double> split_params_u;
    std::vector<double> split_params_v;
    split_params_u.push_back(sf.startparam_u());
    split_params_u.insert(split_params_u.end(), u_kinks.begin(), u_kinks.end());
    split_params_u.push_back(sf.endparam_u());
    split_params_v.push_back(sf.startparam_v());
    split_params_v.insert(split_params_v.end(), v_kinks.begin(), v_kinks.end());
    split_params_v.push_back(sf.endparam_v());
    int num_u = (int)split_params_u.size();
    int num_v = (int)split_params_v.size();
    for (int ki = 0; ki < num_u - 1; ++ki)
	for (int kj = 0; kj < num_v - 1; ++kj)
	  {
// #if _MSC_VER > 0 && _MSC_VER < 1300
// 	    return_sfs.push_back(shared_ptr<SplineSurface>
// 				 (dynamic_cast<SplineSurface*>
// 				  (sf.subSurface(split_params_u[ki], split_params_v[kj],
// 						  split_params_u[ki+1], split_params_v[kj+1]))));
// #else	    
	    return_sfs.push_back(shared_ptr<SplineSurface>
				 (sf.subSurface(split_params_u[ki], split_params_v[kj],
						 split_params_u[ki+1], split_params_v[kj+1])));
// #endif
	  }
    
    return return_sfs;
}

void GeometryTools::curveKinks(const SplineCurve& cv, double tol, double ang_tol,
		vector<double>& c1_disconts, vector<double>& g1_disconts)
{

    // Get candidate discontinuity parameters
    vector<double> cand_disconts;
    cv.basis().cNDiscontinuities(cand_disconts, 1);

    // Compare results of left- and right evaluation in the candidates
    int der = 1;
    int dim = cv.dimension();
    vector<Point> pnt_left(2, Point(dim)), pnt_right(2, Point(dim));
    for (size_t ki=0; ki<cand_disconts.size(); ++ki)
    {
	cv.point(pnt_left, cand_disconts[ki], der, false);
	cv.point(pnt_right, cand_disconts[ki], der, true);
	if (pnt_left[1].dist(pnt_right[1]) > tol)
	    c1_disconts.push_back(cand_disconts[ki]);
	
	if (pnt_left[1].angle(pnt_right[1]) > ang_tol ||
	    pnt_left[1]*pnt_right[1] < 0.0)
	    g1_disconts.push_back(cand_disconts[ki]);
    }
}



//-------------------------------------------------------------------
bool GeometryTools::isCoincident(const ParamCurve& cv1, const ParamCurve& cv2, double epsge)
//-------------------------------------------------------------------
{
  // Check endpoints
  double ta1 = cv1.startparam();
  double tb1 = cv1.endparam();
  double ta2 = cv2.startparam();
  double tb2 = cv2.endparam();

  Point pnt1 = cv1.point(ta1); 
  Point pnt2 = cv1.point(tb1);
  Point pnt3 = cv2.point(ta2);
  Point pnt4 = cv2.point(tb2);
  
  if (std::min(pnt1.dist(pnt3), pnt1.dist(pnt4)) > epsge ||
      std::min(pnt2.dist(pnt3), pnt2.dist(pnt4)) > epsge)
    return false;

  // Endpoints matching. Check the inner.
  bool same = (pnt1.dist(pnt3) < pnt1.dist(pnt4));

  double length = cv1.estimatedCurveLength();
  int nmb_test = std::min(500, (int)(length/epsge));
  double del = (tb1 - ta1)/(int)(nmb_test + 1);
  double del2 = (tb2 - ta2)/(int)(nmb_test + 1);
  double tpar, tpar2, seed, dist;
  int ki;
  for (ki=0, tpar=ta1+del; ki<nmb_test; ki++, tpar+=del)
    {
      seed = (same) ? ta2 + (ki+1)*del2 : tb2 - (ki+1)*del2;
      pnt1 = cv1.point(tpar);
      cv2.closestPoint(pnt1, ta2, tb2, tpar2, pnt2, dist, &seed);
      if (dist > epsge)
	return false;
    }

  return true;
}



//-------------------------------------------------------------------
void GeometryTools::translateSplineSurf(const Point& trans_vec, SplineSurface& sf)
//-------------------------------------------------------------------
{
    int ki;
    ASSERT(trans_vec.dimension() == 3); // We're working in 3D space.
    int dim = 3 + sf.rational();
    std::vector<double>::iterator iter = sf.rational() ? sf.rcoefs_begin() : sf.coefs_begin();
    std::vector<double>::iterator end_iter = sf.rational() ? sf.rcoefs_end() : sf.coefs_end();
    std::vector<double>::iterator coef_iter = sf.coefs_begin();
    while (iter != end_iter) {
	double weight = sf.rational() ? iter[3] : 1.0;
	for (ki = 0; ki < 3; ++ki)
	    iter[ki] = (iter[ki]/weight + trans_vec[ki])*weight;
	if (sf.rational()) {
	    for (ki = 0; ki < 3; ++ki)
		coef_iter[ki] = iter[ki]/weight;
	    coef_iter += dim - 1;
	}
	iter += dim;
    }
}

//-------------------------------------------------------------------
void GeometryTools::translateSplineCurve(const Point& trans_vec, SplineCurve& cv)
//-------------------------------------------------------------------
{
    int ki;
    int dim0 = trans_vec.dimension();
    ASSERT(dim0 == 2 || dim0 == 3); // We're working in 2D or 3D space.
    int dim = dim0 + cv.rational();
    std::vector<double>::iterator iter = cv.rational() ? cv.rcoefs_begin() : cv.coefs_begin();
    std::vector<double>::iterator end_iter = cv.rational() ? cv.rcoefs_end() : cv.coefs_end();
    std::vector<double>::iterator coef_iter = cv.coefs_begin();
    while (iter != end_iter) { // @@ A faster approach would be to use rotation matrix directly.
	double weight = cv.rational() ? iter[dim0] : 1.0;
	for (ki = 0; ki < dim0; ++ki)
	    iter[ki] = (iter[ki]/weight + trans_vec[ki])*weight;
	if (cv.rational()) {
	    for (ki = 0; ki < dim0; ++ki)
		coef_iter[ki] = iter[ki]/weight;
	    coef_iter += dim - 1;
	}
	iter += dim;
    }
}


//-------------------------------------------------------------------
void GeometryTools::translateLineCloud(const Point& trans_vec, LineCloud& lc)
//-------------------------------------------------------------------
{
    int ki;
    int nmb_pts = lc.numLines()*2;
    int dim = 3;
    for (ki = 0; ki < nmb_pts; ++ki) {
	Vector3D vec = lc.point(ki);
	Point pt(vec.begin(), vec.end());
	pt += trans_vec;
	copy(pt.begin(), pt.end(), lc.rawData() + ki*dim);
    }
}


//-------------------------------------------------------------------
void GeometryTools::rotateSplineSurf(Point rot_axis, double alpha, SplineSurface& sf)
//-------------------------------------------------------------------
{
    rot_axis.normalize();
    ASSERT(rot_axis.dimension() == 3); // We're working in 3D space.
    int dim = 3 + sf.rational();
    std::vector<double>::iterator iter = sf.rational() ? sf.rcoefs_begin() : sf.coefs_begin();
    std::vector<double>::iterator end_iter = sf.rational() ? sf.rcoefs_end() : sf.coefs_end();
    std::vector<double>::iterator coef_iter = sf.coefs_begin();
    while (iter != end_iter) { // @@ A faster approach would be to use rotation matrix directly.
	rotatePoint(rot_axis, alpha, &*iter);
	if (sf.rational()) {
	    for (int ki = 0; ki < 3; ++ki)
		coef_iter[ki] = iter[ki]/iter[3];
	    coef_iter += dim - 1;
	}
	iter += dim;
    }
}


//-------------------------------------------------------------------
void GeometryTools::rotateSplineCurve(Point rot_axis, double alpha, SplineCurve& cv)
//-------------------------------------------------------------------
{
    rot_axis.normalize();
    int dim0 = rot_axis.dimension();
    ASSERT(dim0 == 2 || dim0 == 3); // We're working in 2D or 3D space.
    int dim = dim0 + cv.rational();
    std::vector<double>::iterator iter = cv.rational() ? cv.rcoefs_begin() : cv.coefs_begin();
    std::vector<double>::iterator end_iter = cv.rational() ? cv.rcoefs_end() : cv.coefs_end();
    std::vector<double>::iterator coef_iter = cv.coefs_begin();
    while (iter != end_iter) { // @@ A faster approach would be to use rotation matrix directly.
	rotatePoint(rot_axis, alpha, &*iter);
	if (cv.rational()) {
	    for (int ki = 0; ki < dim0; ++ki)
		coef_iter[ki] = iter[ki]/iter[3];
	    coef_iter += dim - 1;
	}
	iter += dim;
    }
}


//-------------------------------------------------------------------
void GeometryTools::rotateLineCloud(Point rot_axis, double alpha, LineCloud& lc)
//-------------------------------------------------------------------
{
    rot_axis.normalize();
    ASSERT(rot_axis.dimension() == 3); // We're working in 3D space.
    int ki;
    int nmb_pts = lc.numLines()*2;
    int dim = 3;
    for (ki = 0; ki < nmb_pts; ++ki) {
	Vector3D vec = lc.point(ki);
	Point pt(vec.begin(), vec.end());
	rotatePoint(rot_axis, alpha, pt.begin());
	copy(pt.begin(), pt.end(), lc.rawData() + ki*dim);
    }
}


//-------------------------------------------------------------------
void GeometryTools::rotatePoint(Point rot_axis, double alpha, double* space_pt)
//-------------------------------------------------------------------
{
  int dim = rot_axis.dimension();
  ASSERT(dim == 2 || dim == 3);
  // We're working in 2D or 3D space. The axis content is used only in 3D
  if (dim == 3)
    rot_axis.normalize();

  std::vector<double> rot_mat = GeometryTools::getRotationMatrix(rot_axis, alpha);
  int ki, kj;
  Point rotated_pt(dim);
  for (ki = 0; ki < dim; ++ki)
    for (kj = 0; kj < dim; ++kj)
      rotated_pt[ki] += rot_mat[ki*dim+kj]*space_pt[kj];
  
  for (ki = 0; ki < dim; ++ki)
    space_pt[ki] = rotated_pt[ki];
}


//-------------------------------------------------------------------
void GeometryTools::rotatePoint(Point rot_axis, double alpha, Point& space_pt)
//-------------------------------------------------------------------
{
    rotatePoint(rot_axis, alpha, space_pt.begin());
}


//-------------------------------------------------------------------
std::vector<double> GeometryTools::getRotationMatrix(const Point& unit_rot_axis, double alpha)
//-------------------------------------------------------------------
{
  int dim = unit_rot_axis.dimension();
  ASSERT(dim == 2 || dim == 3);
  // We're working in 2D or 3D space. The axis content is used only in 3D
  std::vector<double> return_matrix;
  double cos_a = cos(alpha);
  double sin_a = sin(alpha);
  if (dim == 2)
    {
      return_matrix.resize(4);
      return_matrix[0] = cos_a;
      return_matrix[1] = -sin_a;
      return_matrix[2] = sin_a;
      return_matrix[3] = cos_a;
    }
  else
    {
      return_matrix.resize(9);
      // Using the notation in 'Computational Gemoetry for Design and Manufacture'.
      double u1 = unit_rot_axis[0];
      double u2 = unit_rot_axis[1];
      double u3 = unit_rot_axis[2];

      return_matrix[0] = u1*u1 + cos_a*(1 - u1*u1);
      return_matrix[1] = u1*u2*(1 - cos_a) - u3*sin_a;
      return_matrix[2] = u3*u1*(1 - cos_a) + u2*sin_a;
      return_matrix[3] = u1*u2*(1 - cos_a) + u3*sin_a;
      return_matrix[4] = u2*u2 + cos_a*(1 - u2*u2);
      return_matrix[5] = u2*u3*(1 - cos_a) - u1*sin_a;
      return_matrix[6] = u3*u1*(1 - cos_a) - u2*sin_a;
      return_matrix[7] = u2*u3*(1 - cos_a) + u1*sin_a;
      return_matrix[8] = u3*u3 + cos_a*(1 - u3*u3);
    }
    return return_matrix;
}


//===========================================================================
void 
GeometryTools::insertKnotsEvenly(BsplineBasis& basis, int num_knots)
//===========================================================================
{
    for (int i=0; i<num_knots; ++i) {
	double new_knot = GeometryTools::getKnotAtLargestInterval(basis);
	basis.insertKnot(new_knot);
    }
}


//===========================================================================
void GeometryTools::insertKnotsEvenly(BsplineBasis& basis, double tmin, double tmax,
		       int num_knots, double knot_diff_tol)
//===========================================================================
{
    BsplineBasis sub_basis = basis.subBasis(tmin, tmax, knot_diff_tol);
    vector<double> new_knots;
    for (int i=0; i<num_knots; ++i) {
	double new_knot = getKnotAtLargestInterval(sub_basis);
	sub_basis.insertKnot(new_knot);
	new_knots.push_back(new_knot);
    }

    basis.insertKnot(new_knots);
}


//===========================================================================
double GeometryTools::getKnotAtLargestInterval(const BsplineBasis& basis)
//===========================================================================
{
    const pair<double, double> interval = 
	GeometryTools::getLargestParameterInterval(basis);
    double new_knot = 0.5*(interval.first+interval.second);

    return new_knot;
}


//===========================================================================
pair<double, double>
GeometryTools::getLargestParameterInterval(const BsplineBasis& basis)
//===========================================================================
{
    double largest = 0.0;
    std::vector<double>::const_iterator index;
    for (std::vector<double>::const_iterator knot = basis.begin();
	 knot != basis.end()-1; ++knot) {
	if ( (*(knot+1) - *knot) > largest) {
	    largest = *(knot+1) - *knot;
	    index = knot;
	}
    }

    return pair<double, double>(*index, *(index+1));
}

//===========================================================================
void
GeometryTools::averageCoefsAtSeam(shared_ptr<SplineSurface>& srf, int dir,
				  bool c1_cont)
//===========================================================================
{
  if (dir == 1)
    srf->swapParameterDirection();

  int dim = srf->dimension();
  int in1 = srf->numCoefs_u();
  int in2 = srf->numCoefs_v();
  vector<double>::iterator c1 = srf->coefs_begin();
  vector<double>::iterator c2 = srf->coefs_begin();
  c2 += (in2-1)*in1*dim;
  for (int ki=0; ki<in1*dim; ++ki, ++c1, ++c2)
    {
      double tmid = 0.5*((*c1)+(*c2));
      *c1 = tmid;
      *c2 = tmid;
    }

  if (c1_cont)
    {
      c1 = srf->coefs_begin();
      c2 = srf->coefs_begin();
      c2 += (in2-1)*in1*dim;
      vector<double>::iterator c1_2 = c1 + in1*dim;
      vector<double>::iterator c2_2 = c2 - in1*dim;
      double vpar1 = srf->startparam_v();
      double vpar2 = srf->endparam_v();
      for (int ki=0; ki<in1; ++ki, c1+=dim, c2+=dim)
	{
	  double upar = srf->basis_u().grevilleParameter(ki);
	  vector<Point> pts1(3), pts2(3);
	  srf->point(pts1, upar, vpar1, 1);
	  srf->point(pts2, upar, vpar2, 1, false, true);
	  Point norm1 = pts1[1].cross(pts1[2]);
	  Point norm2 = pts2[1].cross(pts2[2]);
	  norm1.normalize();
	  norm2.normalize();
	  Point norm = 0.5*(norm1+norm2);
	  norm.normalize();

	  Point c3(c1_2, c1_2+dim);
	  Point c4(c2_2, c2_2+dim);
	  Point c0(c1, c1+dim);
	  
	  c3 -= (c3-c0)*norm*norm;
	  c4 -= (c4-c0)*norm*norm;
	  int ka;
	  for (ka=0; ka<dim; ++ka, c1_2++, c2_2++)
	    {
	      *c1_2 = c3[ka];
	      *c2_2 = c4[ka];
	    }
	  int stop_break = 1;
	}
   }

  if (srf->rational())
    {
	// This fix will not work for all configuration of rational surfaces.
	// Update the rational coefficients with respect to the divided
	// ones
	c1 = srf->coefs_begin();
	vector<double>::iterator r1 = srf->rcoefs_begin();
	int kn = in1*in2;
	for (int ki=0; ki<kn; ++ki)
	  {
	    for (int kr=0; kr<dim; ++kr)
	      r1[kr] = c1[kr]*r1[dim];
	    c1 += dim;
	    r1 += (dim+1);
	  }
      }

  if (dir == 1)
    srf->swapParameterDirection();
}

//===========================================================================
void
GeometryTools::averageBoundaryCoefs(shared_ptr<SplineSurface>& srf1, int bd1, bool keep_first,
		     shared_ptr<SplineSurface>& srf2, int bd2, bool keep_second,
		     bool found_corner1, Point corner1, bool found_corner2,
		     Point corner2, bool opposite)
//===========================================================================
{
    // Make sure that the parameter directions of the two surfaces correspond.
    // This means that we only need to average coefs at start or end, i.e.
    // we swap such that the matching edges are numer 2 or 3.
    // Let the coefficients with constant v-parameter be the ones to average
    if (bd1 <= 1)
	srf1->swapParameterDirection();
    if (bd2 <= 1)
	srf2->swapParameterDirection();

    if (opposite)
	srf2->reverseParameterDirection(true);

    // Make sure that the parameter interval in u-direction is unique
    double umin1 = srf1->startparam_u();
    double umax1 = srf1->endparam_u();
    double umin2 = srf2->startparam_u();
    double umax2 = srf2->endparam_u();
    double ptol = 1.0e-12;
    if (fabs(umin1-umin2) > ptol || fabs(umax1-umax2) > ptol)
    {
	double umin = 0.5*(umin1 + umin2);
	double umax = 0.5*(umax1 + umax2);
    
	srf1->setParameterDomain(umin, umax, srf1->startparam_v(), srf1->endparam_v());
	srf2->setParameterDomain(umin, umax, srf2->startparam_v(), srf2->endparam_v());
    }

    // Ensure the same spline space in the u-direction
    vector<shared_ptr<SplineSurface> > sfs(2);
    sfs[0] = srf1;
    sfs[1] = srf2;
    GeometryTools::unifySurfaceSplineSpace(sfs, ptol, 1);
    srf1 = sfs[0];
    srf2 = sfs[1];

    // Check degeneracy
    double deg_tol = 1.0e-6;
    bool b1, r1, t1, l1, b2, r2, t2, l2;
    bool degen1 = srf1->isDegenerate(b1, r1, t1, l1, deg_tol);
    bool degen2 = srf2->isDegenerate(b2, r2, t2, l2, deg_tol);
    degen1 = (bd1 == 1 || bd1 == 3) ? t1 : b1;
    degen2 = (bd2 == 1 || bd2 == 3) ? t2 : b2;
    // Replace the specified boundary coefficients with the average ones. It is either the first
    // or the last row of coefficients
    // Be careful not to destroy degenerate boundaries
    int dim = srf1->dimension();
    vector<double>::iterator c1 = srf1->coefs_begin();
    int in1 = srf1->numCoefs_u();
    vector<double>::iterator c2 = srf2->coefs_begin();
    if (bd1 == 1 || bd1 == 3)
	c1 += (srf1->numCoefs_v()-1)*in1*dim;
    if (bd2 == 1 || bd2 == 3)
	c2 += (srf2->numCoefs_v()-1)*in1*dim;

    // Check if the corner information corresponds to the boundary orientation
    vector<double> d1(dim), d2(dim);
    vector<double>::iterator c3 = c1;
    vector<double>::iterator c4 = c1 + (in1-1)*dim;
    for (int kr=0; kr<dim; kr++, c3++, c4++)
    {
	d1[kr] = *c3;
	d2[kr] = *c4;
    }
    Point pt1(d1.begin(), d1.end());
    Point pt2(d2.begin(), d2.end());
    if ((found_corner1 && found_corner2 && pt1.dist(corner1) > pt1.dist(corner2)) ||
	(!found_corner1 && pt1.dist(corner2) < pt2.dist(corner2)) ||
	(!found_corner2 && pt2.dist(corner1) < pt1.dist(corner1)))
    {
	// Switch 
	std::swap(found_corner1,found_corner2);
	std::swap(corner1,corner2);
    }

    if (!(degen1 && degen2))
    {
	for (int ki=0; ki<in1*dim; ki++, c1++, c2++)
	{
	    bool start = (ki==0 || ki==1 || ki==2);
	    bool end = (ki==in1*dim-3 || ki==in1*dim-2 || ki==in1*dim-1);
	    if (start && l1 && l2)
		continue;
	    if (end && r1 && r2)
		continue;
	    double tmid = 0.5*((*c1) + (*c2));
	    if (start && found_corner1)
		tmid = corner1[ki%dim];
	    if (end && found_corner2)
		tmid = corner2[ki%dim];
	    if (keep_first || degen1 || (start && l1) || (end && r1))
		tmid = *c1;
	    if (keep_second || degen2 || (start && l2) || (end && r2))
		tmid = *c2;
	    *c1 = tmid;
	    *c2 = tmid;
	}
    }

    if (srf1->rational())
      {
	// This fix will not work for all configuration of rational surfaces.
	// Update the rational coefficients with respect to the divided
	// ones
	c1 = srf1->coefs_begin();
	vector<double>::iterator r1 = srf1->rcoefs_begin();
	int kn = srf1->numCoefs_u()*srf1->numCoefs_v();
	for (int ki=0; ki<kn; ++ki)
	  {
	    for (int kr=0; kr<dim; ++kr)
	      r1[kr] = c1[kr]*r1[dim];
	    c1 += dim;
	    r1 += (dim+1);
	  }
      }

    if (srf2->rational())
      {
	// This fix will not work for all configuration of rational surfaces.
	// Update the rational coefficients with respect to the divided
	// ones
	c2 = srf2->coefs_begin();
	vector<double>::iterator r2 = srf2->rcoefs_begin();
	int kn = srf2->numCoefs_u()*srf2->numCoefs_v();
	for (int ki=0; ki<kn; ++ki)
	  {
	    for (int kr=0; kr<dim; ++kr)
	      r2[kr] = c2[kr]*r2[dim];
	    c2 += dim;
	    r2 += (dim+1);
	  }
      }

     // Set the surfaces back to the initial parameter domain
    if (fabs(umin1-umin2) > ptol || fabs(umax1-umax2) > ptol)
    {
	srf1->setParameterDomain(umin1, umax1, srf1->startparam_v(), srf1->endparam_v());
	srf2->setParameterDomain(umin2, umax2, srf2->startparam_v(), srf2->endparam_v());
    }
    
    if (opposite)
	srf2->reverseParameterDirection(true);

    if (bd1 <= 1)
	srf1->swapParameterDirection();
    if (bd2 <= 1)
	srf2->swapParameterDirection();


}

//==========================================================================
void GeometryTools::setParameterDomain(vector<shared_ptr<BoundedSurface> >& sfs,
				       double u1, double u2, 
				       double v1, double v2)
//==========================================================================
{
  // Fetch information about underlying surfaces
  vector<shared_ptr<ParamSurface> > under_sfs;
  for (size_t ki=0; ki<sfs.size(); ++ki)
    {
      shared_ptr<ParamSurface> surf = sfs[ki]->underlyingSurface();
      size_t kj;
      for (kj=0; kj<under_sfs.size(); ++kj)
	if (under_sfs[kj].get() == surf.get())
	  break;
      if (kj == under_sfs.size())
	under_sfs.push_back(surf);
    }

  // Set new parameter domain for the trimming loop
  for (size_t ki=0; ki<sfs.size(); ++ki)
    sfs[ki]->setParameterDomainBdLoops(u1, u2, v1, v2);

  // Set new parameter domain for the underlying surfaces
  for (size_t ki=0; ki<under_sfs.size(); ++ki)
    under_sfs[ki]->setParameterDomain(u1, u2, v1, v2);
  
}
} // namespace Go
