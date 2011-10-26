//===========================================================================
//                                                                           
// File: SplineSurfaceInt.C 
//                                                                           
// Created: 
//                                                                           
// Author: oan
//                                                                           
// Revision: $Id: SplineSurfaceInt.C,v 1.49 2009-03-03 11:32:05 vsk Exp $
//                                                                           
// Description:
//                                                                           
//===========================================================================


#include "GoTools/intersections/SplineSurfaceInt.h"
#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/intersections/SplineCurveInt.h"
#include "GoTools/geometry/GeometryTools.h"
#include "GoTools/intersections/AlgObj3DInt.h"
#include "GoTools/implicitization/ImplicitizeSurfaceAlgo.h"
#include "GoTools/utils/RotatedBox.h"
#include "GoTools/geometry/Utils.h"
#include <fstream> // For debugging


using std::vector;
using std::ofstream;
using std::endl;
using std::min;
using std::cout;
using std::shared_ptr;
using std::dynamic_pointer_cast;


namespace Go {


//===========================================================================
SplineSurfaceInt::SplineSurfaceInt(shared_ptr<ParamSurface> surf)
    : ParamSurfaceInt(surf)
//===========================================================================
{
    spsf_ = dynamic_pointer_cast<SplineSurface, ParamSurface>(surf);

    if (parent_ && parent_->numParams() == 2);  
            // K-regularity must be ensured at an earlier recursion level
    else if (checkPeriodicity(0) >= 1 || checkPeriodicity(1) >= 1) {
	// Make surface k-regular
	spsf_ = (shared_ptr<SplineSurface>)
	    (spsf_->subSurface(spsf_->startparam_u(),
			       spsf_->startparam_v(),
			       spsf_->endparam_u(),
			       spsf_->endparam_v()));
	surf_ = (shared_ptr<ParamSurface>)(spsf_);
    }

    setImplicitDeg();
}


//===========================================================================
SplineSurfaceInt::SplineSurfaceInt(shared_ptr<ParamSurface> surf,
				   ParamGeomInt *parent)
  : ParamSurfaceInt(surf, parent)
//===========================================================================
{
    spsf_ = dynamic_pointer_cast<SplineSurface, ParamSurface>(surf);
    
    if (parent_ && parent_->numParams() == 2);  
            // K-regularity must be ensured at an earlier recursion level
    else if (checkPeriodicity(0) >= 1 || checkPeriodicity(1) >= 1) {
	// Make surface k-regular
	spsf_ = (shared_ptr<SplineSurface>)
	    (spsf_->subSurface(spsf_->startparam_u(),
			       spsf_->startparam_v(),
			       spsf_->endparam_u(),
			       spsf_->endparam_v()));
	surf_ = (shared_ptr<ParamSurface>)(spsf_);
    }

    // Pick part of normal surface if this one exists
    ParamSurfaceInt *parentsf = parent->getParamSurfaceInt();
    if (parentsf && parentsf->isSpline()) {
	SplineSurfaceInt *parentInt
	    = dynamic_cast<SplineSurfaceInt*>(parentsf);
	if (parentInt->normalsf_.get() != 0) {
	    SplineSurface *normalsf
		= parentInt->normalsf_->subSurface(spsf_->startparam_u(),
						   spsf_->startparam_v(),
						   spsf_->endparam_u(),
						   spsf_->endparam_v());
	    normalsf_ = shared_ptr<SplineSurface>(normalsf);
	}
    }

    setImplicitDeg();  
}


//===========================================================================
std::shared_ptr<ParamSurfaceInt> 
SplineSurfaceInt::makeIntObject(shared_ptr<ParamSurface> surf)
//===========================================================================
{
    shared_ptr<SplineSurfaceInt> surf_int =
	shared_ptr<SplineSurfaceInt>(new SplineSurfaceInt(surf, this));
    return surf_int;
}


//===========================================================================
std::shared_ptr<ParamCurveInt> 
SplineSurfaceInt::makeIntCurve(shared_ptr<ParamCurve> crv, 
			       ParamGeomInt *parent)
//===========================================================================
{
    shared_ptr<SplineCurveInt> curve_int =
	shared_ptr<SplineCurveInt>(new SplineCurveInt(crv, parent));
    return curve_int;
}


//===========================================================================
int SplineSurfaceInt::checkPeriodicity(int pardir) const
//===========================================================================
{
    ASSERT(pardir == 0 || pardir == 1);
    int per = analyzePeriodicity(*(spsf_.get()), pardir);  
    return per;
}


//===========================================================================
void SplineSurfaceInt::
getBoundaryObjects(std::vector<std::shared_ptr<BoundaryGeomInt> >& bd_objs)
//===========================================================================
{
    // Compute number of boundaries
    int nmb_bd = 0;
    int per1 = -1; //analyzePeriodicity(*(spsf_.get()), 0);  
    int per2 = -1; //analyzePeriodicity(*(spsf_.get()), 1);
    if (per1 < 1)
	nmb_bd += 2;
    if (per2 < 1)
	nmb_bd += 2;
    if (int(boundary_obj_.size()) == nmb_bd) {
	bd_objs.insert(bd_objs.begin(), boundary_obj_.begin(), 
		       boundary_obj_.end());
    } else {
	boundary_obj_.clear();

	// First parameter direction
	int ki;
	double par[2];
	if (per1 < 1) {
	    par[0] = spsf_->startparam_u();
	    par[1] = spsf_->endparam_u();
	    for (ki=0; ki<2; ki++) {
		shared_ptr<ParamCurve> curr_cv = 
		    shared_ptr<ParamCurve>
		    (spsf_->constParamCurve(par[ki], false));
		shared_ptr<ParamCurveInt> curr_cv_int = 
		    shared_ptr<ParamCurveInt>
		    (new SplineCurveInt(curr_cv, this));
		shared_ptr<BoundaryGeomInt> bd = 
		    shared_ptr<BoundaryGeomInt>
		    (new BoundaryGeomInt(curr_cv_int, 0, par[ki]));
		boundary_obj_.push_back(bd);
		bd_objs.push_back(bd);
	    }
	}
  
	// Second parameter direction
	if (per2 < 1) {
	    par[0] = spsf_->startparam_v();
	    par[1] = spsf_->endparam_v();
	    for (ki=0; ki<2; ki++) {
		shared_ptr<ParamCurve> curr_cv = 
		    shared_ptr<ParamCurve>
		    (spsf_->constParamCurve(par[ki], true));
		shared_ptr<ParamCurveInt> curr_cv_int = 
		    shared_ptr<ParamCurveInt>
		    (new SplineCurveInt(curr_cv, this));
		shared_ptr<BoundaryGeomInt> bd = 
		  shared_ptr<BoundaryGeomInt>
		    (new BoundaryGeomInt(curr_cv_int, 1, par[ki]));
	      boundary_obj_.push_back(bd);
	      bd_objs.push_back(bd);
	    }
	}
    }
  
}


//===========================================================================
int SplineSurfaceInt::getMeshSize(int dir)
//===========================================================================
{
    int meshsize = 5;
    if (spsf_->numCoefs_u() < meshsize || spsf_->numCoefs_v()< meshsize)
	return ParamSurfaceInt::getMeshSize(dir);
    else if (dir == 0)
	return spsf_->numCoefs_u();
    else if (dir == 1)
	return spsf_->numCoefs_v();
    else
	return 1;
}


//===========================================================================
bool SplineSurfaceInt::hasInnerKnots(int pardir) const
//===========================================================================
{
    ASSERT(pardir == 0 || pardir == 1);
    return (pardir == 0) ? (spsf_->numCoefs_u() > spsf_->order_u())
	: (spsf_->numCoefs_v() > spsf_->order_v());
}


//===========================================================================
bool SplineSurfaceInt::hasCriticalValsOrKnots(int pardir) const 
//===========================================================================
{
    ASSERT(pardir == 0 || pardir == 1);
    bool critical = hasCriticalVals(pardir);
    return (segment_[pardir].size() > 0 || critical);
}


//===========================================================================
vector<double> SplineSurfaceInt::getCriticalValsAndKnots(int pardir) const
//===========================================================================
{
    vector<double> critical = getCriticalVals(pardir);
    vector<double> knots = getInnerKnotVals(pardir, false);
    vector<double> vals;
    double ta = startParam(pardir);
    double tb = endParam(pardir);
    vals.push_back(ta);
    int ki, kj;
    for (ki=0, kj=0; ;)
    {
	if (ki >= int(critical.size()) && kj >= int(knots.size()))
	    break;
	else if (ki >= int(critical.size()))
	{
	    if (knots[kj] > vals[vals.size()-1])
		vals.push_back(knots[kj]);
	    kj++;
	}
	else if (kj >= int(knots.size()))
	{
	    if (critical[ki] > vals[vals.size()-1])
		vals.push_back(critical[ki]);
	    ki++;
	}
	else if (critical[ki] < knots[kj])
	{
	    if (critical[ki] > vals[vals.size()-1])
		vals.push_back(critical[ki]);
	    ki++;
	}
	else
	{
	    if (critical[ki] > vals[vals.size()-1])
		vals.push_back(critical[ki]);
	    ki++;
	}
    }
    if (tb > vals[vals.size()-1])
	vals.push_back(tb);

    return vals;
}


//===========================================================================
struct sort_distance
//===========================================================================
{
    // Functor used to sort elements of a vector in getInnerKnotVals()

    double mid;
    sort_distance(double start, double end)
    { mid = 0.5*(start+end); }

    bool operator()(double a, double b) const
    {
	return (fabs(a-mid) < fabs(b-mid));
    }

};


//===========================================================================
vector<double> SplineSurfaceInt::
getInnerKnotVals(int pardir, bool sort) const
//===========================================================================
{
    ASSERT(pardir == 0 || pardir == 1);

    // First fetch all inner knots
    vector<double> vals;
    int kk = (pardir == 0) ? spsf_->order_u() : spsf_->order_v();
    int kn = (pardir == 0) ? spsf_->numCoefs_u() : spsf_->numCoefs_v();
    if (kk == kn)
	return vals;

    std::vector<double>::const_iterator et = spsf_->basis(pardir).begin();
    vals.push_back(et[kk]);
    for (int ki = kk+1; ki < kn; ki++) {
	if (et[ki] > vals[vals.size()-1]) {
	    vals.push_back(et[ki]);
	}
    }
    if (sort) {
	// Sort knot vector with respect to the distance from the
	// midpoint of the current parameter interval
	sort_distance compare(startParam(pardir), endParam(pardir));
	std::sort(vals.begin(), vals.end(), compare);
    }

    return vals;
}


//===========================================================================
vector<double>::iterator SplineSurfaceInt::getMesh()
//===========================================================================
{
    int meshsize = 5;
    if (spsf_->numCoefs_u() < meshsize || spsf_->numCoefs_v()< meshsize)
	return ParamSurfaceInt::getMesh();
    else
	return spsf_->coefs_begin();
}


//===========================================================================
double SplineSurfaceInt::paramFromMesh(int dir, int idx)
//===========================================================================
{
    int meshsize = 5;
    if (spsf_->numCoefs_u() < meshsize || spsf_->numCoefs_v()< meshsize) {
	return ParamSurfaceInt::paramFromMesh(dir, idx);
    } else {
	return (dir == 0 || dir == 1) ? 
	    spsf_->basis(dir).grevilleParameter(idx) : 0.0;
    }
}


//===========================================================================
shared_ptr<ParamSurfaceInt> SplineSurfaceInt::getNormalSurface() const
//===========================================================================
{
    if (normalsf_.get() == 0)
	normalsf_ = (shared_ptr<SplineSurface>)(spsf_->normalSurface());

    return (shared_ptr<ParamSurfaceInt>)(new SplineSurfaceInt(normalsf_));
}


//===========================================================================
DirectionCone SplineSurfaceInt::reducedDirectionCone(bool reduce_at_bd[4],
						     double epsge) const
//===========================================================================
{
    // Make exact direction cone. Only for Bezier case
    if (spsf_->numCoefs_u() > 3*spsf_->order_u() || 
	spsf_->numCoefs_v() > 3*spsf_->order_v())
    {
	return ParamSurfaceInt::directionCone();
    }

    // Make sure that a normal surface is computed
    if (normalsf_.get() == 0)
	normalsf_ = (shared_ptr<SplineSurface>)(spsf_->normalSurface());

    // Find parameter intervals of reduced normal surface
    double param[4];
    for (int ki=0; ki<2; ki++)
    {
	if (reduce_at_bd[2*ki])
	    param[2*ki] = getParOffBd(ki, true, epsge);
	else
	    param[2*ki] = startParam(ki);

	if (reduce_at_bd[2*ki+1])
	    param[2*ki+1] = getParOffBd(ki, false, epsge);
	else
	    param[2*ki+1] = endParam(ki);
    }
	
    DirectionCone cone2;
    // Make reduced normal surface
    if (param[1] <= param[0] || param[3] <= param[2])
	return cone2;  // Too small surface piece. Makes no sense

    shared_ptr<SplineSurface> red_sf =
	(shared_ptr<SplineSurface>)(normalsf_->subSurface(param[0], param[2],
							  param[1], param[3]));

    // Make cone from reduced normal surface
    vector<double>::iterator coefs_start = red_sf->coefs_begin();
    vector<double>::iterator coefs_end = red_sf->coefs_end();
    try {
        int nmb_elem = (int)(coefs_end - coefs_start);
	cone2.setFromArray(&coefs_start[0], &coefs_start[0] + nmb_elem, dim_);
    } catch (...)
    {
	cone2 = ParamSurfaceInt::directionCone();
    }

     return cone2;
}


//===========================================================================
DirectionCone SplineSurfaceInt::directionCone() const
//===========================================================================
{
    // Make exact direction cone. Only for Bezier case
    if (spsf_->numCoefs_u() > 3*spsf_->order_u() || 
	spsf_->numCoefs_v() > 3*spsf_->order_v())
    {
	return ParamSurfaceInt::directionCone();
    }

    DirectionCone cone2 = spsf_->normalCone();
    if (cone_.greaterThanPi() < 0 || normalsf_.get() == 0)
    {
	// Make sure that a normal surface is computed
	if (normalsf_.get() == 0)
	    normalsf_ = (shared_ptr<SplineSurface>)(spsf_->normalSurface());

	// Make cone
	vector<double>::iterator coefs_start = normalsf_->coefs_begin();
	vector<double>::iterator coefs_end = normalsf_->coefs_end();
	int nmb_elem = (int)(coefs_end - coefs_start);
	try {
	cone_.setFromArray(&coefs_start[0], &coefs_start[0] + nmb_elem, dim_); //0], dim_);
	} catch (...)
	{
	    cone_ = cone2;
	}

	if (cone_.angle() > cone2.angle() ||
	    (cone_.greaterThanPi() && !cone2.greaterThanPi()))
	{
	    ofstream debug("cone_sf.g2");
	    spsf_->writeStandardHeader(debug);
	    spsf_->write(debug);
	    ofstream debug2("normal_sf.g2");
	    normalsf_->writeStandardHeader(debug2);
	    normalsf_->write(debug2);
	    cone_ = cone2;
	}
    }
    return cone_;
}


//===========================================================================
bool SplineSurfaceInt::isSimple()
//===========================================================================
{
    // Estimates if the current surface is simple enough for a singularity
    // iteration. Checks the span of the normal cone and the size of the
    // surface
    int fac1 = 5;
    int fac2 = 20;

    if (spsf_->numCoefs_u() == spsf_->order_u() &&
	spsf_->numCoefs_v() == spsf_->order_v())
	return true;           // Bezier case

    if (spsf_->numCoefs_u() > fac1*spsf_->order_u() &&
	spsf_->numCoefs_v() > fac1*spsf_->order_v())  
	return false;

    if (spsf_->numCoefs_u() > fac2*spsf_->order_u() ||
	spsf_->numCoefs_v() > fac2*spsf_->order_v())  
	return false;

    return ParamSurfaceInt::isSimple();
}


//===========================================================================
bool SplineSurfaceInt::isSpline()
//===========================================================================
{
    return true;
}


//===========================================================================
bool SplineSurfaceInt::isIsoParametric(ParamCurveInt *curve, int dir, double par,
				       double ptol, double tol)
//===========================================================================
{
    if (!curve->isSpline())
	return false;  // Don't know

    const SplineCurve *spcv = curve->getSpline();

    // Check order and number of coefficients
    int sf_order = (dir == 1) ? spsf_->order_u() : spsf_->order_v();
    int sf_ncoef = (dir == 1) ? spsf_->numCoefs_u() : spsf_->numCoefs_v();
    if (spcv->order() != sf_order || spcv->numCoefs() != sf_ncoef)
	return false;

    // Check if none or both objects are rational
    if ((spcv->rational() && !(spsf_->rational())) ||
	(spsf_->rational() && !(spcv->rational())))
	return false;

    // Check knot vector
    vector<double>::const_iterator cv_knots = spcv->knotsBegin();
    vector<double>::const_iterator sf_knots = spsf_->basis(1-dir).begin();
    for (int ki=0; ki<sf_order+sf_ncoef; ki++)
	if (fabs(cv_knots[ki]-sf_knots[ki]) > ptol)
	    return false;

    // Check coefficients
    shared_ptr<SplineCurve> sf_cv(spsf_->constParamCurve(par, (dir == 1)));
    vector<double>::const_iterator cv_coefs = (spcv->rational()) ? spcv->rcoefs_begin() :
	spcv->coefs_begin();
    vector<double>::const_iterator sf_coefs = (spcv->rational()) ? sf_cv->rcoefs_begin() :
	sf_cv->coefs_begin();
    int dim = spcv->dimension();
    if (spcv->rational())
	dim++;

    double tol2 = tol*tol;
    for (int ki=0; ki<sf_order; ki++)
    {
	double d2 = distance_squared(&(cv_coefs+ki*dim)[0], 
				     &(cv_coefs+(ki+1)*dim)[0],
				     &(sf_coefs+ki*dim)[0]);
	if (d2 > tol2)
	    return false;
    }

    return true;
}


//===========================================================================
double SplineSurfaceInt::getOptimizedConeAngle(Point& axis1, Point& axis2)
//===========================================================================
{
    // We are making a cone surrounding the orientating surface on the
    // unit sphere. The cone is representated with centre coordinates
    // and an angle. The orientation is computed from aproximation of
    // the normal to the surface.  Based on the sisl function s1795.

    // Use exact surface normals if a normal surface exist
    double angle = 0.0;
    if (normalsf_.get() == 0)
    {
	int in1 = spsf_->numCoefs_u();
	int in2 = spsf_->numCoefs_v();
	vector<double>::iterator coefs = spsf_->coefs_begin();

	Point corner[4];  // The coefficients making the corner of
			  // each patch
	Point diff[4];    // Difference vector between corner
			  // coefficients
	Point norm;       // Estimated surface normal (cross product
	                  // between difference vectors)
	int kver, khor;   // The index to the vertice in the upper
	                  // left corner to the patch to treat.
	vector<double>::iterator it1;
	int ki, kj;

	// Here we are treating each patch in the control polygon
	// separately.
	for (it1=coefs, kver=0; kver < (in2-1); kver++, it1+=dim_)
	    for (khor=0; khor < (in1-1); khor++, it1+=dim_)
	    {
		// Here we make the tangents in each corner of the
		// patch, and in direction with the clock. The first
		// and the last vector contains both the first
		// tangent.
		corner[0].resize(dim_);
		corner[0].setValue(&it1[0]);
		corner[1].resize(dim_);
		corner[1].setValue(&it1[dim_]);
		corner[2].resize(dim_);
		corner[2].setValue(&it1[(in1+1)*dim_]);
		corner[3].resize(dim_);
		corner[3].setValue(&it1[in1+dim_]);
		for (ki=0; ki<4; ki++)
		{
		    kj = ((ki+1) % 4);
		    diff[ki] = corner[kj] - corner[ki];
		}
	
		// Here we makes the normales in each corner of the
		// patch.  We are using a cross product between two
		// tangents.  The normals ar also normalized.
		for (ki=0; ki<4; ki++)
		{
		    kj = (ki == 0) ? 3 : ki-1;
		    norm = diff[kj].cross(diff[ki]);
		    double len = norm.normalize_checked();
		    if (len == 0.0)
			norm = axis1;

		    double t1 = axis1*norm;
		    double t2 = axis2*norm;
		    double tang = t1/sqrt(t1*t1 + t2*t2);
		    tang = std::min(tang, 1.0);
		    tang = std::max(tang, -1.0);
		    tang = acos(tang);
		    angle = std::max(angle, tang);
		}
	    }
    }
    else
    {
	//int in1 = normalsf_->numCoefs_u();
	//int in2 = normalsf_->numCoefs_v();
	vector<double>::iterator coefs_start = normalsf_->coefs_begin();
	vector<double>::iterator coefs_end = normalsf_->coefs_end();
	vector<double>::iterator it1;
	Point norm(dim_);       // Estimated surface normal (cross
				// product between
	for (it1=coefs_start; it1<coefs_end; it1+=dim_)
	{
	    norm.setValue(&it1[0]);
	    double len = norm.normalize_checked();
	    if (len == 0.0)
		norm = axis1;

	    double t1 = axis1*norm;
	    double t2 = axis2*norm;
	    double tang = t1/sqrt(t1*t1 + t2*t2);
	    tang = std::min(tang, 1.0);
	    tang = std::max(tang, -1.0);
	    tang = acos(tang);
	    angle = std::max(angle, tang);
	}
	    

    }	
    return angle;
}


//===========================================================================
RotatedBox  SplineSurfaceInt::getRotatedBox(std::vector<Point>& axis) const
//===========================================================================
{
    RotatedBox box(spsf_->coefs_begin(), dimension(), spsf_->numCoefs_u(),
		   spsf_->numCoefs_v(), &axis[0]);
    return box;
}


//===========================================================================
void SplineSurfaceInt::
knotIntervalFuzzy(double& u, double&v, double utol, double vtol) const
//===========================================================================
{
    spsf_->basis_u().knotIntervalFuzzy(u, utol);
    spsf_->basis_v().knotIntervalFuzzy(v, vtol);
}


//===========================================================================
double SplineSurfaceInt::
nextSegmentVal(int dir, double par, bool forward, double tol) const
//===========================================================================
{
    return spsf_->nextSegmentVal(dir, par, forward, tol);
}


//===========================================================================
void SplineSurfaceInt::
splitAtG0(double angtol,
	  std::vector<std::shared_ptr<ParamSurfaceInt> >& subG1)
//===========================================================================
{
    // Split surface at G1 discontinuities

    vector<double> g1_disc_u;
    vector<double> g1_disc_v;
    surfaceKinks(*spsf_, angtol, g1_disc_u, g1_disc_v);

    vector<shared_ptr<SplineSurface> > subspline;
    subspline = splitInKinks(*spsf_, g1_disc_u, g1_disc_v);

    for (size_t ki=0; ki<subspline.size(); ki++)
    {
	subG1.push_back(shared_ptr<ParamSurfaceInt>
			(new SplineSurfaceInt(subspline[ki])));
	subG1[ki]->setParent(this);
    }
}


//===========================================================================
bool SplineSurfaceInt::canImplicitize()
//===========================================================================
{
    // There are at present no restrictions on the possibility to
    // implicitize: Rational is OK. Not Bezier is OK.

    return true;
}


//===========================================================================
bool SplineSurfaceInt::implicitize(double tol)
//===========================================================================
{
    // NOTE: 'tol' is not in use at the moment. Neither is
    // 'implicit_tol_'. @@@jbt

    // Check to see if we can implicitize
    bool can_impl = canImplicitize();
    if (!can_impl) {
	return can_impl;
    }

    // Initialize with implicit degree = 1
    int deg = 1;
    impl_sf_algo_ = shared_ptr<ImplicitizeSurfaceAlgo>
	(new ImplicitizeSurfaceAlgo(*spsf_, deg));
    impl_sf_algo_->setTolerance(tol); // At present no effect. @@@jbt
    impl_sf_algo_->perform();

    // Check if we have a plane. exact_eps is the tolerance we check
    // against to see if we have exact implicitization.
    const double exact_eps = 1.0e-12;
    BernsteinTetrahedralPoly impl;
    BaryCoordSystem3D bc;
    impl_sf_algo_->getResultData(impl, bc, implicit_err_);
    if (implicit_err_ <= exact_eps) {
	if (getenv("DEBUG_IMPL") && (*getenv("DEBUG_IMPL"))=='1') {
	cout << "Implicit degree = 1" << endl;
	}
	implicit_obj_ = shared_ptr<AlgObj3DInt>(new AlgObj3DInt(impl, bc));
	return true;
    }

    // Not a plane - we need to loop through a series of
    // degrees. maxdeg is our choice of maximal degree.
    const int maxdeg = 4;
    for (int ideg = 2; ideg <= maxdeg; ++ideg) {
	impl_sf_algo_->setDegree(ideg);
	impl_sf_algo_->perform();
	impl_sf_algo_->getResultData(impl, bc, implicit_err_);
	if (implicit_err_ <= exact_eps) {
	    break;
	}
    }
    if (getenv("DEBUG_IMPL") && (*getenv("DEBUG_IMPL"))=='1') {
	cout << "Implicit degree = " << impl.degree()
	     << "  Implicit error = " << implicit_err_ << endl;
	ofstream os("implicit.dat");
	os << impl << endl
	   << bc << endl;
    }
    implicit_obj_ = shared_ptr<AlgObj3DInt>(new AlgObj3DInt(impl, bc));

    return true;
}


//===========================================================================
bool SplineSurfaceInt::getImplicit(double tol, double& tol2,
				   AlgObj3DInt& alg_obj_3d_int)
//===========================================================================
{
    // NOTE: 'tol' is not in use at the moment. Neither is
    // 'implicit_tol_'. @@@jbt

    // If we haven't implicitized, return 'false'.
    if (implicit_obj_.get() == 0) {
	return false;
    }

    tol2 = implicit_err_;
    alg_obj_3d_int = *implicit_obj_;

    return true;
}


//===========================================================================
SplineCurve* SplineSurfaceInt::constParamCurve(double parameter,
					       bool pardir_is_u) const
//===========================================================================
{
    return spsf_->constParamCurve(parameter, pardir_is_u);
}


//===========================================================================
void SplineSurfaceInt::setImplicitDeg()
//===========================================================================
{
    // The order of sf may be lower than default order.
    int deg_u = spsf_->order_u() - 1;
    int deg_v = spsf_->order_v() - 1;
    int exact_alg_deg = 2*deg_u*deg_v;
    if (deg_u == 1 && deg_v == 1) {
      // The sf need not be planar. We compute the normalcone of the
      // spsf_.
      DirectionCone normal_cone;
      bool failed = false;
      try {
	normal_cone = spsf_->normalCone();
      }
      catch (...)
	{
	  failed = true; // Degenerate surface. Stay with initial algebraic degree
	}
      if (!failed)
	{
	  double ang_tol = 1e-10;
	  if (normal_cone.angle() < ang_tol) {
	    exact_alg_deg = 1;
	  }
	}
    }

    // If the exact alg degree is lower than the default we alter the
    // degree.
    impl_deg_ = min(impl_deg_, exact_alg_deg);
}


//===========================================================================


} // namespace Go
