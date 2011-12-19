//==========================================================================
//                                                                          
// File: SurfaceOfRevolution.C                                              
//                                                                          
// Created: Wed Oct 21 17:49:45 2009                                         
//                                                                          
// Author: Jan B. Thomassen <jbt@sintef.no>
//                                                                          
// Revision: $Id$
//                                                                          
// Description:
//                                                                          
//==========================================================================


#include "GoTools/geometry/SurfaceOfRevolution.h"
#include "GoTools/geometry/SplineCurve.h"
#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/geometry/GeometryTools.h"
#include "GoTools/geometry/SweepSurfaceCreator.h"
#include <vector>
#include <limits>


using std::vector;
using std::cout;
using std::endl;


namespace Go
{


// Constructor
//===========================================================================
SurfaceOfRevolution::SurfaceOfRevolution(Point location, Point axis_dir,
					 shared_ptr<SplineCurve> curve)
    : location_(location), axis_dir_(axis_dir), curve_(curve)
//===========================================================================
{
    if (location_.dimension() != 3) {
	THROW("Dimension must be 3.");
	return;
    }

    axis_dir_.normalize();
    setDefaultDomain();
}


// Destructor
//===========================================================================
SurfaceOfRevolution::~SurfaceOfRevolution()
//===========================================================================
{
}


//===========================================================================
void SurfaceOfRevolution::read (std::istream& is)
//===========================================================================
{
    bool is_good = is.good();
    if (!is_good) {
	THROW("Invalid geometry file!");
    }

    int dim;
    is >> dim;
    if (dim != 3)
	THROW("Dimension must be 3.");
    location_.resize(dim);
    axis_dir_.resize(dim);
    is >> location_
       >> axis_dir_;
    curve_ = shared_ptr<SplineCurve>(new SplineCurve);
    curve_->read(is);

    axis_dir_.normalize();
    setDefaultDomain();
}


//===========================================================================
void SurfaceOfRevolution::write(std::ostream& os) const
//===========================================================================
{
    os << dimension() << endl
       << location_ << endl
       << axis_dir_ << endl;
    curve_->write(os);
}


//===========================================================================
int SurfaceOfRevolution::dimension() const
//===========================================================================
{
    return location_.dimension();
}


//===========================================================================
ClassType SurfaceOfRevolution::instanceType() const
//===========================================================================
{
    return classType();
}


//===========================================================================
ClassType SurfaceOfRevolution::classType()
//===========================================================================
{
    return Class_SurfaceOfRevolution;
}


//===========================================================================
BoundingBox SurfaceOfRevolution::boundingBox() const
//===========================================================================
{
    // A rather unefficient hack...
    SurfaceOfRevolution* sor = const_cast<SurfaceOfRevolution*>(this);
    SplineSurface* tmp = sor->geometrySurface();
    BoundingBox box = tmp->boundingBox();
    delete tmp;
    return box;
}

//===========================================================================
SurfaceOfRevolution* SurfaceOfRevolution::clone() const
//===========================================================================
{
    return new SurfaceOfRevolution(location_, axis_dir_, curve_);
}


//===========================================================================
const RectDomain& SurfaceOfRevolution::parameterDomain() const
//===========================================================================
{
    return domain_;
}


//===========================================================================
RectDomain SurfaceOfRevolution::containingDomain() const
//===========================================================================
{
    return parameterDomain();
}


//===========================================================================
bool SurfaceOfRevolution::inDomain(double u, double v) const
//===========================================================================
{
    Array<double, 2> pt(u, v);
    // Using an arbitrary tolerance... @jbt
    double tol = 1.0e-12;
    return domain_.isInDomain(pt, tol);
}


//===========================================================================
Point SurfaceOfRevolution::closestInDomain(double u, double v) const
//===========================================================================
{
    Array<double, 2> pt(u, v);
    Array<double, 2> clo_pt;
    // Using an arbitrary tolerance... @jbt
    double tol = 1.0e-12;
    domain_.closestInDomain(pt, clo_pt, tol);
    return Point(clo_pt[0], clo_pt[1]);
}


//===========================================================================
CurveLoop
SurfaceOfRevolution::outerBoundaryLoop(double degenerate_epsilon) const
//===========================================================================
{
    MESSAGE("Does not make sense. Returns an empty loop.");
    CurveLoop loop;
    return loop;
}


//===========================================================================
std::vector<CurveLoop> 
SurfaceOfRevolution::allBoundaryLoops(double degenerate_epsilon) const
//===========================================================================
{
    MESSAGE("Does not make sense. Returns an empty vector.");
    vector<CurveLoop> loops;
    return loops;
}


//===========================================================================
DirectionCone SurfaceOfRevolution::normalCone() const
//===========================================================================
{
    // A rather unefficient hack...
    SurfaceOfRevolution* sof = const_cast<SurfaceOfRevolution*>(this);
    SplineSurface* tmp = sof->geometrySurface();
    DirectionCone dc = tmp->normalCone();
    delete tmp;
    return dc;
}


//===========================================================================
DirectionCone SurfaceOfRevolution::tangentCone(bool pardir_is_u) const
//===========================================================================
{
    // A rather unefficient hack...
    SurfaceOfRevolution* sof = const_cast<SurfaceOfRevolution*>(this);
    SplineSurface* tmp = sof->geometrySurface();
    DirectionCone dc = tmp->tangentCone(pardir_is_u);
    delete tmp;
    return dc;
}


//===========================================================================
void SurfaceOfRevolution::point(Point& pt, double upar, double vpar) const
//===========================================================================
{
    Point cvpt;
    curve_->point(cvpt, vpar);
    Point vec = cvpt - location_;

    pt = location_
	+ vec * cos(upar) + (vec * location_) * location_ * (1.0 - cos(upar))
	+ location_.cross(vec) * sin(upar);
}


//===========================================================================
void SurfaceOfRevolution::point(std::vector<Point>& pts, 
				double upar, double vpar,
				int derivs,
				bool u_from_right,
				bool v_from_right,
				double resolution) const
//===========================================================================
{
    DEBUG_ERROR_IF(derivs < 0,
		   "Negative number of derivatives makes no sense.");
    int totpts = (derivs + 1)*(derivs + 2)/2;
    int ptsz = (int)pts.size();
    DEBUG_ERROR_IF(ptsz< totpts,
		   "The vector of points must have sufficient size.");

    int dim = dimension();
    for (int i = 0; i < totpts; ++i) {
        if (pts[i].dimension() != dim) {
            pts[i].resize(dim);
        }
	pts[i].setValue(0.0);
    }

    // Zero'th derivative
    point(pts[0], upar, vpar);
    if (derivs == 0)
        return;

    // First derivatives
    double cosu = cos(upar);
    double sinu = sin(upar);
    // double cosv = cos(vpar);
    // double sinv = sin(vpar);
    vector<Point> cvpts;
    curve_->point(cvpts, vpar, 1);
    Point lc = cvpts[0] - location_;
    Point& dl = cvpts[1];
    pts[1] = -lc * sinu + (lc * location_) * location_ * sinu
	+ location_.cross(lc) * cosu;
    pts[2] = dl * cosu + (dl * location_) * location_ * (1.0 - cosu)
	+ location_.cross(dl) * sinu;
    if (derivs == 1)
	return;

    // Second order and higher derivatives.
    MESSAGE("Second order or higher derivatives not yet implemented.");

}


//===========================================================================
void SurfaceOfRevolution::normal(Point& n, double upar, double vpar) const
//===========================================================================
{
    vector<Point> pts;
    point(pts, upar, vpar, 1);
    n = pts[1].cross(pts[2]);
    n.normalize();
}


//===========================================================================
vector<shared_ptr<ParamCurve> >
SurfaceOfRevolution::constParamCurves(double parameter, bool pardir_is_u) const
//===========================================================================
{
    MESSAGE("constParamCurves() not yet implemented");
    vector<shared_ptr<ParamCurve> > res;
    return res;
}


//===========================================================================
SurfaceOfRevolution*
SurfaceOfRevolution::subSurface(double from_upar, double from_vpar,
				double to_upar, double to_vpar,
				double fuzzy) const
//===========================================================================
{
    SurfaceOfRevolution* sor = clone();
    sor->setParameterBounds(from_upar, from_vpar, to_upar, to_vpar);
    return sor;
}


//===========================================================================
vector<shared_ptr<ParamSurface> >
SurfaceOfRevolution::subSurfaces(double from_upar, double from_vpar,
				 double to_upar, double to_vpar,
				 double fuzzy) const
//===========================================================================
{
    vector<shared_ptr<ParamSurface> > res;
    shared_ptr<SurfaceOfRevolution> sor(subSurface(from_upar, from_vpar,
						   to_upar, to_vpar));
    res.push_back(sor);
    return res;
}


//===========================================================================
double 
SurfaceOfRevolution::nextSegmentVal(int dir, double par, bool forward,
				    double tol) const
//===========================================================================
{
    MESSAGE("Does not make sense. Return arbitrarily zero.");
    return 0.0;
}


//===========================================================================
void SurfaceOfRevolution::closestPoint(const Point& pt,
				       double&        clo_u,
				       double&        clo_v, 
				       Point&         clo_pt,
				       double&        clo_dist,
				       double         epsilon,
				       const RectDomain* domain_of_interest,
				       double   *seed) const
//===========================================================================
{
    MESSAGE("May provide incorrect result - use with caution!");

    // We use the closest point function from SplineSurface. This will
    // in general lead to incorrect results!

    // Get the SplineSurface
    SplineSurface* sf = geometrySurface();

    // We first use the given seed to find values
    RectDomain curr_domain_of_interest = domain_;
    if (domain_of_interest != NULL) {
	curr_domain_of_interest.intersectWith(*domain_of_interest);
    }
    double curr_seed[2];
    if (seed == NULL) {
	curr_seed[0] = 0.5 * (curr_domain_of_interest.umin() 
			      + curr_domain_of_interest.umax());
	curr_seed[1] = 0.5 * (curr_domain_of_interest.vmin() 
			      + curr_domain_of_interest.vmax());
    }
    else {
	curr_seed[0] = seed[0];
	curr_seed[1] = seed[1];
    }
    sf->closestPoint(pt, clo_u, clo_v, clo_pt, clo_dist, epsilon,
		     &curr_domain_of_interest, curr_seed);

    // Try to reseed
    double seeds[4][2];
    seeds[0][0] = curr_domain_of_interest.umin();
    seeds[0][1] = curr_domain_of_interest.vmin();
    seeds[1][0] = curr_domain_of_interest.umax();
    seeds[1][1] = curr_domain_of_interest.vmin();
    seeds[2][0] = curr_domain_of_interest.umin();
    seeds[2][1] = curr_domain_of_interest.vmax();
    seeds[3][0] = curr_domain_of_interest.umax();
    seeds[3][1] = curr_domain_of_interest.vmax();
    Point tmppt(3);
    double tmpu, tmpv, tmpdist;
    for (int i = 0; i < 4; ++i) {
	sf->closestPoint(pt, tmpu, tmpv, tmppt, tmpdist, epsilon,
			 &curr_domain_of_interest, seeds[i]);
	if (tmpdist < clo_dist - epsilon) {
	    clo_u = tmpu;
	    clo_v = tmpv;
	    clo_pt = tmppt;
	    clo_dist = tmpdist;
	}
    }

    delete sf;
}


//===========================================================================
void SurfaceOfRevolution::closestBoundaryPoint(const Point& pt,
					       double&        clo_u,
					       double&        clo_v, 
					       Point&       clo_pt,
					       double&        clo_dist,
					       double epsilon,
					       const RectDomain* rd,
					       double *seed) const
//===========================================================================
{
    MESSAGE("May provide incorrect result - use with caution!");

    // This is a bit like cheating...

    SplineSurface* sf = geometrySurface();
    sf->closestBoundaryPoint(pt, clo_u, clo_v, clo_pt, clo_dist, epsilon,
			     rd, seed);
    delete sf;
}


//===========================================================================
void
SurfaceOfRevolution::getBoundaryInfo(Point& pt1, Point& pt2,
				     double epsilon, SplineCurve*& cv,
				     SplineCurve*& crosscv,
				     double knot_tol) const
//===========================================================================
{
    MESSAGE("getBoundaryInfo() not yet implemented");
}


//===========================================================================
void SurfaceOfRevolution::turnOrientation()
//===========================================================================
{
    swapParameterDirection();
}


//===========================================================================
void SurfaceOfRevolution::swapParameterDirection()
//===========================================================================
{
    MESSAGE("swapParameterDirection() not implemented.");
}


//===========================================================================
double SurfaceOfRevolution::area(double tol) const
//===========================================================================
{
    MESSAGE("area() not implemented.");
    return -1.0;
}


//===========================================================================
void SurfaceOfRevolution::reverseParameterDirection(bool direction_is_u)
//===========================================================================
{
    MESSAGE("reverseParameterDirection() not implemented.");
}


//===========================================================================
bool SurfaceOfRevolution::isDegenerate(bool& b, bool& r,
				       bool& t, bool& l,
				       double tolerance) const
//===========================================================================
{
    MESSAGE("isDegenerate() not implemented.");
    return false;
}


//===========================================================================
void SurfaceOfRevolution::getDegenerateCorners(vector<Point>& deg_corners,
					       double tol) const
//===========================================================================
{
    MESSAGE("getDegenerateCorners() not implemented.");
}


//===========================================================================
  void SurfaceOfRevolution::getCornerPoints(std::vector<std::pair<Point,Point> >& corners) const
//===========================================================================
{
    MESSAGE("getCornerPoints() not implemented.");
}

//===========================================================================
void SurfaceOfRevolution::setDefaultDomain()
//===========================================================================
{
    double umin = 0.0;
    double umax = 2.0 * M_PI;
    double vmin = curve_->startparam();
    double vmax = curve_->endparam();
    Array<double, 2> ll(umin, vmin);
    Array<double, 2> ur(umax, vmax);
    domain_ = RectDomain(ll, ur);
}


//===========================================================================
void
SurfaceOfRevolution::setParameterBounds(double from_upar, double from_vpar,
					double to_upar, double to_vpar)
//===========================================================================
{
    if (from_upar >= to_upar )
	THROW("First u-parameter must be strictly less than second.");
    if (from_upar < -2.0 * M_PI || to_upar > 2.0 * M_PI)
	THROW("u-arameters must be in [-2pi, 2pi].");
    if (to_upar - from_upar > 2.0 * M_PI)
	THROW("(to_upar - from_upar) must not exceed 2pi.");

//     curve_->basis().rescale(from_vpar, to_vpar);

    Array<double, 2> ll(from_upar, from_vpar);
    Array<double, 2> ur(to_upar, to_vpar);
    domain_ = RectDomain(ll, ur);
}


//===========================================================================
SplineSurface* SurfaceOfRevolution::geometrySurface() const
//===========================================================================
{
    double umin = domain_.umin();
    double umax = domain_.umax();
    double vmin = domain_.vmin();
    double vmax = domain_.vmax();
    double angle = umax - umin;

    // First clone the swept curve and rotate to umin
    shared_ptr<SplineCurve> curve(curve_->clone());
    rotateSplineCurve(axis_dir_, umin, *curve);

    // Sweep out surface and set correct parameter domain
    SplineSurface* ssof 
	= SweepSurfaceCreator::rotationalSweptSurface(*curve, angle,
						      location_, axis_dir_);
    ssof->basis_u().rescale(umin, umax);
    ssof->basis_v().rescale(vmin, vmax);

    return ssof;
}


//===========================================================================


} // namespace Go
