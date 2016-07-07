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
using std::streamsize;
using std::swap;


namespace Go
{


// Constructor
//===========================================================================
SurfaceOfRevolution::SurfaceOfRevolution(Point location, Point axis_dir,
                                         shared_ptr<SplineCurve> curve,
                                         bool isSwapped)
    : location_(location), axis_dir_(axis_dir), curve_(curve), isSwapped_(false)
//===========================================================================
{
    if (location_.dimension() != 3) {
        THROW("Dimension must be 3.");
        return;
    }

    axis_dir_.normalize();
    setParameterBounds(0.0, curve_->startparam(), 2.0 * M_PI, curve_->endparam());

    if (isSwapped)
        swapParameterDirection();
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

    // "Reset" swapping
    isSwapped_ = false;

    double from_upar, to_upar;
    is >> from_upar >> to_upar;

    // Need to take care of rounding errors: If upars are "roughly"
    // (0, 2*M_PI) it is probably meant *exactly* (0, 2*M_PI).
    const double pareps = 1.0e-4; // This is admittedly arbitrary...
    if (fabs(from_upar) < pareps && fabs(to_upar - 2.0*M_PI) < pareps) {
        from_upar = 0.0;
        to_upar = 2.0 * M_PI;
    }

    setParameterBounds(from_upar, curve_->startparam(), 
        to_upar, curve_->endparam());

    // Swapped flag
    int isSwapped; // 0 or 1
    is >> isSwapped;
    if (isSwapped == 0) {
        // Do nothing
    }
    else if (isSwapped == 1) {
        swapParameterDirection();
    }
    else {
        THROW("Swapped flag must be 0 or 1");
    }
}


//===========================================================================
void SurfaceOfRevolution::write(std::ostream& os) const
//===========================================================================
{
    streamsize prev = os.precision(15);
    os << dimension() << endl
       << location_ << endl
       << axis_dir_ << endl;
    curve_->write(os);

    // Bounds in the v-direction is contained in the curve
    os << domain_.umin() << " " << domain_.umax() << endl;

    if (!isSwapped()) {
        os << "0" << endl;
    }
    else {
        os << "1" << endl;
    }

    os.precision(prev);   // Reset precision to it's previous value
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
    SurfaceOfRevolution* sor 
        = new SurfaceOfRevolution(location_, axis_dir_, curve_, isSwapped_);
    sor->domain_ = domain_;
    return sor;
}


//===========================================================================
const RectDomain& SurfaceOfRevolution::parameterDomain() const
//===========================================================================
{
    if (!isSwapped())
        return domain_;

    // If parameters are swapped, we must make a swapped domain
    Array<double, 2> ll, ur;
    ll[0] = domain_.vmin();
    ll[1] = domain_.umin();
    ur[0] = domain_.vmax();
    ur[1] = domain_.umax();
    orientedDomain_ = RectDomain(ll, ur);
    return orientedDomain_;
}


//===========================================================================
RectDomain SurfaceOfRevolution::containingDomain() const
//===========================================================================
{
    return parameterDomain();
}


//===========================================================================
bool SurfaceOfRevolution::inDomain(double u, double v, double eps) const
//===========================================================================
{
    getOrientedParameters(u, v);
    Array<double, 2> pt(u, v);
    return domain_.isInDomain(pt, eps);
}


//===========================================================================
int SurfaceOfRevolution::inDomain2(double u, double v, double eps) const
//===========================================================================
{
    getOrientedParameters(u, v);
    Array<double, 2> pt(u, v);
    return domain_.isInDomain2(pt, eps);
}


//===========================================================================
bool SurfaceOfRevolution::onBoundary(double u, double v, double eps) const 
//===========================================================================
{
    getOrientedParameters(u, v);
    Array<double, 2> pnt(u, v);
    return parameterDomain().isOnBoundary(pnt, eps);
}

//===========================================================================
Point SurfaceOfRevolution::closestInDomain(double u, double v) const
//===========================================================================
{
    getOrientedParameters(u, v);
    Array<double, 2> pt(u, v);
    Array<double, 2> clo_pt;
    // Using an arbitrary tolerance... @jbt
    double tol = 1.0e-12;
    domain_.closestInDomain(pt, clo_pt, tol);
    if (!isSwapped())
        return Point(clo_pt[0], clo_pt[1]);
    else
        return Point(clo_pt[1], clo_pt[0]);
}


//===========================================================================
CurveLoop
SurfaceOfRevolution::outerBoundaryLoop(double degenerate_epsilon) const
//===========================================================================
{
    vector<shared_ptr<ParamCurve> > loop_cvs(4);

    loop_cvs[0] = getCircle(curve_->startparam());

    loop_cvs[1] = shared_ptr<ParamCurve>(curve_->clone());

    loop_cvs[2] = getCircle(curve_->endparam());
    loop_cvs[2]->reverseParameterDirection();

    loop_cvs[3] = shared_ptr<ParamCurve>(curve_->clone());
    loop_cvs[3]->reverseParameterDirection();

    CurveLoop loop(loop_cvs, degenerate_epsilon);

    return loop;
}


//===========================================================================
std::vector<CurveLoop> 
SurfaceOfRevolution::allBoundaryLoops(double degenerate_epsilon) const
//===========================================================================
{
    vector<CurveLoop> loops(1);
    loops[0] = outerBoundaryLoop(degenerate_epsilon);

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
    getOrientedParameters(upar, vpar); // In case of swapped
    Point cvpt;
    curve_->point(cvpt, vpar);
    Point lc = cvpt - location_;
    double cosu = cos(upar);
    double sinu = sin(upar);

    pt = location_
        + lc * cosu + (lc * axis_dir_) * axis_dir_ * (1.0 - cosu)
        + axis_dir_.cross(lc) * sinu;
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

    // Swap parameters, if needed
    getOrientedParameters(upar, vpar);
    int ind1 = 1;
    int ind2 = 2;
    if (isSwapped())
        swap(ind1, ind2);

    // First derivatives
    double cosu = cos(upar);
    double sinu = sin(upar);
    vector<Point> cvpts(2);
    curve_->point(cvpts, vpar, 1);
    Point lc = cvpts[0] - location_;
    Point& dl = cvpts[1];
    pts[ind1] = -lc * sinu + (lc * axis_dir_) * axis_dir_ * sinu
        + axis_dir_.cross(lc) * cosu;
    pts[ind2] = dl * cosu + (dl * axis_dir_) * axis_dir_ * (1.0 - cosu)
        + axis_dir_.cross(dl) * sinu;
    if (derivs == 1)
        return;

    // Second order and higher derivatives.
    MESSAGE("Second order or higher derivatives not yet implemented.");

}


//===========================================================================
void SurfaceOfRevolution::normal(Point& n, double upar, double vpar) const
//===========================================================================
{
    vector<Point> pts(3);
    point(pts, upar, vpar, 1);
    n = pts[1].cross(pts[2]);
    if (n.length() != 0.0)
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
    // We use the closest point function from SplineSurface. This will
    // in general lead to incorrect results!

    // Get the SplineSurface
    SplineSurface* sf = geometrySurface();

    // We first use the given seed to find values
    RectDomain curr_domain_of_interest = parameterDomain();
    if (domain_of_interest != NULL) {
        curr_domain_of_interest.intersectWith(*domain_of_interest);
    }
    double umin = curr_domain_of_interest.umin();
    double umax = curr_domain_of_interest.umax();
    double vmin = curr_domain_of_interest.vmin();
    double vmax = curr_domain_of_interest.vmax();
    double curr_seed[2];
    if (seed == NULL) {
        curr_seed[0] = 0.5 * (umin + umax);
        curr_seed[1] = 0.5 * (vmin + vmax);
    }
    else {
        curr_seed[0] = seed[0];
        curr_seed[1] = seed[1];
    }
    sf->closestPoint(pt, clo_u, clo_v, clo_pt, clo_dist, epsilon,
                     &curr_domain_of_interest, curr_seed);

    // Try to reseed
    double seeds[4][2];
    seeds[0][0] = umin;
    seeds[0][1] = vmin;
    seeds[1][0] = umax;
    seeds[1][1] = vmin;
    seeds[2][0] = umin;
    seeds[2][1] = vmax;
    seeds[3][0] = umax;
    seeds[3][1] = vmax;
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

    // Fix incorrect parametrization in u-direction by extracting
    // a circle at the current v value and call closestPoint on
    // the circle.
    getOrientedParameters(clo_u, clo_v);
    getOrientedParameters(umin, vmin);
    getOrientedParameters(umax, vmax);
    shared_ptr<Circle> circle = getCircle(clo_v);
    circle->closestPoint(pt, umin, umax, clo_u, clo_pt, clo_dist);
    getOrientedParameters(clo_u, clo_v);

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
    // This is a bit like cheating...

    SplineSurface* sf = geometrySurface();
    sf->closestBoundaryPoint(pt, clo_u, clo_v, clo_pt, clo_dist, epsilon,
                             rd, seed);

    // Fix incorrect parametrization in u-direction by extracting
    // a circle at the current v value and call closestPoint on
    // the circle.
    getOrientedParameters(clo_u, clo_v);
    shared_ptr<Circle> circle = getCircle(clo_v);
    circle->closestPoint(pt, domain_.umin(), domain_.umax(),
                         clo_u, clo_pt, clo_dist);
    getOrientedParameters(clo_u, clo_v);

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
    isSwapped_ = !isSwapped_;
}


//===========================================================================
bool SurfaceOfRevolution::isSwapped() const
//===========================================================================
{
    return isSwapped_;
}


//===========================================================================
void SurfaceOfRevolution::getOrientedParameters(double& u, double& v) const 
//===========================================================================
{
    if (isSwapped_) {
        swap(u, v);
    }
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


////===========================================================================
//bool SurfaceOfRevolution::isDegenerate(bool& b, bool& r,
//                                       bool& t, bool& l,
//                                       double tolerance) const
////===========================================================================
//{
//    MESSAGE("isDegenerate() not implemented.");
//    return false;
//}


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
void
SurfaceOfRevolution::setParameterBounds(double from_upar, double from_vpar,
                                        double to_upar, double to_vpar)
//===========================================================================
{
    if (from_upar >= to_upar )
        THROW("First u-parameter must be strictly less than second.");
    if (from_vpar >= to_vpar )
	THROW("First v-parameter must be strictly less than second.");

    getOrientedParameters(from_upar, from_vpar);
    getOrientedParameters(to_upar, to_vpar);

    // NOTE: If parameters are swapped, from_upar and from_vpar are swapped.
    // Ditto for to_upar/to_vpar.
    if (from_upar < -2.0 * M_PI || to_upar > 2.0 * M_PI)
        THROW("u-parameters must be in [-2pi, 2pi].");
    if (to_upar - from_upar > 2.0 * M_PI)
        THROW("(to_upar - from_upar) must not exceed 2pi.");

    Array<double, 2> ll(from_upar, from_vpar);
    Array<double, 2> ur(to_upar, to_vpar);
    domain_ = RectDomain(ll, ur);
}


//===========================================================================
SplineSurface* SurfaceOfRevolution::geometrySurface() const
//===========================================================================
{
    return createSplineSurface();
}


//===========================================================================
SplineSurface* SurfaceOfRevolution::createSplineSurface() const
//===========================================================================
{
    double umin = domain_.umin();
    double umax = domain_.umax();
    double vmin = domain_.vmin();
    double vmax = domain_.vmax();
    double angle = umax - umin;

    // First clone the swept curve and rotate to umin
    shared_ptr<SplineCurve> curve(curve_->clone());
    GeometryTools::rotateSplineCurve(axis_dir_, umin, *curve);

    // Sweep out surface and set correct parameter domain
    SplineSurface* ssof 
        = SweepSurfaceCreator::rotationalSweptSurface(*curve, angle,
                                                      location_, axis_dir_);
    ssof->basis_u().rescale(umin, umax);
    ssof->basis_v().rescale(vmin, vmax);

    if (isSwapped())
        ssof->swapParameterDirection();

    return ssof;
}


//===========================================================================
shared_ptr<Circle> SurfaceOfRevolution::getCircle(double vpar) const
//===========================================================================
{
    Point lambda;
    curve_->point(lambda, vpar);
    Point lc = lambda - location_;
    double lcv = lc * axis_dir_;
    Point lcvv = lcv * axis_dir_;

    Point centre = location_ + lcvv;
    double radius = sqrt(lc * lc - lcv * lcv);
    Point x_axis = lc - lcvv; // Not normalized
    if (x_axis.length() == 0.0) {
        x_axis = Point(1.0, 0.0, 0.0);
        if (x_axis.cross(axis_dir_).length() == 0.0) {
            x_axis = Point(0.0, 1.0, 0.0);
        }
    }
    shared_ptr<Circle> circle(new Circle(radius, centre, axis_dir_, x_axis));
    double umin = domain_.umin();
    double umax = domain_.umax();
    circle->setParamBounds(umin, umax);
    return circle;
}


//===========================================================================


} // namespace Go
