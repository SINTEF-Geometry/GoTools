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


#include "GoTools/geometry/SurfaceOfLinearExtrusion.h"
#include "GoTools/geometry/SplineCurve.h"
#include "GoTools/geometry/Line.h"

#include <limits>
#include <fstream>


using std::vector;
using std::cout;
using std::endl;
using std::numeric_limits;
using std::streamsize;
using std::swap;


namespace Go
{


// Constructor
//===========================================================================
SurfaceOfLinearExtrusion::SurfaceOfLinearExtrusion(shared_ptr<SplineCurve> curve,
                                                   Point axis_dir,                                                   
                                                   bool isSwapped)
    : curve_(curve), axis_dir_(axis_dir), isSwapped_(false)
//===========================================================================
{
    if (axis_dir_.dimension() != 3) {
        THROW("Dimension must be 3.");
        return;
    }

    double inf = numeric_limits<double>::infinity();
    setParameterBounds(curve_->startparam(), -inf, curve_->endparam(), inf);

    if (isSwapped) {
        swapParameterDirection();
    }
}


// Destructor
//===========================================================================
SurfaceOfLinearExtrusion::~SurfaceOfLinearExtrusion()
//===========================================================================
{
}


//===========================================================================
void SurfaceOfLinearExtrusion::read (std::istream& is)
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

    curve_ = shared_ptr<SplineCurve>(new SplineCurve);
    curve_->read(is);

    axis_dir_.resize(dim);
    is >> axis_dir_;

    // "Reset" swapping
    isSwapped_ = false;

    int isBounded; 
    is >> isBounded;
    if (isBounded == 0) {
        // Unbounded in v direction

        // NB: See comment on parameter sequence above!
        double from_upar, to_upar;
        is >> from_upar >> to_upar;

        double inf = numeric_limits<double>::infinity();
        setParameterBounds(from_upar, -inf, to_upar, inf);
    }
    else if (isBounded == 1) {
        // NB: See comment on parameter sequence above!
        double from_upar, from_vpar, to_upar, to_vpar;
        is >> from_upar >> to_upar
            >> from_vpar >> to_vpar;

	try
	{
	    setParameterBounds(from_upar, from_vpar, to_upar, to_vpar);
	}
	catch (...)
	{ // We want the read routine to continue reading data, not a good strategy to throw before all object data is parsed.
	    MESSAGE("Failed setting parameter bounds.");
	}
    }
    else {
        THROW("Bounded flag must be 0 or 1");
    }

    // double inf = numeric_limits<double>::infinity();
    // setParameterBounds(curve_->startparam(), -inf,
    //                    curve_->endparam(), inf);

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
void SurfaceOfLinearExtrusion::write(std::ostream& os) const
//===========================================================================
{
    streamsize prev = os.precision(15);
    os << dimension() << endl;
    curve_->write(os); // Assuming spline curve.
    os << axis_dir_ << endl;

    // NB: Mind the parameter sequence!
    if (!isBounded()) {
        os << "0" << endl
           << domain_.umin() << " " << domain_.umax() << endl;
    }
    else {
        os << "1" << endl
           << domain_.umin() << " " << domain_.umax() << endl
           << domain_.vmin() << " " << domain_.vmax() << endl;
    }

    if (!isSwapped()) {
        os << "0" << endl;
    }
    else {
        os << "1" << endl;
    }

    os.precision(prev);   // Reset precision to it's previous value
}


//===========================================================================
int SurfaceOfLinearExtrusion::dimension() const
//===========================================================================
{
    return axis_dir_.dimension();
}


//===========================================================================
ClassType SurfaceOfLinearExtrusion::instanceType() const
//===========================================================================
{
    return classType();
}


//===========================================================================
ClassType SurfaceOfLinearExtrusion::classType()
//===========================================================================
{
    return Class_SurfaceOfLinearExtrusion;
}


//===========================================================================
BoundingBox SurfaceOfLinearExtrusion::boundingBox() const
//===========================================================================
{
    // A rather inefficient hack...
    SurfaceOfLinearExtrusion* sor = const_cast<SurfaceOfLinearExtrusion*>(this);
    SplineSurface* tmp = sor->geometrySurface();
    if (tmp == NULL) {
        MESSAGE("Method not supported yet!");
        BoundingBox dummy_box;
        return dummy_box;
    }
    BoundingBox box = tmp->boundingBox();
    delete tmp;
    return box;
}

//===========================================================================
SurfaceOfLinearExtrusion* SurfaceOfLinearExtrusion::clone() const
//===========================================================================
{
    SurfaceOfLinearExtrusion* sor 
        = new SurfaceOfLinearExtrusion(curve_, axis_dir_, isSwapped_);
    sor->domain_ = domain_;
    return sor;
}


//===========================================================================
const RectDomain& SurfaceOfLinearExtrusion::parameterDomain() const
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
RectDomain SurfaceOfLinearExtrusion::containingDomain() const
//===========================================================================
{
    return parameterDomain();
}


//===========================================================================
bool SurfaceOfLinearExtrusion::inDomain(double u, double v, double eps) const
//===========================================================================
{
    getOrientedParameters(u, v);
    Array<double, 2> pt(u, v);
    return domain_.isInDomain(pt, eps);
}


//===========================================================================
int SurfaceOfLinearExtrusion::inDomain2(double u, double v, double eps) const
//===========================================================================
{
    getOrientedParameters(u, v);
    Array<double, 2> pt(u, v);
    return domain_.isInDomain2(pt, eps);
}


//===========================================================================
bool SurfaceOfLinearExtrusion::onBoundary(double u, double v, double eps) const 
//===========================================================================
{
    getOrientedParameters(u, v);
    Array<double, 2> pnt(u, v);
    return parameterDomain().isOnBoundary(pnt, eps);
}

//===========================================================================
Point SurfaceOfLinearExtrusion::closestInDomain(double u, double v) const
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
SurfaceOfLinearExtrusion::outerBoundaryLoop(double degenerate_epsilon) const
//===========================================================================
{
#if 0
    MESSAGE("outerBoundaryLoop() not implemented, returning an empty vector");
    CurveLoop dummy_loop;

    return dummy_loop;
#else
    double umin = domain_.umin();
    double umax = domain_.umax();
    double vmin = domain_.vmin();
    double vmax = domain_.vmax();
    if (!isBounded()) {
        if (isSwapped_) {
            limitDomain(umin, umax);
        } else {
            limitDomain(vmin, vmax);
        }
    }

    // if (isSwapped_) {
    //     std::cout << "Parameter directions are swapped, expect failure ..." << std::endl;
    // }

    vector<shared_ptr<ParamCurve> > loop_cvs(4);

    shared_ptr<SplineCurve> edge0(curve_->clone());
    const int dim = edge0->dimension();
    // We must then translate the coefs with vmin*axis_dir_.
    for (auto iter = edge0->coefs_begin(); iter != edge0->coefs_end(); iter += dim) {
        for (int kj = 0; kj < dim; ++kj) {
            iter[kj] += vmin*axis_dir_[kj];
        }
    }
    loop_cvs[0] = edge0;

#if 0
    std::cout << "domain: " << domain_.umin() << ", " << domain_.umax() << ", " <<
        domain_.vmin() << ", " << domain_.vmax() << std::endl;
#endif

    Point cv_pt_min, cv_pt_max;
    // getOrientedParameters(umin, vmin);
    // getOrientedParameters(umax, vmax);
    curve_->point(cv_pt_min, umin);
    curve_->point(cv_pt_max, umax);
    Point linear_pt_min = vmin*axis_dir_;
    Point linear_pt_max = vmax*axis_dir_;

    Point lr = cv_pt_max + linear_pt_min;
    Point ur = cv_pt_max + linear_pt_max;
    loop_cvs[1] = shared_ptr<Line>(new Line(lr, ur, vmin, vmax));

    shared_ptr<SplineCurve> edge2(curve_->clone());
    edge2->reverseParameterDirection();
    // We must then translate the coefs with vmin*axis_dir_.
    for (auto iter = edge2->coefs_begin(); iter != edge2->coefs_end(); iter += dim) {
        for (int kj = 0; kj < dim; ++kj) {
            iter[kj] += vmax*axis_dir_[kj];
        }
    }
    loop_cvs[2] = edge2;

    Point ul = cv_pt_min + linear_pt_max;
    Point ll = cv_pt_min + linear_pt_min;
    loop_cvs[3] = shared_ptr<Line>(new Line(ul, ll, vmin, vmax));

    if (isSwapped_) {
        reverse(loop_cvs.begin(), loop_cvs.end());
        for (auto iter = loop_cvs.begin(); iter != loop_cvs.end(); ++iter) {
            (*iter)->reverseParameterDirection();
        }
    }

    CurveLoop loop(loop_cvs, degenerate_epsilon);

    return loop;
#endif
}


//===========================================================================
std::vector<CurveLoop> 
SurfaceOfLinearExtrusion::allBoundaryLoops(double degenerate_epsilon) const
//===========================================================================
{
    vector<CurveLoop> loops(1);
    loops[0] = outerBoundaryLoop(degenerate_epsilon);

    return loops;
}


//===========================================================================
DirectionCone SurfaceOfLinearExtrusion::normalCone() const
//===========================================================================
{
    // A rather inefficient hack...
    SurfaceOfLinearExtrusion* sof = const_cast<SurfaceOfLinearExtrusion*>(this);
    SplineSurface* tmp = sof->geometrySurface();
    if (tmp == NULL) {
        MESSAGE("Method not supported yet!");
        DirectionCone empty_dc;
        return empty_dc;
    }
    DirectionCone dc = tmp->normalCone();
    delete tmp;
    return dc;
}


//===========================================================================
DirectionCone SurfaceOfLinearExtrusion::tangentCone(bool pardir_is_u) const
//===========================================================================
{
    // A rather inefficient hack...
    SurfaceOfLinearExtrusion* sof = const_cast<SurfaceOfLinearExtrusion*>(this);
    SplineSurface* tmp = sof->geometrySurface();
    if (tmp == NULL) {
        MESSAGE("Method not supported yet!");
        DirectionCone empty_dc;
        return empty_dc;
    }
    DirectionCone dc = tmp->tangentCone(pardir_is_u);
    delete tmp;
    return dc;
}


//===========================================================================
void SurfaceOfLinearExtrusion::point(Point& pt, double upar, double vpar) const
//===========================================================================
{
    getOrientedParameters(upar, vpar); // In case of swapped
    Point cvpt;
    curve_->point(cvpt, upar);

    pt = cvpt + vpar*axis_dir_;
}


//===========================================================================
void SurfaceOfLinearExtrusion::point(std::vector<Point>& pts, 
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
    vector<Point> cvpts(2);
    curve_->point(cvpts, upar, 1);
    pts[ind1] = cvpts[1];
    pts[ind2] = axis_dir_;
    if (derivs == 1)
        return;

    // Second order and higher derivatives.
    MESSAGE("Second order or higher derivatives not yet implemented.");

}


//===========================================================================
void SurfaceOfLinearExtrusion::normal(Point& n, double upar, double vpar) const
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
SurfaceOfLinearExtrusion::constParamCurves(double parameter, bool pardir_is_u) const
//===========================================================================
{
    MESSAGE("constParamCurves() not yet implemented");
    vector<shared_ptr<ParamCurve> > res;
    return res;
}



//===========================================================================
bool SurfaceOfLinearExtrusion::isBounded() const
//===========================================================================
{
    // It is enough to check the v-direction, since the u-direction is
    // always bounded.

    return domain_.vmin() > -numeric_limits<double>::infinity() &&
        domain_.vmax() < numeric_limits<double>::infinity();
}


//===========================================================================
SurfaceOfLinearExtrusion*
SurfaceOfLinearExtrusion::subSurface(double from_upar, double from_vpar,
                                double to_upar, double to_vpar,
                                double fuzzy) const
//===========================================================================
{
    SurfaceOfLinearExtrusion* sor = clone();
    sor->setParameterBounds(from_upar, from_vpar, to_upar, to_vpar);
    return sor;
}


//===========================================================================
vector<shared_ptr<ParamSurface> >
SurfaceOfLinearExtrusion::subSurfaces(double from_upar, double from_vpar,
                                 double to_upar, double to_vpar,
                                 double fuzzy) const
//===========================================================================
{
    vector<shared_ptr<ParamSurface> > res;
    shared_ptr<SurfaceOfLinearExtrusion> sor(subSurface(from_upar, from_vpar,
                                                   to_upar, to_vpar));
    res.push_back(sor);
    return res;
}


//===========================================================================
double 
SurfaceOfLinearExtrusion::nextSegmentVal(int dir, double par, bool forward,
                                    double tol) const
//===========================================================================
{
    MESSAGE("Does not make sense. Return arbitrarily zero.");
    return 0.0;
}


//===========================================================================
void SurfaceOfLinearExtrusion::closestPoint(const Point& pt,
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
    // in general lead to incorrect results! @@sbr201601 Why?

    // We first use the given seed to find values
    RectDomain curr_domain_of_interest = parameterDomain();
    if (domain_of_interest != NULL) {
        curr_domain_of_interest.intersectWith(*domain_of_interest);
    }

    double umin = curr_domain_of_interest.umin();
    double umax = curr_domain_of_interest.umax();
    double vmin = curr_domain_of_interest.vmin();
    double vmax = curr_domain_of_interest.vmax();
    if (!isBounded()) {
        if (isSwapped_) {
            limitDomain(umin, umax);
        } else {
            limitDomain(vmin, vmax);
        }

        Array<double, 2> ll, ur;
        ll[0] = umin;
        ll[1] = vmin;
        ur[0] = umax;
        ur[1] = vmax;
        curr_domain_of_interest = RectDomain(ll, ur);
    }
    double curr_seed[2];
    if (seed == NULL) {
        curr_seed[0] = 0.5 * (umin + umax);
        curr_seed[1] = 0.5 * (vmin + vmax);
    }
    else {
        curr_seed[0] = seed[0];
        curr_seed[1] = seed[1];
    }

#if 0
    MESSAGE("Currently using the generic closest point algorithm from ParamSurface.");
    ParamSurface::closestPoint(pt, clo_u, clo_v, clo_pt, clo_dist, epsilon,
                               &curr_domain_of_interest, curr_seed);
    return;
#else
    // Get the SplineSurface. If isSwapped_ is true then the spline surafce is flipped and all is good.
    SplineSurface* sf = geometrySurface();

    if (sf == NULL) {
        MESSAGE("Method not available yet!");
        return;
    }

    double u_span = umax - umin;
    double v_span = vmax - vmin;
    // if (std::max(u_span, v_span)* epsilon > 1.0) {
    //     std::cout << "Domain too large, closest point will most likely fail! Consider using smaller epsilon. Prod: " <<
    //         std::max(u_span, v_span)* epsilon << std::endl;
    //     std::cout << "epsilon: " << epsilon << std::endl;
    //     std::cout << "DEFAULT_PARAMETER_EPSILON: " << DEFAULT_PARAMETER_EPSILON << std::endl;
    // }
    // @@sbr201601 We need to convert from epsgeo to epspar.
    sf->closestPoint(pt, clo_u, clo_v, clo_pt, clo_dist, epsilon,
                     &curr_domain_of_interest, curr_seed);

    // // Try to reseed
    // double seeds[4][2];
    // seeds[0][0] = umin;
    // seeds[0][1] = vmin;
    // seeds[1][0] = umax;
    // seeds[1][1] = vmin;
    // seeds[2][0] = umin;
    // seeds[2][1] = vmax;
    // seeds[3][0] = umax;
    // seeds[3][1] = vmax;
    // Point tmppt(3);
    // double tmpu, tmpv, tmpdist;
    // for (int i = 0; i < 4; ++i) {
    //     sf->closestPoint(pt, tmpu, tmpv, tmppt, tmpdist, epsilon,
    //                      &curr_domain_of_interest, seeds[i]);
    //     if (tmpdist < clo_dist - epsilon) {
    //         clo_u = tmpu;
    //         clo_v = tmpv;
    //         clo_pt = tmppt;
    //         clo_dist = tmpdist;
    //     }
    // }

    delete sf;
#endif
}


//===========================================================================
void SurfaceOfLinearExtrusion::closestBoundaryPoint(const Point& pt,
                                                    double&        clo_u,
                                                    double&        clo_v, 
                                                    Point&       clo_pt,
                                                    double&        clo_dist,
                                                    double epsilon,
                                                    const RectDomain* rd,
                                                    double *seed) const
//===========================================================================
{
#if 0
    // The (proper) loop based method needs handling of bd cv orientation & handling if isSwapped_ == true.
    CurveLoop cv_loop = outerBoundaryLoop(epsilon);
    double global_clo_dist = numeric_limits<double>::infinity();
    for (size_t ki = 0; ki < cv_loop.size(); ++ki) {
        double* cv_seed = NULL;
        if (seed != NULL) {
            cv_seed = &seed[ki%2];
        }
        Point cv_clo_pt;
        double cv_clo_par, cv_clo_dist;
        cv_loop[ki]->closestPoint(pt, cv_clo_par, cv_clo_pt, cv_clo_dist);
        if (cv_clo_dist < global_clo_dist) {
            global_clo_dist = cv_clo_dist;
            clo_pt = cv_clo_pt;
            if (ki == 0) {
                clo_u = cv_clo_par;
                clo_v = domain_.vmin();
            } else if (ki == 1) {
                clo_u = domain_.umax();
                clo_v = cv_clo_par;
            } else if (ki == 2) {
                clo_u = domain_.umin() + domain_.umax() - cv_clo_par;
                clo_v = domain_.vmax();
            } else {
                clo_u = domain_.umin();
                clo_v = domain_.vmin() + domain_.vmax() - cv_clo_par;
            }
        }
    }

    clo_dist = global_clo_dist;

#else
    // A rather inefficient hack...
    SplineSurface* sf = geometrySurface();
    if (sf == NULL) {
        MESSAGE("Method not implemented yet!");
        return;
    }
    sf->closestBoundaryPoint(pt, clo_u, clo_v, clo_pt, clo_dist, epsilon,
                             rd, seed);

    delete sf;
#endif

}


//===========================================================================
void
SurfaceOfLinearExtrusion::getBoundaryInfo(Point& pt1, Point& pt2,
                                     double epsilon, SplineCurve*& cv,
                                     SplineCurve*& crosscv,
                                     double knot_tol) const
//===========================================================================
{
    MESSAGE("getBoundaryInfo() not yet implemented");
}


//===========================================================================
void SurfaceOfLinearExtrusion::turnOrientation()
//===========================================================================
{
    swapParameterDirection();
}


//===========================================================================
void SurfaceOfLinearExtrusion::swapParameterDirection()
//===========================================================================
{
//    std::cout << "SurfaceOfLinearExtrusion: Swapping parameter directions!" << std::endl;
    isSwapped_ = !isSwapped_;
}


//===========================================================================
bool SurfaceOfLinearExtrusion::isSwapped() const
//===========================================================================
{
    return isSwapped_;
}


//===========================================================================
void SurfaceOfLinearExtrusion::getOrientedParameters(double& u, double& v) const 
//===========================================================================
{
    if (isSwapped_) {
        swap(u, v);
    }
}


//===========================================================================
double SurfaceOfLinearExtrusion::area(double tol) const
//===========================================================================
{
    MESSAGE("area() not implemented.");
    return -1.0;
}


//===========================================================================
void SurfaceOfLinearExtrusion::reverseParameterDirection(bool direction_is_u)
//===========================================================================
{
    MESSAGE("reverseParameterDirection() not implemented.");
}


//===========================================================================
void SurfaceOfLinearExtrusion::getDegenerateCorners(vector<Point>& deg_corners,
                                               double tol) const
//===========================================================================
{
    MESSAGE("getDegenerateCorners() not implemented.");
}


//===========================================================================
  void SurfaceOfLinearExtrusion::getCornerPoints(std::vector<std::pair<Point,Point> >& corners) const
//===========================================================================
{
    MESSAGE("getCornerPoints() not implemented.");
}

//===========================================================================
void
SurfaceOfLinearExtrusion::setParameterBounds(double from_upar, double from_vpar,
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
    if (from_upar < curve_->startparam() || to_upar > curve_->endparam())
        THROW("u-parameters must be in domain of curve.");

    Array<double, 2> ll(from_upar, from_vpar);
    Array<double, 2> ur(to_upar, to_vpar);
    domain_ = RectDomain(ll, ur);
}


//===========================================================================
SplineSurface* SurfaceOfLinearExtrusion::geometrySurface() const
//===========================================================================
{
    return createSplineSurface();
}


//===========================================================================
SplineSurface* SurfaceOfLinearExtrusion::createSplineSurface() const
//===========================================================================
{
    // First handle the case if not bounded
    // @@sbr201601 Why? Is it not better to require that the domain is bounded
    // to allow this function to be called? This is misleading ...
    double umin = domain_.umin();
    double umax = domain_.umax();
    double vmin = domain_.vmin();
    double vmax = domain_.vmax();
#if 0
    std::cout << "umin: " << umin << ", umax: " << umax << ", vmin: " << vmin << ", vmax: " << vmax << std::endl;
#endif
    if (!isBounded()) {
        if (isSwapped_) {
            limitDomain(umin, umax);
        } else {
            limitDomain(vmin, vmax);
        }
    }

    Point transl_vmin = vmin*axis_dir_;
    Point transl_vmax = vmax*axis_dir_;

    vector<double> coefs_vmin(curve_->coefs_begin(), curve_->coefs_end());
    vector<double> coefs_vmax(curve_->coefs_begin(), curve_->coefs_end());
    const int num_cv_coefs = curve_->numCoefs();
    const int dim = curve_->dimension();
    for (int ki = 0; ki < num_cv_coefs; ++ki) {
        for (int kj = 0; kj < dim; ++kj) {
            coefs_vmin[ki*dim+kj] += transl_vmin[kj];
            coefs_vmax[ki*dim+kj] += transl_vmax[kj];
        }
    }

    vector<double> sf_coefs = coefs_vmin;
    sf_coefs.insert(sf_coefs.end(), coefs_vmax.begin(), coefs_vmax.end());

    vector<double> knots_vdir(4, vmin);
    knots_vdir[2] = vmax;
    knots_vdir[3] = vmax;
    SplineSurface* spline_sf = new SplineSurface(num_cv_coefs, 2,
                                                 curve_->order(), 2,
                                                 curve_->basis().begin(), knots_vdir.begin(),
                                                 sf_coefs.begin(),
                                                 dim);


    if (isSwapped_) {
//        std::cout << "SurfaceOfLinearExtrusion: Flipping the spline surface!" << std::endl;
        spline_sf->swapParameterDirection();
    }

#if 0
    MESSAGE("Debugging!");
    std::ofstream debug("tmp/sf_debug.g2");
    spline_sf->writeStandardHeader(debug);
    spline_sf->write(debug);
#endif

    return spline_sf;
}


//===========================================================================
void SurfaceOfLinearExtrusion::limitDomain(double& vmin, double& vmax) const
//===========================================================================
{
    double max = 1.0e8;//6;//8; // "Large" number...
    // if (isSwapped_) {
    //     if (umin == -numeric_limits<double>::infinity())
    //         umin = -max;
    //     if (umax == numeric_limits<double>::infinity())
    //         umax = max;

    // } else {
    // Assuming that a params have been swapped if surface is swapped.
    if (vmin == -numeric_limits<double>::infinity())
        vmin = -max;
    if (vmax == numeric_limits<double>::infinity())
        vmax = max;
//    }
}


} // namespace Go
