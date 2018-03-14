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


#include "GoTools/creators/OffsetSurface.h"
#include "GoTools/geometry/Factory.h"
#include "GoTools/geometry/ClassType.h"
#include "GoTools/creators/CreatorsOffsetUtils.h"
#include "GoTools/creators/CurveCreators.h"


using std::vector;
using std::endl;


namespace Go
{


// Constructor
//===========================================================================
OffsetSurface::OffsetSurface(shared_ptr<ParamSurface> param_sf,
                             double offset_dist, double epsgeo, bool self_int)
    : surface_(param_sf), offset_dist_(offset_dist), epsgeo_(epsgeo), self_int_(self_int)
//===========================================================================
{
    try
    {
        createOffsetOuterBdLoop();
    }
    catch (...)
    {
        THROW("Failed creating boundary looops.");
    }
}


// Destructor
//===========================================================================
OffsetSurface::~OffsetSurface()
//===========================================================================
{
}


//===========================================================================
void OffsetSurface::read(std::istream& is)
//===========================================================================
{
    is >> offset_dist_;
    int self_int_val;
    is >> self_int_val;
    if ((self_int_val != 0) && (self_int_val != 1))
    {
        MESSAGE("Unexpected value!");
    }
    self_int_ = (self_int_val == 1) ? true : false;

    int instance_type;
    is >> instance_type;
    ClassType type = ClassType(instance_type); // Needs this conversion

    shared_ptr<GeomObject> goobject(Factory::createObject(type));
    shared_ptr<ParamSurface> tmp_srf 
	= dynamic_pointer_cast<ParamSurface, GeomObject>(goobject);
    ALWAYS_ERROR_IF(tmp_srf.get() == 0,
		    "Can not read this instance type");

    try
    {
	tmp_srf->read(is);
	surface_ = tmp_srf;
    }
    catch (...)
    { // We want the read routine to continue reading data, not a good strategy to throw before all object data is parsed.
	MESSAGE("Failed reading the surface.");
    }

}

//===========================================================================
void OffsetSurface::write(std::ostream& os) const
//===========================================================================
{
    std::streamsize prev = os.precision(15);

    os << offset_dist_ << endl;

    if (!self_int_)
	os << "0";
    else
	os << "1";
    os << endl;

    os << surface_->instanceType() << std::endl;
    surface_->write(os);
    os << endl;

    os.precision(prev);   // Reset precision to it's previous value
}


//===========================================================================
BoundingBox OffsetSurface::boundingBox() const
//===========================================================================
{
    // We create an over-estimate by using the base surface boundingbox and enlarging it with
    // fabs(offset_dist) in all both directions along all axis.
    BoundingBox sf_box = surface_->boundingBox();
    Point offset_pt(surface_->dimension());
    offset_pt.setValue(fabs(offset_dist_));

    Point box_low = sf_box.low() - offset_pt;
    Point box_high = sf_box.high() + offset_pt;
    BoundingBox bd_box(box_low, box_high);

    return bd_box;
}


//===========================================================================
int OffsetSurface::dimension() const
//===========================================================================
{
    return surface_->dimension();
}


//===========================================================================
ClassType OffsetSurface::instanceType() const
//===========================================================================
{
    return classType();
}


//===========================================================================
SplineSurface* OffsetSurface::asSplineSurface()
//===========================================================================
{
    MESSAGE("asSplineSurface() not implemented.");

    return nullptr;  // Default behaviour
}


//===========================================================================
SplineSurface* OffsetSurface::getSplineSurface() 
//===========================================================================
{
    MESSAGE("getSplineSurface() not implemented.");

    return nullptr;  // Default behaviour
}


//===========================================================================
const Domain& OffsetSurface::parameterDomain() const
//===========================================================================
{
    return surface_->parameterDomain();
}


//===========================================================================
RectDomain OffsetSurface::containingDomain() const
//===========================================================================
{
    return surface_->containingDomain();
}


//===========================================================================
bool OffsetSurface::isBounded() const
//===========================================================================
{
    return surface_->isBounded();
}


//===========================================================================
bool OffsetSurface::inDomain(double u, double v, double eps) const
//===========================================================================
{
    return surface_->inDomain(u, v, eps);
}


//===========================================================================
int OffsetSurface::inDomain2(double u, double v, double eps) const
//===========================================================================
{
    return surface_->inDomain2(u, v, eps);
}


//===========================================================================
bool OffsetSurface::onBoundary(double u, double v, double eps) const
//===========================================================================
{
    return surface_->onBoundary(u, v, eps);
}


//===========================================================================
Point OffsetSurface::closestInDomain(double u, double v) const
//===========================================================================
{
    return surface_->closestInDomain(u, v);

}


//===========================================================================
void OffsetSurface::setParameterDomain(double u1, double u2, double v1, double v2)
//===========================================================================
{
    surface_->setParameterDomain(u1, u2, v1, v2);
}


//===========================================================================
CurveLoop OffsetSurface::outerBoundaryLoop(double degenerate_epsilon) const
//===========================================================================
{
    return offset_outer_bd_loop_;
}


//===========================================================================
std::vector<CurveLoop> OffsetSurface::allBoundaryLoops(double degenerate_epsilon) const
//===========================================================================
{
    vector<CurveLoop> loops;
    MESSAGE("allBoundaryLoops() not implemented");

    return loops;
}


//===========================================================================    
DirectionCone OffsetSurface::normalCone() const
//===========================================================================    
{
    DirectionCone dir_cone;
    MESSAGE("normalCone() not implemented");

    return dir_cone;
}
    

//===========================================================================
DirectionCone OffsetSurface::tangentCone(bool pardir_is_u) const
//===========================================================================    
{
    DirectionCone dir_cone;
    MESSAGE("tangentCone() not implemented");

    return dir_cone;
}

   
//===========================================================================
void OffsetSurface::point(Point& pt, double upar, double vpar) const
//===========================================================================
{
    Point sf_pt = surface_->ParamSurface::point(upar, vpar);
    Point normal;
    surface_->normal(normal, upar, vpar);
    pt = sf_pt + normal*offset_dist_;
}


//===========================================================================
void OffsetSurface::point(std::vector<Point>& pts, 
                          double upar, double vpar,
                          int derivs,
                          bool u_from_right,
                          bool v_from_right,
                          double resolution) const
//===========================================================================
{
    DEBUG_ERROR_IF(derivs < 0, "Negative number of derivatives makes no sense.");
    int totpts = (derivs + 1)*(derivs + 2)/2;
    DEBUG_ERROR_IF((int)pts.size() < totpts, "The vector of points must have sufficient size.");

    // Calling blend_s1421.
    vector<Point> offset_pt(((derivs+1)*(derivs+2)/2) + 1); // Derivs & normal in the exact surface.
    vector<Point> base_pt(((derivs+1)*(derivs+2)/2) + 1); // Derivs & normal.

    if (surface_->instanceType() == Class_SplineSurface)
    {
        shared_ptr<SplineSurface> spline_sf = dynamic_pointer_cast<SplineSurface>(surface_);
        Point epar(upar, vpar);
        int ind_u=0;             /* Pointer into knot vector                       */
        int ind_v=0;             /* Pointer into knot vector                       */
        int kstat = 0;
        OffsetUtils::blend_s1421(spline_sf.get(), offset_dist_, derivs, epar, ind_u, ind_v,
                                 offset_pt, base_pt, &kstat);
        if (kstat != 0)
        {
            MESSAGE("WARNING: The returned status value was not 0 as expected!");
        }
        pts = offset_pt;
    }
    else
    {
        THROW("Only supported for SplineSurface's!");
    }
}


//===========================================================================
void OffsetSurface::normal(Point& n, double upar, double vpar) const
//===========================================================================
{
    vector<Point> pts(4);
    point(pts, upar, vpar, 1);

    n = pts[1]%pts[2];
    n.normalize();
}


//===========================================================================
void OffsetSurface::evalGrid(int num_u, int num_v, 
                             double umin, double umax, 
                             double vmin, double vmax,
                             std::vector<double>& points,
                             double nodata_val) const
//===========================================================================
{
    MESSAGE("evalGrid() not implemented");
}


//===========================================================================
Point OffsetSurface::getInternalPoint(double& u, double& v) const
//===========================================================================
{
    Point sf_pt = getInternalPoint(u, v);
    Point offset_pt;
    point(offset_pt, u, v);

    return offset_pt;
}


//===========================================================================
std::vector<shared_ptr<ParamCurve> >
OffsetSurface::constParamCurves(double parameter, bool pardir_is_u) const
//===========================================================================
{
    vector<shared_ptr<ParamCurve> > iso_cvs;
    MESSAGE("constParamCurves() not implemented");

    return iso_cvs;
}


//===========================================================================
std::vector<shared_ptr<ParamSurface> >
OffsetSurface::subSurfaces(double from_upar, double from_vpar,
                           double to_upar, double to_vpar,
                           double fuzzy) const
//===========================================================================
{
    vector<shared_ptr<ParamSurface> > sub_sfs;
    MESSAGE("constParamCurves() not implemented");

    return sub_sfs;
}


//===========================================================================
double OffsetSurface::nextSegmentVal(int dir, double par, bool forward, double tol) const
//===========================================================================
{
    return surface_->nextSegmentVal(dir, par, forward, tol);
}


//===========================================================================
void OffsetSurface::closestBoundaryPoint(const Point& pt,
                                         double&        clo_u,
                                         double&        clo_v, 
                                         Point&       clo_pt,
                                         double&        clo_dist,
                                         double epsilon,
                                         const RectDomain* rd,
                                         double *seed) const
//===========================================================================
{
    RectDomain domain = containingDomain();
    if (!rd)
	rd = &domain;

    CurveLoop curve_loop = outerBoundaryLoop(epsilon);
    double loop_clo_dist = MAXDOUBLE;
    for (auto cv : curve_loop)
    {
        double cv_clo_t, cv_clo_dist;
        Point cv_clo_pt;
        cv->closestPoint(pt, cv_clo_t, cv_clo_pt, cv_clo_dist);
        if (cv_clo_dist < loop_clo_dist)
        {
            clo_pt = cv_clo_pt;
            clo_dist = cv_clo_dist;
        }
    }

    // We must then find the (u, v) pair in the surface corresponding to clo_pt.
    // We need to get parameter values for clo_pt (suppose we could use seed if set in above routine...).
    double sf_u, sf_v, sf_clo_dist;
    Point sf_clo_pt;
    surface_->closestPoint(clo_pt, sf_u, sf_v, sf_clo_pt, sf_clo_dist,
			   epsilon, rd, seed);
    // Now, the parameter point (tmp_u, tmp_v) should be in the parameter
    // domain of the surface. If so, we return happily.
    // Otherwise, we find the closest point in the domain.

    // VSK, 0902. First check if the point found on the boundary and the point in the surface
    // is the same. In that case, we are done
    if (clo_pt.dist(sf_clo_pt) <= epsilon)
    {
	clo_u = sf_u;
	clo_v = sf_v;
	clo_pt = sf_clo_pt;
	clo_dist = clo_pt.dist(pt); 
    }
    else
    {
	//	clo_dist = tmp_cld;
//	const CurveBoundedDomain& dom = parameterDomain();
	const Domain& dom = parameterDomain();
	bool is_in_domain = false;
	double domain_tol = std::max(epsilon, 1.0e-7);
	try
        {
            // Test is rather unstable when point is on/near boundary.
            is_in_domain = dom.isInDomain(Vector2D(sf_u, sf_v), domain_tol);
	} catch (...)
        {
            // 	MESSAGE("Failed deciding whether point was in domain.");
            is_in_domain = true;
	}

	if (is_in_domain)
        {
            clo_u = sf_u;
            clo_v = sf_v;
            clo_pt = sf_clo_pt;
            clo_dist = clo_pt.dist(pt); //	clo_dist = tmp_cld;
	}
        else
        {
            Vector2D par_pt;
            // @afr: Again, I use (spatial) epsilon for domain comparisons...
            dom.closestInDomain(Vector2D(sf_u, sf_v), par_pt, epsilon);
            point(clo_pt, par_pt[0], par_pt[1]);
            clo_u = par_pt[0];
            clo_v = par_pt[1];
            clo_dist = clo_pt.dist(pt);	
	}
    }
}


//===========================================================================
void OffsetSurface::getBoundaryInfo(Point& pt1, Point& pt2,
                                    double epsilon, SplineCurve*& cv,
                                    SplineCurve*& crosscv, double knot_tol) const
//===========================================================================
{
    MESSAGE("getBoundaryInfo() not implemented");
}


//===========================================================================
void OffsetSurface::turnOrientation()
//===========================================================================
{
    surface_->swapParameterDirection();
}


//===========================================================================
void OffsetSurface::reverseParameterDirection(bool direction_is_u)
//===========================================================================
{
    surface_->reverseParameterDirection(direction_is_u);
    offset_dist_ *= -1.0;
}


//===========================================================================
void OffsetSurface::swapParameterDirection()
//===========================================================================
{
    surface_->swapParameterDirection();
    offset_dist_ *= -1.0;
}


//===========================================================================
double OffsetSurface::area(double tol) const
//===========================================================================
{
    double area = -1.0;
    MESSAGE("area() not implemented");

    return area;
}


//===========================================================================
void OffsetSurface::getDegenerateCorners(std::vector<Point>& deg_corners, double tol) const
//===========================================================================
{
    MESSAGE("getDegenerateCorners() not implemented");
}


//===========================================================================
void OffsetSurface::getCornerPoints(std::vector<std::pair<Point,Point> >& corners) const
//===========================================================================
{
    MESSAGE("getCornerPoints() not implemented");
}


//===========================================================================
bool OffsetSurface::isPlanar(Point& normal, double tol)
//===========================================================================
{
    return surface_->isPlanar(normal, tol);
}


//===========================================================================
int OffsetSurface::ElementOnBoundary(int elem_ix, double eps)
//===========================================================================
{
    return surface_->ElementOnBoundary(elem_ix, eps);
}


//===========================================================================
int OffsetSurface::ElementBoundaryStatus(int elem_ix, double eps)
//===========================================================================
{
    return surface_->ElementBoundaryStatus(elem_ix, eps);
}


//===========================================================================
void OffsetSurface::createOffsetOuterBdLoop()
//===========================================================================
{
    MESSAGE("createOffsetOuterBdLoop() under construction!");

    if (offset_outer_bd_loop_.size() == 0)
    {
        const double deg_eps = epsgeo_;
        vector<shared_ptr<ParamCurve> > offset_loop_cvs;
        CurveLoop sf_outer_bd_loop = surface_->outerBoundaryLoop(deg_eps);
        for (auto cv_iter = sf_outer_bd_loop.begin(); cv_iter != sf_outer_bd_loop.end(); ++cv_iter)
        {
            shared_ptr<ParamCurve> par_cv;
            shared_ptr<SplineCurve> offset_cv =
                CurveCreators::offsetCurveNormalDir(par_cv,
                                                    *cv_iter, surface_,
                                                    epsgeo_, offset_dist_);
            if (offset_cv.get() == nullptr)
            {
                THROW("Offset curve was not created!");
            }
            offset_loop_cvs.push_back(offset_cv);
        }

        const bool allow_fix = false;
        offset_outer_bd_loop_ = CurveLoop(offset_loop_cvs, epsgeo_, allow_fix);        
    }
}


} // namespace Go
