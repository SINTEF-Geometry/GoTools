//===========================================================================
//                                                                           
// File: IntersectorFuncConst.C                                                
//                                                                           
// Created: Wed Sep 22 09:59:10 2004                                         
//                                                                           
// Author: Sverre Briseid <Sverre.Briseid@sintef.no>
//                                                                           
// Revision: $Id: IntersectorFuncConst.C,v 1.23 2006-11-03 14:43:22 jbt Exp $
//                                                                           
// Description:
//                                                                           
//===========================================================================


#include "GoTools/intersections/IntersectorFuncConst.h"
#include "GoTools/utils/CompositeBox.h"
#include "GoTools/utils/DirectionCone.h"
#include "GoTools/intersections/ParamPointInt.h"
#include "GoTools/intersections/Param1FunctionInt.h"
#include "GoTools/intersections/Param2FunctionInt.h"
#include "GoTools/intersections/IntersectionUtils.h"
#include "GoTools/intersections/IntersectionPool.h"
#include "GoTools/geometry/SplineUtils.h"
#include "GoTools/intersections/Param0FunctionInt.h"
#include "GoTools/intersections/PlaneInt.h"
#include "GoTools/geometry/PointCloud.h"
#include "GoTools/geometry/GeometryTools.h"
#include "GoTools/creators/CreatorsUtils.h"
#include "GoTools/geometry/SplineDebugUtils.h"
#include "GoTools/creators/SurfaceCreators.h"
#include "GoTools/creators/CurveCreators.h"

using namespace Go;


//===========================================================================
IntersectorFuncConst::IntersectorFuncConst(shared_ptr<ParamFunctionInt> func,
					   shared_ptr<ParamFunctionInt> C,
					   shared_ptr<GeoTol> epsge, Intersector *prev,
					   int eliminated_parameter,
					   double eliminated_value)
//===========================================================================
    : Intersector(epsge, prev), func_int_(func), C_(C)
{
    shared_ptr<IntersectionPool> parent_pool;
    if (prev) {
	parent_pool = prev->getIntPool();
    }
    int_results_ =
	shared_ptr<IntersectionPool>(new IntersectionPool(func,
							  C,
							  parent_pool, 
							  eliminated_parameter,
							  eliminated_value));
}

///////////////////////////////////////////////////////////////////
//
//  Purpose    : Destructor
//
//  Written by : Sverre Briseid, 0404
//
/////////////////////////////////////////////////////////////////

IntersectorFuncConst::~IntersectorFuncConst()
{
  // Currently empty
}

///////////////////////////////////////////////////////////////////
//
//  Purpose    : Do interseption test between two parametric objects
//               Currently only box testing is performed
//
//  Written by : Sverre Briseid, SINTEF, Sep 2004
//
/////////////////////////////////////////////////////////////////

int IntersectorFuncConst::performInterception()
{
  int do_intercept;

  // Currently only box testing is done.
  // @@@ vsk, Note that only expanded boxes are used. This is not good enough.

  // Get the boxes belonging to the objects
  CompositeBox box = func_int_->compositeBox();

  Point c_pt;
  C_->point(c_pt, NULL);
  // Check overlap
  if ((box.low()[0] - epsge_->getEpsge() < c_pt[0]) &&
      (box.high()[0] + epsge_->getEpsge() > c_pt[0]))
    {
      // Check if the boxes both lie in an area of size equal to the
      // tolerance

      // Some check for this should be included in the bounding box
//       if (box.zeroSize(box2, epsge_)
      if ((box.low()[0] > c_pt[0] - epsge_->getEpsge()) &&
        (box.high()[0] < c_pt[0] + epsge_->getEpsge()))
	{ // The box is epsge-degenerate.
	  do_intercept = 2; // microcase.
	}
      else
      // or something like this
      {
	do_intercept = 1;
      }
    }
  else 
    do_intercept = 0;

  return do_intercept;
}

///////////////////////////////////////////////////////////////////
//
//  Purpose    : Do simple case test between parametric function
//               and constant. We want to know if more than one
//               solution is possible (checking monotonicity).
//
//  Written by : Sverre Briseid, SINTEF, Sep 2004
//
/////////////////////////////////////////////////////////////////

int IntersectorFuncConst::simpleCase()
{

    // If both objects are pts it is a simple case.
    if (func_int_->getParam0FunctionInt() == 0) {
	Point mon_dir;
	bool monotone = func_int_->monotone(mon_dir); // strictly monotone, i.e.

	if (monotone) { // At most one solution possible.
	    return 1;
	}

	Point c_pt;
	C_->point(c_pt, NULL);
	// We then check the boundingbox.
	// If constant C_ lies inside range we may have multiple solutions.
	CompositeBox box = func_int_->compositeBox(); // 1d box.
	if ((box.low()[0] - epsge_->getEpsge() < c_pt[0]) &&
	    (box.high()[0] + epsge_->getEpsge() > c_pt[0])) {
	    return 0; // Possibly more than one solution (at least one).
	}
    }

    // Otherwise it is a simple case
    return 1;
}

///////////////////////////////////////////////////////////////////
//
//  Purpose    : Check linearity. Not relly applicable in 1D
//
//  Written by : Vibeke Skytt, 12.04
//
/////////////////////////////////////////////////////////////////
bool IntersectorFuncConst::isLinear()
{
    return false;
}

///////////////////////////////////////////////////////////////////
//
//  Purpose    : Action in case of linearity
//
//  Written by : Vibeke Skytt, 12.04
//
/////////////////////////////////////////////////////////////////
int IntersectorFuncConst::linearCase()
{
    return 0;
}


void IntersectorFuncConst::print_objs()
{
    std::ofstream debug("data/debug.g2");
    Param2FunctionInt *par2_func1 = func_int_->getParam2FunctionInt();
    shared_ptr<SplineSurface> sf_3d;
    if (par2_func1)
	{
	    shared_ptr<ParamSurface> srf = par2_func1->getParamSurface();
	    shared_ptr<SplineSurface> spline_sf =
		dynamic_pointer_cast<SplineSurface, ParamSurface>(srf);
	    sf_3d = SurfaceCreators::insertParamDomain(*spline_sf);
	    sf_3d->writeStandardHeader(debug);
	    sf_3d->write(debug);
	}
    shared_ptr<Param0FunctionInt> par0_func2 =
	dynamic_pointer_cast<Param0FunctionInt, ParamFunctionInt>(C_);
    if (par0_func2) {
	double c = par0_func2->getValue(); // The constant z-value for the plane.
	shared_ptr<PlaneInt> plane_int(new PlaneInt(0.0, 0.0, 1.0, -c));
	BoundingBox bd_box_sf = sf_3d->boundingBox();
	Point mid_pt = 0.5*(bd_box_sf.low() + bd_box_sf.high());
	double length_x = 2.0*(bd_box_sf.low().dist(bd_box_sf.high()));
	double length_y = length_x;
	shared_ptr<SplineSurface> plane_sf_repr(plane_int->surface(mid_pt, length_x, length_y));
	plane_sf_repr->writeStandardHeader(debug);
	plane_sf_repr->write(debug);
    }
}


///////////////////////////////////////////////////////////////////
//
//  Purpose    : 
//
//  Written by : Sverre Briseid, SINTEF, Sep 2004
//
/////////////////////////////////////////////////////////////////

int IntersectorFuncConst::getBoundaryIntersections()
{

  // Check if a previous intersector at the same level exist.
  // In that case the boundary intersections are computed already,
  // and it is nothing to do.
  if (prev_intersector_ != 0 &&
      prev_intersector_->numParams() == numParams())
    return 0;  // Not necessary to compute

  // The boundary intersections must be computed.
  // Make the appropriate intersectors
  // Get the boundary objects
  // Keep the boundary intersectors until this function goes out of 
  // scope to be able to search through already defined intersector for 
  // relevant intersection results.  
  shared_ptr<Intersector> bd_intersector;
  std::vector<shared_ptr<BoundaryFunctionInt> > bd_objs;
  func_int_->getBoundaryObjects(bd_objs);

  int ki;
  for (ki=0; ki < int(bd_objs.size()); ki++)
    {
      bd_intersector = lowerOrderIntersector(bd_objs[ki]->getObject(),
					     C_, this,
					     bd_objs[ki]->getDir(),
					     bd_objs[ki]->getPar());
      bd_intersector->compute();
      int_results_->includeReducedInts(bd_intersector->getIntPool());
    }

  return 1;  // Boundary intersections computed
}


//===========================================================================
void IntersectorFuncConst::printDebugInfo()
//===========================================================================
{
    int npar1 = func_int_->numParams();
    int npar2 = C_->numParams();
    double rel_par_res = epsge_->getRelParRes();
    int ki, kj;
    if (npar1 > 0 &&
	(prev_intersector_ == 0
	 || prev_intersector_->numParams() > numParams())) {

	// Write objects to file
	std::ofstream debug("geom_out.g2");
	Param1FunctionInt *par1_func1 = func_int_->getParam1FunctionInt();
	Param2FunctionInt *par2_func1 = func_int_->getParam2FunctionInt();
	if (par1_func1) {
	    SplineCurve *cv = par1_func1->getParamCurve()->geometryCurve();
	    // As the cv lives in 1d we add the parameter domain and
	    // let it live in z=0.0.
	    shared_ptr<SplineCurve> cv_2d
		= CurveCreators::insertParamDomain(*cv, rel_par_res);
	    writeSpaceParamCurve(*cv_2d, debug, 0.0);
	}
	if (par2_func1) {
	    shared_ptr<ParamSurface> srf = par2_func1->getParamSurface();
	    shared_ptr<SplineSurface> spline_sf
		= dynamic_pointer_cast<SplineSurface, ParamSurface>(srf);
	    // As the sf lives in 1d we add the parameter domain.
	    shared_ptr<SplineSurface> sf_3d
		= SurfaceCreators::insertParamDomain(*spline_sf);
	    sf_3d->writeStandardHeader(debug);
	    sf_3d->write(debug);
	    // We then write the plane.
	    shared_ptr<Param0FunctionInt> p0_func_int
		= dynamic_pointer_cast<Param0FunctionInt, ParamFunctionInt>(C_);
	    double c = p0_func_int->getValue(); // Our z value.
	    shared_ptr<PlaneInt> plane_int(new PlaneInt(0.0, 0.0, 1.0, -c));
	    BoundingBox bd_box = sf_3d->boundingBox();
	    Point mid_pt = 0.5*(bd_box.low() + bd_box.high());
	    double length_x = 2.0*(bd_box.low().dist(bd_box.high()));
	    double length_y = length_x;
	    shared_ptr<SplineSurface>
		plane_sf(plane_int->surface(mid_pt, length_x, length_y));
	    plane_sf->writeStandardHeader(debug);
	    plane_sf->write(debug);
	}
    }

    if (npar1 == 1 && 
	(prev_intersector_ == 0
	 || prev_intersector_->numParams() > numParams())) {
 	std::ofstream debug("geom_crv_out.g2");
	Param1FunctionInt *par1_func1 = func_int_->getParam1FunctionInt();
	if (par1_func1) {
	    SplineCurve *cv = par1_func1->getParamCurve()->geometryCurve();
	    // @@sbr Object may be a bd cv for a sf; in that case we
	    // must insert constant parameter value in missing
	    // parameter direction, and add function values.
	    shared_ptr<SplineCurve> cv_2d
		= CurveCreators::insertParamDomain(*cv);
	    writeSpaceParamCurve(*cv_2d, debug, 0.0);
	}
    }

    if (npar1 == 2) {
	bool write_sfs = false;
	int nmb_rec = nmbRecursions();
	if (nmb_rec == 1) {
	    write_sfs = true;
	}
	if (write_sfs) {
	    std::ofstream debug("geom_sfs_out.g2");
	    Param2FunctionInt *par2_func1 = func_int_->getParam2FunctionInt();
	    shared_ptr<SplineSurface> sf_3d;
	    if (par2_func1) {
		shared_ptr<ParamSurface> srf = par2_func1->getParamSurface();
		shared_ptr<SplineSurface> spline_sf
		    = dynamic_pointer_cast<SplineSurface, ParamSurface>(srf);
		sf_3d = SurfaceCreators::insertParamDomain(*spline_sf);
		sf_3d->writeStandardHeader(debug);
		sf_3d->write(debug);
	    }
	    shared_ptr<Param0FunctionInt> par0_func2
		= dynamic_pointer_cast<Param0FunctionInt, ParamFunctionInt>(C_);
	    if (par0_func2) {
		double c = par0_func2->getValue(); // The constant
						   // z-value for the
						   // plane.
		shared_ptr<PlaneInt> plane_int(new PlaneInt(0.0, 0.0,
							    1.0, -c));
		BoundingBox bd_box_sf = sf_3d->boundingBox();
		Point mid_pt = 0.5*(bd_box_sf.low() + bd_box_sf.high());
		double length_x = 2.0*(bd_box_sf.low().dist(bd_box_sf.high()));
		double length_y = length_x;
		shared_ptr<SplineSurface>
		    plane_sf_repr(plane_int->surface(mid_pt,
						     length_x, length_y));
		plane_sf_repr->writeStandardHeader(debug);
		plane_sf_repr->write(debug);
	    }
	}

	double ta1[2], tb1[2];
	std::cout << "================================================"
		  << std::endl;
	std::cout << "Domain 1: ";
	for (ki=0; ki<npar1; ki++) {
	    ta1[ki] = func_int_->startParam(ki);
	    tb1[ki] = func_int_->endParam(ki);
	    std::cout << ta1[ki] << " ";
	    std::cout << tb1[ki] << " ";
	}
	std::cout << std::endl;
	std::cout << "Domain 2: ";
	std::cout << std::endl;

	std::vector<shared_ptr<IntersectionPoint> > ipoint;
	int_results_->getIntersectionPoints(ipoint);
	std::cout << "Intersection points : " << ipoint.size() << std::endl;
	// We also write to file the 3d-visualization of the
	// intersection point(s).
	std::vector<double> pts;
	for (ki=0; ki < int(ipoint.size()); ki++) {
	    for (kj=0; kj<npar1+npar2; kj++) {
		std::cout << ipoint[ki]->getPar(kj) << "  ";
	    }
	    std::cout << std::endl;
	    Point space_pt = ipoint[ki]->getPoint();
	    if (space_pt.size() == 3) {
		pts.insert(pts.end(), space_pt.begin(), space_pt.end());
	    } else if ((space_pt.size() == 1)
		       && (ipoint[ki]->numParams1() == 2)) {
		const double* par1 = ipoint[ki]->getPar1();
		pts.push_back(par1[0]);
		pts.push_back(par1[1]);
		pts.push_back(space_pt[0]);
	    }
	}
	
	bool write_pts = (pts.size() > 0);
	if (write_pts) {
	    std::ofstream debug("int_pts.g2");
	    int dim = 3;
	    int nmb_pts = (int)pts.size()/dim;
	    Go::PointCloud3D pt_cloud(pts.begin(), nmb_pts);
	    pt_cloud.writeStandardHeader(debug);
	    pt_cloud.write(debug);
	}
    }

    if (npar1 > 0) {
	if (getenv("DEBUG_PAR") && (*getenv("DEBUG_PAR"))=='1') {
	    int_results_->writeDebug();
	}
    }
}


//===========================================================================
