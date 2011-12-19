//===========================================================================
//                                                                           
// File: IntersectorAlgPar.C                                                 
//                                                                           
// Created: Wed Jan 26 15:44:50 2005                                         
//                                                                           
// Author: Sverre Briseid <Sverre.Briseid@sintef.no>
//                                                                           
// Revision: $Id: IntersectorAlgPar.C,v 1.15 2006-11-03 14:43:22 jbt Exp $
//                                                                           
// Description:
//                                                                           
//===========================================================================


#include "GoTools/intersections/IntersectorAlgPar.h"
#include "GoTools/intersections/Line2DInt.h"
#include "GoTools/intersections/PlaneInt.h"
#include "GoTools/intersections/ParamCurveInt.h"
#include "GoTools/intersections/ParamSurfaceInt.h"
#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/intersections/Spline1FunctionInt.h"
#include "GoTools/intersections/Spline2FunctionInt.h"
#include "GoTools/intersections/Param0FunctionInt.h"
#include "GoTools/intersections/Par1FuncIntersector.h"
#include "GoTools/intersections/Par2FuncIntersector.h"
#include "GoTools/geometry/GeometryTools.h"
#include "GoTools/intersections/IntersectionUtils.h"
#include "GoTools/creators/CreatorsUtils.h"
#include "GoTools/creators/SurfaceCreators.h"

#include <vector>


using namespace Go;
using std::vector;


//===========================================================================
IntersectorAlgPar::IntersectorAlgPar(shared_ptr<AlgObjectInt> alg_obj,
				     shared_ptr<ParamObjectInt> param_obj,
				     shared_ptr<GeoTol> epsge, Intersector* prev,
				     int eliminated_parameter,
				     double eliminated_value)
    : Intersector(epsge, prev), param_int_(param_obj), algobj_int_(alg_obj)
//===========================================================================
{
    // @@sbr So far we expect the input param_obj to be either a
    // ParamCurveInt or a ParamSurfaceInt.
    if (dynamic_pointer_cast<ParamCurveInt, ParamObjectInt>(param_int_)) {
	shared_ptr<ParamCurveInt> cv_int =
	    dynamic_pointer_cast<ParamCurveInt, ParamObjectInt>(param_int_);
	shared_ptr<AlgObj2DInt> alg_obj2d_int =
	    dynamic_pointer_cast<AlgObj2DInt, AlgObjectInt>(alg_obj);
	ASSERT(alg_obj2d_int.get() != 0);

	// We further expect the curve to be a 2d spline curve.
	shared_ptr<SplineCurve> spline_cv_2d =
	    dynamic_pointer_cast<SplineCurve, ParamCurve>(cv_int->getParamCurve());
	ASSERT(spline_cv_2d.get() != 0 && spline_cv_2d->dimension() == 2);

	// We then plot the parametric curve into the implicit curve to get a functional.
       int_func_const_ = shared_ptr<IntersectorFuncConst>
	   (insertCurveInAlgobj(spline_cv_2d.get(), alg_obj2d_int.get(), epsge, prev,
				eliminated_parameter, eliminated_value));
    } else if (dynamic_pointer_cast<ParamSurfaceInt, ParamObjectInt>(param_int_)) {
	shared_ptr<ParamSurfaceInt> sf_int =
	    dynamic_pointer_cast<ParamSurfaceInt, ParamObjectInt>(param_int_);
	shared_ptr<AlgObj3DInt> alg_obj3d_int =
	    dynamic_pointer_cast<AlgObj3DInt, AlgObjectInt>(alg_obj);
	ASSERT(alg_obj3d_int.get() != 0);

	// We further expect the curve to be a 3d spline surface.
	shared_ptr<SplineSurface> spline_sf_3d =
	    dynamic_pointer_cast<SplineSurface, ParamSurface>(sf_int->getParamSurface());
	ASSERT(spline_sf_3d.get() != 0 && spline_sf_3d->dimension() == 3);

	// We then plot the parametric surface into the implicit surface to get a functional.
       int_func_const_ = shared_ptr<IntersectorFuncConst>
	   (insertSurfaceInAlgobj(spline_sf_3d.get(), alg_obj3d_int.get(), epsge, prev,
				  eliminated_parameter, eliminated_value));
    } else {
	THROW("Unexpected parametric object!");
    }

}


//===========================================================================
void IntersectorAlgPar::compute(bool compute_at_boundary)
//===========================================================================
{
    return int_func_const_->compute(compute_at_boundary);
}


//===========================================================================
void IntersectorAlgPar::getResult(vector<shared_ptr<IntersectionPoint> >& int_points,
				  vector<shared_ptr<IntersectionCurve> >& int_curves)
//===========================================================================
{
    vector<shared_ptr<IntersectionPoint> > func_int_points;
    vector<shared_ptr<IntersectionCurve> > func_int_curves;
    int_func_const_->getResult(func_int_points, func_int_curves);
    // @@sbr Possibly transform points to format in objects?
    // void transformIntersectionResult()?
    int_points = func_int_points;
    int_curves = func_int_curves;
}


//===========================================================================
int IntersectorAlgPar::numParams() const
//===========================================================================
{
    return int_func_const_->numParams();
}


//===========================================================================
void IntersectorAlgPar::print_objs()
//===========================================================================
{
    return int_func_const_->print_objs();
}


//===========================================================================
int IntersectorAlgPar::getBoundaryIntersections()
//===========================================================================
{
    return int_func_const_->getBoundaryIntersections();
}


//===========================================================================
int IntersectorAlgPar::performInterception()
//===========================================================================
{
    return int_func_const_->performInterception();
}


//===========================================================================
int IntersectorAlgPar::simpleCase()
//===========================================================================
{
    return int_func_const_->simpleCase();
}


//===========================================================================
bool IntersectorAlgPar::isLinear()
//===========================================================================
{
    return int_func_const_->isLinear();
}


//===========================================================================
bool IntersectorAlgPar::complexityReduced()
//===========================================================================
{
    return int_func_const_->complexityReduced();
}


//===========================================================================
void IntersectorAlgPar::handleComplexity()
//===========================================================================
{
    return int_func_const_->handleComplexity();
}


//===========================================================================
int IntersectorAlgPar::checkCoincidence()
//===========================================================================
{
    return int_func_const_->checkCoincidence();
}
    

//===========================================================================
void IntersectorAlgPar::microCase()
//===========================================================================
{
    return int_func_const_->microCase();
}


//===========================================================================    
int IntersectorAlgPar::updateIntersections()
//===========================================================================
{
    return int_func_const_->updateIntersections();
}


//===========================================================================
int IntersectorAlgPar::linearCase()
//===========================================================================
{
    return int_func_const_->linearCase();
}


//===========================================================================
int IntersectorAlgPar::doSubdivide()
//===========================================================================
{
    return int_func_const_->doSubdivide();
}


//===========================================================================
void IntersectorAlgPar::printDebugInfo()
//===========================================================================
{
    return int_func_const_->printDebugInfo();
}


//===========================================================================
shared_ptr<IntersectorFuncConst>
IntersectorAlgPar::insertCurveInAlgobj(SplineCurve* cv,
				       AlgObj2DInt* alg_obj2d_int,
				       shared_ptr<GeoTol> epsge,
				       Intersector* intersector,
				       int eliminated_parameter,
				       double eliminated_value)
//===========================================================================
{
    shared_ptr<SplineCurve> spline_cv(IntersectionUtils::
				      insertCvInAlgcv(*cv, alg_obj2d_int));

    shared_ptr<Spline1FunctionInt>
	spline1_func_int(new Spline1FunctionInt(spline_cv));
    // We intersect the spline functional with 0.
    shared_ptr<Param0FunctionInt> C(new Param0FunctionInt(0.0));
    shared_ptr<IntersectorFuncConst>
	int_func_const(new Par1FuncIntersector(spline1_func_int, C, 
					       epsge, intersector,
					       eliminated_parameter,
					       eliminated_value));

    return int_func_const;
}


//===========================================================================
shared_ptr<IntersectorFuncConst>
IntersectorAlgPar::insertSurfaceInAlgobj(SplineSurface* sf,
					 AlgObj3DInt* alg_obj3d_int,
					 shared_ptr<GeoTol> epsge,
					 Intersector* intersector,
					 int eliminated_parameter,
					 double eliminated_value)
//===========================================================================
{
    // s(u,v) = (s_1(u,v), s_2(u,v), s_3(u,v))
    shared_ptr<SplineSurface> spline_sf(IntersectionUtils::
					insertSfInAlgsf(*sf, alg_obj3d_int));

#ifdef INTERSECTIONS_DEBUG
    std::ofstream debug("data/debug.g2");
    shared_ptr<SplineSurface> sf_3d(SurfaceCreators::
				    insertParamDomain(*spline_sf));
    sf_3d->writeStandardHeader(debug);
    sf_3d->write(debug);

    double c = 0.0;
    shared_ptr<PlaneInt> plane_int(new PlaneInt(0.0, 0.0, 1.0, -c));
    BoundingBox bd_box_sf = sf_3d->boundingBox();
    Point mid_pt = 0.5*(bd_box_sf.low() + bd_box_sf.high());
    double length_x = 2.0*(bd_box_sf.low().dist(bd_box_sf.high()));
    double length_y = length_x;
    shared_ptr<SplineSurface> plane_sf_repr(plane_int->surface(mid_pt,
							       length_x,
							       length_y));
    plane_sf_repr->writeStandardHeader(debug);
    plane_sf_repr->write(debug);
#endif // INTERSECTIONS_DEBUG

    shared_ptr<Spline2FunctionInt>
	spline2_func_int(new Spline2FunctionInt(spline_sf));
    // We intersect the spline functional with 0.
    shared_ptr<Param0FunctionInt> C(new Param0FunctionInt(0.0));
    shared_ptr<IntersectorFuncConst> 
	int_func_const(new Par2FuncIntersector(spline2_func_int, C,
					       epsge, intersector,
					       eliminated_parameter,
					       eliminated_value));

    return int_func_const;
}
