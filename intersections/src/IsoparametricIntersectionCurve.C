//===========================================================================
//                                                                           
// File: IsoparametricIntersectionCurve.C                                    
//                                                                           
// Created: Wed Apr 20 16:15:33 2005                                         
//                                                                           
// Author: Odd A. Andersen <Odd.Andersen@sintef.no>
//                                                                           
// Revision:
// $Id: IsoparametricIntersectionCurve.C,v 1.9 2006-11-03 14:15:12 jbt Exp $
//                                                                           
// Description:
//                                                                           
//===========================================================================

#include "GoTools/intersections/IntersectionCurve.h"
#include "GoTools/intersections/ParamSurfaceInt.h"
#include "GoTools/intersections/Param2FunctionInt.h"
#include "GoTools/geometry/SplineCurve.h"
#include <stdexcept>


using std::vector;
using std::swap;
using std::runtime_error;
using std::shared_ptr;


namespace Go {


//===========================================================================
shared_ptr<ParamCurve> IsoparametricIntersectionCurve::getCurve() const
//===========================================================================
{
    ASSERT(isopar_geom_curve_.get() != 0); // this curve should be precalculated!
    return isopar_geom_curve_;
}

//=========================================================================== 
shared_ptr<ParamCurve> IsoparametricIntersectionCurve::getParamCurve(int obj_nmb) const
//===========================================================================
{
    shared_ptr<ParamCurve> result;
    switch (obj_nmb) {
    case 1:
	result = isopar_param_curve_1_;
	break;
    case 2:
	result = isopar_param_curve_2_;
	break;
    default:
	throw std::logic_error("Argument to getParamCurve() should be 1 or 2.");
    }
    if (result.get() == 0) {
	MESSAGE("Warning;  Returned isocurve is a zero pointer.  It should have been "
		"precalculated, but this functionality is not yet implemented. "
		"Complete precalculate_iso_curves() to get rid of this problem.");
    }
    return result;
}


//===========================================================================
void IsoparametricIntersectionCurve::precalculate_iso_curves(int isopar)
//===========================================================================
{
    const shared_ptr<IntersectionPoint>& pt_1 = ipoints_.front();
    const shared_ptr<IntersectionPoint>& pt_2 = ipoints_.back();

    const bool use_obj_2 = pt_1->numParams1() <= isopar;

    const ParamObjectInt* const ref_obj = use_obj_2 ? pt_1->getObj2() : pt_1->getObj1();
    ASSERT(ref_obj->numParams() == 2); // only makes sense on 2-manifolds.

    double isoval = pt_1->getPar(isopar);
    bool running_u = isopar != 0 && isopar != pt_1->numParams1();
    int running_ix = running_u ? isopar - 1 : isopar + 1;
    
    double endpar_1 = pt_1->getPar(running_ix);
    double endpar_2 = pt_2->getPar(running_ix);
    if (endpar_1 > endpar_2) {
	swap(endpar_1, endpar_2);
    }
    if (endpar_2 - endpar_1 < pt_1->parameterTolerance(running_ix)) {
	// paramspan too short to define a curve
	throw Zero_Parameter_Span_Error();
    }

    const ParamSurfaceInt* const ref_obj_surf = dynamic_cast<const ParamSurfaceInt*>(ref_obj);
    if (!ref_obj_surf) {
	const Param2FunctionInt* const ref_obj_func = dynamic_cast<const Param2FunctionInt*>(ref_obj);
	if (!ref_obj_func) {
	    throw runtime_error("IntersectionCurve::precalculate_iso_curves() called "
				"on unexpected object.");
	} else {
	    isopar_geom_curve_ = ref_obj_func->getIsoCurve(endpar_1, endpar_2, isoval, running_u);
	}
    } else {
	isopar_geom_curve_ = ref_obj_surf->getIsoCurve(endpar_1, endpar_2, isoval, running_u);
    }

    shared_ptr<ParamCurve>& linear_param_curve = 
	use_obj_2 ? isopar_param_curve_2_ : isopar_param_curve_1_;
    shared_ptr<ParamCurve>& other_param_curve = 
	use_obj_2 ? isopar_param_curve_1_ : isopar_param_curve_2_;

    // constructing isoparametric param curve
    vector<double> coefs;
    if (!running_u) {
	coefs.push_back(isoval);
    }
    coefs.push_back(endpar_1);
    coefs.push_back(isoval);
    coefs.push_back(endpar_2);
    coefs.push_back(isoval);

    vector<double> param(4);
    param[0] = param[1] = endpar_1;
    param[2] = param[3] = endpar_2;
    linear_param_curve = 
	shared_ptr<ParamCurve>(new SplineCurve(2, 2, &param[0], &coefs[0], 2));
    
    MESSAGE("WARNING: Computing of parametric curve of an isocurve in the non-isoparametric "
	    "object is not yet implemented in 'precalculate_iso_curves()'");
    other_param_curve = shared_ptr<ParamCurve>(); // @@ IMPLEMENT THIS
}

}; // end namespace Go
