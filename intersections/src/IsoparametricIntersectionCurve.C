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

#include "GoTools/intersections/IntersectionCurve.h"
#include "GoTools/intersections/ParamSurfaceInt.h"
#include "GoTools/intersections/Param2FunctionInt.h"
#include "GoTools/geometry/SplineCurve.h"
#include <stdexcept>


using std::vector;
using std::swap;
using std::runtime_error;


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
