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

#include <algorithm>
#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/geometry/GeometryTools.h"

//#ifdef __BORLANDC__
#include <iterator> // For back_inserter.  Should be required by VC++ and GCC as well...
//#endif

using std::vector;
using std::back_inserter;

namespace Go {


//===========================================================================
SplineSurface* SplineSurface::subSurface(double from_upar,
					 double from_vpar,
					 double to_upar,
					 double to_vpar,
					 double fuzzy) const
//===========================================================================
{
    if (from_upar >= to_upar) {
	THROW("First u-parameter must be strictly less than second.");
    }
    if (from_vpar >= to_vpar) {
	THROW("First v-parameter must be strictly less than second.");
    }
    if (from_upar < startparam_u()-fuzzy || from_vpar < startparam_v()-fuzzy) {
	THROW("Subsurface defined outside surface.");
    }

    // Check for periodic case.
    bool uper = to_upar > endparam_u() + fuzzy;
    bool vper = to_vpar > endparam_v() + fuzzy;
    if (uper || vper) {
	// Handle the periodic case by representing surface
	// as a curve, then calling subCurve(), which handles
	// the situation properly.
	// First in the u direction.
	shared_ptr<SplineCurve> temp_cv = GeometryTools::representSurfaceAsCurve(*this, 1);
	shared_ptr<SplineCurve> temp_sub(temp_cv->subCurve(from_upar, to_upar));
	shared_ptr<SplineSurface> surf
	    = GeometryTools::representCurveAsSurface(*temp_sub, 1, basis_v_, rational());
	BsplineBasis new_ubasis = surf->basis_u();
	// Then in the v direction.
	temp_cv = GeometryTools::representSurfaceAsCurve(*surf, 2);
	temp_sub.reset(temp_cv->subCurve(from_vpar, to_vpar));
	surf = GeometryTools::representCurveAsSurface(*temp_sub, 2, new_ubasis, rational());
	// We have to clone the return value because we cannot return a
	// shared pointer from this function.
	return surf->clone();
    }

    // If boundaries are close to existing knots, we snap.
    // Otherwise insertKnot() will not perform very well.
    basis_u().knotIntervalFuzzy(from_upar, fuzzy);
    basis_u().knotIntervalFuzzy(to_upar, fuzzy);
    basis_v().knotIntervalFuzzy(from_vpar, fuzzy);
    basis_v().knotIntervalFuzzy(to_vpar, fuzzy);

    int ord_u = order_u(); // u-order of the curve
    int ord_v = order_v(); // v-order of the curve

    vector<double> knots_u(ord_u, from_upar);
    knots_u.insert(knots_u.end(), ord_u, to_upar);
    vector<double> knots_v(ord_v, from_vpar);
    knots_v.insert(knots_v.end(), ord_v, to_vpar);

//     int i;
//     for (i = 0; i < ord_u; ++i)
// 	knots_u.push_back(from_upar);
//     for (i = 0; i < ord_u; ++i)
// 	knots_u.push_back(to_upar);
//     for (i = 0; i < ord_v; ++i)
// 	knots_v.push_back(from_vpar);
//     for (i = 0; i < ord_v; ++i)
// 	knots_v.push_back(to_vpar);

    vector<double> new_knots_u, new_knots_v;

    set_difference(knots_u.begin(), knots_u.end(),
		   basis_u().begin(), basis_u().end(),
		   back_inserter(new_knots_u));

    set_difference(knots_v.begin(), knots_v.end(),
		   basis_v().begin(), basis_v().end(),
		   back_inserter(new_knots_v));

    SplineSurface surface_copy;
    bool need_to_insert_knots
	= (new_knots_u.size() != 0) || (new_knots_v.size() != 0);
    if (need_to_insert_knots) {
	surface_copy = *this;
	surface_copy.insertKnot_u(new_knots_u);
	surface_copy.insertKnot_v(new_knots_v);	
    }
    const SplineSurface& the_surface
	= need_to_insert_knots ? surface_copy : (*this);

    typedef vector<double>::const_iterator iter;
    // Iterator to first occurence of from_upar
    iter bu = find(the_surface.basis_u().begin(),
		   the_surface.basis_u().end(), from_upar);
    // Iterator to one after last occurence of to_upar
    iter eu = find(bu + ord_u,
		   the_surface.basis_u().end(), to_upar) + ord_u;

    // Iterator to first occurence of from_vpar
    iter bv = find(the_surface.basis_v().begin(),
		   the_surface.basis_v().end(), from_vpar);
    // Iterator to one after last occurence of to_vpar
    iter ev = find(bv + ord_v,
		   the_surface.basis_v().end(), to_vpar) + ord_v;

    // Constructing the subsurface

    // Intermediate calculation: construct the coefficient vector.
    // This step is necessary for surfaces due to the way the coefficients
    // are stored in a SplineSurface.
    int coefs_dim
	= rational_ ? the_surface.dimension()+1 : the_surface.dimension();
    vector<double> newcoefs;
    int num_u = (int)(eu - bu) - ord_u;
    int num_v = (int)(ev - bv) - ord_v;
    int offset_u = (int)(bu - the_surface.basis_u().begin());
    int offset_v = (int)(bv - the_surface.basis_v().begin());

    iter tempiter;
    for (int iv = offset_v; iv < num_v+offset_v; ++iv) {
	for (int iu = offset_u; iu < num_u+offset_u; ++iu) {
	    tempiter = !rational_
		? the_surface.coefs_begin()
		: the_surface.rcoefs_begin();
	    tempiter += coefs_dim * (iv*the_surface.numCoefs_u() + iu);
	    for (int dd = 0; dd < coefs_dim; ++dd) {
		newcoefs.push_back(*(tempiter+dd));
	    }
	}
    }

    SplineSurface* the_subSurface
	= new SplineSurface(num_u, num_v,
			      ord_u, ord_v,
			      bu, bv,
			      newcoefs.begin(),
			      the_surface.dimension(),
			      the_surface.rational());

    if (elementary_surface_.get())
      the_subSurface->setElementarySurface(elementary_surface_);
    
    return the_subSurface;
}


//===========================================================================
vector<shared_ptr<ParamSurface> >
SplineSurface::subSurfaces(double from_upar,
			   double from_vpar,
			   double to_upar,
			   double to_vpar,
			   double fuzzy) const
//===========================================================================
{
    vector<shared_ptr<ParamSurface> > sub_sf;
    sub_sf.push_back(shared_ptr<ParamSurface>(subSurface(from_upar, from_vpar,
							 to_upar, to_vpar, fuzzy)));

    return sub_sf;
}


} // namespace Go;
