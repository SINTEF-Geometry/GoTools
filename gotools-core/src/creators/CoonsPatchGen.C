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

#include "GoTools/creators/CoonsPatchGen.h"
#include "GoTools/geometry/SplineCurve.h"
#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/geometry/GeometryTools.h"
#include "GoTools/geometry/GoTools.h"
#include <algorithm>
#include <fstream> // for debugging
#include <iterator>

#ifdef __BORLANDC__
#include <iterator> // for back_inserter
#endif

using namespace Go;
using std::vector;
using std::max;
using std::back_inserter;

namespace {

SplineSurface* doCreatePatch(SplineCurve edge[]);

};// end anonymous namespace


//===========================================================================
SplineSurface*
Go::CoonsPatchGen::createCoonsPatch(const CurveLoop& boundary)
//===========================================================================
{
    // Check that the curve loop has four curves.
    ALWAYS_ERROR_IF(boundary.size() != 4, "Boundary must have four curves");


    // Extract and copy curves, checking that they are really SplineCurves.
    SplineCurve edge[4];
    int i;
    for (i = 0; i < 4; ++i) {
	SplineCurve* cv = dynamic_cast<SplineCurve*>(boundary[i].get());
	ALWAYS_ERROR_IF(cv == 0,
			"Curves must be of type SplineCurve.");

	// Assigning to a SplineCurve copies the contents
	edge[i] = *cv;
	// Make the edge curves k-regular
	edge[i].makeKnotStartRegular();
	edge[i].makeKnotEndRegular();
    }
   
    // Check that all curves have the same dimension, and that all curves
    // are nonrational.
    int dim = edge[0].dimension();
    for (i = 1; i < 4; ++i) {
	ALWAYS_ERROR_IF(edge[i].dimension() != dim,
			"Dimension mismatch.");

	// @@@ It may be possible to support rational surfaces
	// in the future.
	ALWAYS_ERROR_IF(edge[i].rational(),
			"Rational curves not supported.");

    }

    // Check that the orders of opposite curves are the same.
    // @@@ We could use order raising instead.
//      for (i = 0; i < 2; ++i) {
//  	ALWAYS_ERROR_IF(edge[i].order() != edge[i+2].order(),
//  		    "Order mismatch.",
//  		    InputError());
//      }
    return doCreatePatch(edge);
}




//===========================================================================
SplineSurface*
Go::CoonsPatchGen::createCoonsPatch(vector<shared_ptr<ParamCurve> >&
				      bd_curves,
				      vector<shared_ptr<ParamCurve> >&
				      cross_curves,
				      double epsge,
				      double kink_tol)
//===========================================================================
{
    double knot_diff_tol = GoTools::knotEpsilon(); // Default is 1e-05

    // Check that the bd_curves and cross_curves have four elements.
    ALWAYS_ERROR_IF((bd_curves.size() != 4) || (cross_curves.size() != 4),
		"Boundary must have four pos & cross tangent curves.");

    // We put all curves on a common vector.
    vector<shared_ptr<SplineCurve> > curves;
    for (int i = 0; i < 4; ++i) {
      shared_ptr<SplineCurve> cv1 = 
	    dynamic_pointer_cast<SplineCurve, ParamCurve>(bd_curves[i]);
      ALWAYS_ERROR_IF(cv1.get() == 0,
                  "Curves must be of type SplineCurve.");
      curves.push_back(cv1);
      curves.push_back(dynamic_pointer_cast<SplineCurve, ParamCurve>
		       (cross_curves[i])); // may be 0
    }

    // Check that all curves have the same dimension and are nonrational.
    int dim = curves[0]->dimension();
    for (size_t i = 0; i < curves.size(); ++i) {
	ALWAYS_ERROR_IF((curves[i].get() != 0) && (curves[i]->dimension() != dim),
		    "Dimension mismatch.");
	ALWAYS_ERROR_IF((curves[i].get() != 0) && curves[i]->rational(),
		    "Rational curves not supported.");
    }

    // Reparameterize curves
    const double tangent_ratio = 1/4.0;
    vector<shared_ptr<SplineCurve> > dummy_vec(1);
    for (int i=0; i<4; i++) {
	dummy_vec[0] = curves[2*i];
	CoonsPatchGen::reparamBoundaryCurve(dummy_vec, tangent_ratio);
    }

    // Let corresponding curves share parameter interval.
    double u_min = curves[0]->startparam();
    double u_max = u_min + max(curves[0]->estimatedCurveLength(),
			       curves[4]->estimatedCurveLength());
    double v_min = curves[2]->startparam();
    double v_max = v_min + max(curves[2]->estimatedCurveLength(),
			       curves[6]->estimatedCurveLength());
    int iind = 0;
    for (int i = 0, j = 0; j < 4; i=(i+4)%7, ++j) {
	if (curves[i].get())
	    curves[i]->setParameterInterval(u_min, u_max);
	if (curves[i+2].get())
	    curves[i+2]->setParameterInterval(v_min, v_max);
	iind = i;
    }

    // We modify existing cross_curves to satisfy tangent and twist requirements.
    vector<shared_ptr<SplineCurve> > mod_cross_curves;

#ifdef CREATORS_DEBUG
    {
	ofstream debug("data/debug.g2");
	for (size_t k = 0; k < bd_curves.size(); ++k) {
	    curves[2*k]->writeStandardHeader(debug);
	    curves[2*k]->write(debug);
	    if (curves[2*k+1].get() != 0) {
		try {
		    shared_ptr<SplineCurve> cv_sum(curveSum(*curves[2*k], 1.0,
							    *curves[2*k+1], 1.0));
		    cv_sum->writeStandardHeader(debug);
		    cv_sum->write(debug);
		} catch (...) {
		    ;
		}
	    }
	}
    }
#endif

    getCrossTangs(curves, mod_cross_curves, epsge, kink_tol);

#ifdef CREATORS_DEBUG
    {
	ofstream debug("data/debug.g2");
	for (size_t k = 0; k < bd_curves.size(); ++k) {
	    curves[2*k]->writeStandardHeader(debug);
	    curves[2*k]->write(debug);
	    if (mod_cross_curves[iind].get() != 0) {
		try {
		    shared_ptr<SplineCurve>
			cv_sum(curveSum(*curves[2*k], 1.0,
					*mod_cross_curves[k], 1.0));
		    cv_sum->writeStandardHeader(debug);
		    cv_sum->write(debug);
		} catch (...) {
		    ;
		}
	    }
	}
    }
#endif

    // Extract boundary curves to a separate vector.
    vector<shared_ptr<SplineCurve> > bound_curves;
    for (int i = 0; i < 4; ++i)
	bound_curves.push_back(curves[2*i]);

    // We make sure that corr. curves live in the same spline space.
    // As curves are given in order, we reverse parameterdirection for side 2 & 3.
    for (int i = 2; i < 4; ++i) {
	bound_curves[i]->reverseParameterDirection();
	if (mod_cross_curves[i].get() != 0)
	    mod_cross_curves[i]->reverseParameterDirection();
    }
    vector<shared_ptr<SplineCurve> > dummy_vector_u(4), dummy_vector_v(4);
    for (int i = 0; i < 2; ++i) {
	dummy_vector_u[i*2] = bound_curves[i*2];
	dummy_vector_u[i*2+1] = mod_cross_curves[i*2];
	dummy_vector_v[i*2] = bound_curves[i*2+1];
	dummy_vector_v[i*2+1] = mod_cross_curves[i*2+1];
    }
    GeometryTools::unifyCurveSplineSpace(dummy_vector_u, knot_diff_tol);
    GeometryTools::unifyCurveSplineSpace(dummy_vector_v, knot_diff_tol);
    // As objects may have changed, we must extract the curves.
    for (int i = 0; i < 2; ++i) {
	bound_curves[i*2] = dummy_vector_u[i*2];
	mod_cross_curves[i*2] = dummy_vector_u[i*2+1];
	bound_curves[i*2+1] = dummy_vector_v[i*2];
        mod_cross_curves[i*2+1] = dummy_vector_v[i*2+1];
    }
    // As curves are to be given in order, reverse parameterdirection for side 2 & 3.
    for (int i = 2; i < 4; ++i) {
	bound_curves[i]->reverseParameterDirection();
	if (mod_cross_curves[i].get() != 0)
	    mod_cross_curves[i]->reverseParameterDirection();
    }

    return createCoonsPatch(bound_curves, mod_cross_curves);
}


// All input data are assumed to exist and satisfy requirements.
// bd_cv[] should form a loop, and cross_cv[] should have the same parametrization.
// Missing cross curves are indicated by a null pointer.
//===========================================================================
SplineSurface*
Go::CoonsPatchGen::createCoonsPatch(vector<shared_ptr<SplineCurve> >&
				      bd_curves,
				      vector<shared_ptr<SplineCurve> >&
				      cross_curves)   
//===========================================================================
{
//     // TEST INPUT
//     if (true) {
// 	cout << "Evaluating the corners of coons patch. " << endl;
// 	vector<Point> pt(2);
// 	for (i=0; i<4; i++) {
// 	    bd_curves[i]->point(pt, bd_curves[i]->startparam(), 1);
// 	    pt[0].write(cout);
// 	    pt[1].write(cout);
// 	    if (cross_curves[i].get() != 0) {
// 		cross_curves[i]->point(pt, cross_curves[i]->startparam(), 1);
// 		pt[0].write(cout);
// 		pt[1].write(cout);
// 	    }
// 	    bd_curves[i]->point(pt, bd_curves[i]->endparam(), 1);
// 	    pt[0].write(cout);
// 	    pt[1].write(cout);
// 	    if (cross_curves[i].get() != 0)  {
// 		cross_curves[i]->point(pt, cross_curves[i]->endparam(), 1);
// 		pt[0].write(cout);
// 		pt[1].write(cout);
// 	    }
// 	    cout << " " << endl;
// 	}
// 	cout << "End corner evaluation." << endl;
//     }

    double knot_diff_tol = 1e-05;

    // The curves are given in order, hence we reverse parameterdirection for 2 & 3.
    bd_curves[2]->reverseParameterDirection();
    bd_curves[3]->reverseParameterDirection();
    if (cross_curves[2].get() != 0)
	cross_curves[2]->reverseParameterDirection();
    if (cross_curves[3].get() != 0)
     cross_curves[3]->reverseParameterDirection();

    // As all cross_curves point inwards, we must turn two (1 & 2).
    for (int i = 1; i < 3; ++i) {
	if (cross_curves[i].get() == 0) continue;
   	for (vector<double>::iterator iter = cross_curves[i]->coefs_begin();
    	     iter != cross_curves[i]->coefs_end(); ++iter)
	    iter[0] *= -1.0;
    }

    vector<double> params;
    params.push_back(bd_curves[1]->startparam());
    params.push_back(bd_curves[1]->endparam());
    params.push_back(bd_curves[0]->startparam());
    params.push_back(bd_curves[0]->endparam());
    vector<int> cross_index_u, cross_index_v, cross_index;
    vector<shared_ptr<SplineCurve> > mesh_curves, bd_cross_curves;
    for (int i = 0; i < 2; ++i)
	for (int j = 0; j < 2; ++j) {
	    // The loft and tp-routines assume u-curves first, then v-curves.
	    mesh_curves.push_back(bd_curves[(3*i+2*j)%4]);
	    if (cross_curves[(3*i+2*j)%4].get() != 0) {
		bd_cross_curves.push_back(cross_curves[(3*i+2*j)%4]);
		if (i == 0)
		    cross_index_u.push_back((int)mesh_curves.size() - 1);
		else
		    cross_index_v.push_back((int)mesh_curves.size() - 1);
	    }
	}

    set_union(cross_index_u.begin(), cross_index_u.end(),
	      cross_index_v.begin(), cross_index_v.end(),
	      std::back_inserter(cross_index));

    // Lofting in v-direction. loft_u_sf indicates we're lofting u-curves.
    shared_ptr<SplineSurface> loft_u_sf =
      shared_ptr<SplineSurface>(loftSurface(mesh_curves.begin(),
					      params.begin(), 2,
					      bd_cross_curves.begin(),
					      cross_index_u));
    // Creating the lofted surface in the u-direction, we must bear in mind
    // that the surface created needs to swap parameter directions.
    for (size_t i = 0; i < cross_index_v.size(); ++i)
	cross_index_v[i] -= 2;
    shared_ptr<SplineSurface> loft_v_sf =
      shared_ptr<SplineSurface>(loftSurface(mesh_curves.begin() + 2,
					      params.begin() + 2, 2,
					      bd_cross_curves.begin() +
					      cross_index_u.size(),
					      cross_index_v));
    loft_v_sf->swapParameterDirection();

    shared_ptr<SplineSurface> tp_sf =
      shared_ptr<SplineSurface>(tpSurface(mesh_curves, params, 2,
					    bd_cross_curves, cross_index));

    vector<shared_ptr<SplineSurface> > surfaces;
    surfaces.push_back(loft_u_sf);
    surfaces.push_back(loft_v_sf);
    surfaces.push_back(tp_sf);
    GeometryTools::unifySurfaceSplineSpace(surfaces, knot_diff_tol);

    for (int i = 0; i < surfaces[0]->dimension() *
	     surfaces[0]->numCoefs_u() * surfaces[0]->numCoefs_v(); ++i) {
	surfaces[0]->coefs_begin()[i] += surfaces[1]->coefs_begin()[i];
	surfaces[0]->coefs_begin()[i] -= surfaces[2]->coefs_begin()[i];
    }

    // As surfaces[0] is controlled by a smart pointer, we must make a copy.
#if ((_MSC_VER > 0) && (_MSC_VER < 1300))
    SplineSurface* return_sf = static_cast<SplineSurface*>(surfaces[0]->clone());
#else
    SplineSurface* return_sf = surfaces[0]->clone();
#endif


//     // TESTING
//       double mindist, meddist, maxdist;
//       double minang, medang, maxang;
//       SplineCurve *poscv, *dercv1, *dercv2;
//       bool udir;
//       double par;
//       int dir, ki;
//       for (ki=0; ki<4; ki++)
// 	{
// 	  if (cross_curves[ki].get() == 0)
// 	    continue;
// 	  udir = (ki == 0 || ki == 2);
// 	  dir = (ki == 0 || ki == 2) ? 1 : 0;
// 	  par = (ki==0 || ki==3) ? return_sf.basis(dir).startparam() 
// 	    : return_sf.basis(dir).endparam();
// 	  return_sf.constParamCurve(par, udir, poscv, dercv2);
// 	  dercv1 = poscv->derivCurve(1);
	      
// 	  checkContinuity(poscv, dercv1, 
// 			  dercv2, cross_curves[ki].get(), 
// 			  minang, medang, maxang);
// 	  delete poscv;
// 	  delete dercv1;
// 	  delete dercv2;
// 	  cout << "Boundary no " << ki << endl;
// 	  cout << "minang : " << minang << ", medang : " << medang;
// 	  cout << ", maxang : " << maxang << endl;
// 	}

    // That should be it. Returning our promised Coons patch.
    return return_sf;
}

namespace {

//===========================================================================
SplineSurface* doCreatePatch(SplineCurve edge[])
//===========================================================================
{
    // We need the corner points
    int i;
    Point corner[4];
    for (i = 0; i < 4; ++i) {
	edge[i].point(corner[i], edge[i].startparam());
    }

    // Turn curves 2 and 3 so that they match 0 and 1
    edge[2].reverseParameterDirection();
    edge[3].reverseParameterDirection();

    // Give all curves a knot interval of [0,1]
    double av[2];
    for (i = 0; i < 2; i+=1) {
      // double len1 = edge[i].estimatedCurveLength();
      // double len2 = edge[i+2].estimatedCurveLength();
      double len1 = edge[i].endparam() - edge[i].startparam();
      double len2 = edge[i+2].endparam() - edge[i+2].startparam();
      av[i] = 0.5*(len1+len2);
      edge[i].setParameterInterval(0.0, 1.0);
      edge[i+2].setParameterInterval(0.0, 1.0);
    }

    // Give opposing curves the same order
    for (i = 0; i < 2; ++i){
      int maxorder = max(edge[i].order(), edge[i+2].order());
      if (edge[i].order() < maxorder)
	edge[i].raiseOrder(maxorder-edge[i].order());
      if (edge[i+2].order() < maxorder)
	edge[i+2].raiseOrder(maxorder-edge[i+2].order());
    }
      

    // Give opposing curves the same knotvectors
    vector<shared_ptr<SplineCurve> > u_curves, v_curves;
    for (i = 0; i < 2; ++i) {
	u_curves.push_back(shared_ptr<SplineCurve>
			   (new SplineCurve(edge[2*i].numCoefs(),
					      edge[2*i].order(),
					      edge[2*i].basis().begin(),
					      edge[2*i].coefs_begin(),
					      edge[2*i].dimension())));
	v_curves.push_back(shared_ptr<SplineCurve>
			   (new SplineCurve(edge[2*i+1].numCoefs(),
					      edge[2*i+1].order(),
					      edge[2*i+1].basis().begin(),
					      edge[2*i+1].coefs_begin(),
					      edge[2*i+1].dimension())));
    }
    double knot_diff_tol = 1e-05;
    GeometryTools::unifyCurveSplineSpace(u_curves, knot_diff_tol);
    GeometryTools::unifyCurveSplineSpace(v_curves, knot_diff_tol);
    for (i = 0; i < 2; ++i) {
	edge[i*2] = SplineCurve(u_curves[i]->numCoefs(), u_curves[i]->order(),
				  u_curves[i]->basis().begin(),
				  u_curves[i]->coefs_begin(),
				  u_curves[i]->dimension());
	edge[i*2+1] = SplineCurve(v_curves[i]->numCoefs(), v_curves[i]->order(),
				    v_curves[i]->basis().begin(),
				    v_curves[i]->coefs_begin(),
				    v_curves[i]->dimension());
    }
//     vector<double> new_knots[4];
//     for (i = 0; i < 2; ++i) {
// 	set_difference(edge[i].basis().begin(), edge[i].basis().end(),
// 		       edge[i+2].basis().begin(), edge[i+2].basis().end(),
// 		       back_inserter(new_knots[i+2]));
// 	set_difference(edge[i+2].basis().begin(), edge[i+2].basis().end(),
// 		       edge[i].basis().begin(), edge[i].basis().end(),
// 		       back_inserter(new_knots[i]));
//     }
//     for (i = 0; i < 4; ++i) {
// 	edge[i].insertKnot(new_knots[i]);
//     }
#ifdef DEBUG
    std::ofstream of("coonsgen.g2");
    for (i=0; i<2; ++i)
      {
	u_curves[i]->writeStandardHeader(of);
	u_curves[i]->write(of);
      }
    for (i=0; i<2; ++i)
      {
	v_curves[i]->writeStandardHeader(of);
	v_curves[i]->write(of);
      }
#endif


    // Store the number of coefficients and the orders of
    // both parameter directions.
    int dim = edge[0].dimension();
    int n0 = edge[0].numCoefs();
    int n1 = edge[1].numCoefs();
    int k0 = edge[0].order();
    int k1 = edge[1].order();
    const BsplineBasis& basu = edge[0].basis();
    const BsplineBasis& basv = edge[1].basis();

    // Make coefficients for three surfaces:
    // co[0] is interpolated in the first (u) parameter direction
    // co[1] is interpolated in the second (v) parameter direction
    // co[2] is the correction surface (to be subtracted from)
    vector<double> co[3];

    vector<double> knu(2*k0, 0.0);
    vector<double> knv(2*k1, 0.0);
    fill(knu.begin()+k0, knu.end(), 1.0);
    fill(knv.begin()+k1, knv.end(), 1.0);
    BsplineBasis minbasu(k0, k0, &(knu[0]));
    BsplineBasis minbasv(k1, k1, &(knv[0]));

    co[0].resize(k0*n1*dim);
    for (i = 0; i < k0; ++i) {
	double fac = minbasu.grevilleParameter(i);
	for (int j = 0; j < n1; ++j) {
	    for (int dd = 0; dd < dim; ++dd) {
		co[0][j*k0*dim + i*dim + dd]
		    = (1.0-fac)*(edge[3].coefs_begin()[j*dim+dd])
		    + fac*(edge[1].coefs_begin()[j*dim+dd]);
	    }
	}
    }

    co[1].resize(k1*n0*dim);
    for (i = 0; i < k1; ++i) {
	double fac = minbasv.grevilleParameter(i);
	for (int dd = 0; dd < n0*dim; ++dd) {
	    co[1][i*n0*dim + dd]
		= (1.0-fac)*(edge[0].coefs_begin()[dd])
		+ fac*(edge[2].coefs_begin()[dd]);
	}
    }

    co[2].resize(k0*k1*dim);
    for (i = 0; i < k0; ++i) {
	double fac0 = minbasu.grevilleParameter(i);
	for (int j = 0; j < k1; ++j) {
	    double fac1 = minbasv.grevilleParameter(j);
	    for (int dd = 0; dd < dim; ++dd) {
		co[2][(j*k0+i)*dim + dd]
		    = (1.0-fac0)*(1.0-fac1)*corner[0][dd]
		    + fac0*(1.0-fac1)*corner[1][dd]
		    + fac0*fac1*corner[2][dd]
		    + (1.0-fac0)*fac1*corner[3][dd];
	    }
	}
    }

    // Put all surfaces on the same knot vectors
    vector<double> newu, newv;
    set_difference(basu.begin(), basu.end(),
			minbasu.begin(), minbasu.end(),
			std::back_inserter(newu));
    set_difference(basv.begin(), basv.end(),
			minbasv.begin(), minbasv.end(),
			std::back_inserter(newv));
    SplineSurface s1(minbasu, basv,    co[0].begin(), dim);
    SplineSurface s2(basu,    minbasv, co[1].begin(), dim);
    SplineSurface s3(minbasu, minbasv, co[2].begin(), dim);


    s1.insertKnot_u(newu);
    s2.insertKnot_v(newv);
    s3.insertKnot_u(newu);
    s3.insertKnot_v(newv);

//      s1.writeStandardHeader(cout);
//      s1.write(cout);

    // Now we add s2 and subtract s3 from s1
    for (i = 0; i < n0*n1*dim; ++i) {
	s1.coefs_begin()[i] += s2.coefs_begin()[i];
	s1.coefs_begin()[i] -= s3.coefs_begin()[i];
    }
    s1.setParameterDomain(0.0, av[0], 0.0, av[1]);
    //s1.setParameterDomain(0.0, 1.0, 0.0, 1.0);
#ifdef DEBUG
    s1.writeStandardHeader(of);
    s1.write(of);
#endif

#if ((_MSC_VER > 0) && (_MSC_VER < 1300))
    return dynamic_cast<SplineSurface*>(s1.clone());
    
#else
    return s1.clone();
#endif
}

}; // end anonymous namespace
