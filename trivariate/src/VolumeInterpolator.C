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

#include "GoTools/geometry/SplineInterpolator.h"
#include "GoTools/trivariate/VolumeInterpolator.h"
#include "GoTools/trivariate/SplineVolume.h"
#include "GoTools/geometry/SplineCurve.h"
#include "GoTools/geometry/SplineSurface.h"
//#include "sislP.h"

using std::vector;

namespace Go
{
  SplineVolume* 
  VolumeInterpolator::regularInterpolation(const BsplineBasis& basis_u,
					   const BsplineBasis& basis_v,
					   const BsplineBasis& basis_w,
					   vector<double>& par_u,
					   vector<double>& par_v,
					   vector<double>& par_w,
					   vector<double>& points,
					   int dimension,
					   bool rational,
					   vector<double>& weights)
  {
    // Check input
    ASSERT(par_u.size()*par_v.size()*par_w.size() == points.size()/dimension);
    ASSERT(basis_u.numCoefs() == (int)par_u.size());
    ASSERT(basis_v.numCoefs() == (int)par_v.size());
    ASSERT(basis_w.numCoefs() == (int)par_w.size());

    vector<double> points2;
    if (rational)
      {
	ASSERT(weights.size() == points.size()/dimension);

	// Include weight information
	// First get weights in interpolation points
	shared_ptr<SplineVolume> denom = 
	  shared_ptr<SplineVolume>(new SplineVolume(basis_u, basis_v,
						    basis_w,
						    weights.begin(), 1,
						    false));
	vector<double> wgtval;
	denom->gridEvaluator(par_u, par_v, par_w, wgtval);
	size_t nmb_pnt = par_u.size()*par_v.size()*par_w.size();
	points2.reserve(nmb_pnt*(dimension+1));
	for (size_t kr=0; kr<nmb_pnt; ++kr)
	  {
	    for (int kh=0; kh<dimension; kh++)
	      points2.push_back(points[kr*dimension+kh]*wgtval[kr]);
	    points2.push_back(wgtval[kr]);
	  }
	dimension++;
      }
    else
      points2 = points;

    // Interpolate surfaces in the second parameter direction and
    // curves in the first parameter direction
    size_t ki, kj;
    vector<double> sf_coefs;
    vector<int> tg_idx;
    vector<double> tg_pnt;
    for (kj=0; kj<par_w.size(); ++kj)
      {
	vector<double> cv_coefs;
	for (ki=0; ki<par_v.size(); ++ki)
	  {
	    // Interpolate
	    vector<double> coefs;
	    SplineInterpolator u_interpolator;
	    vector<double> pnts;
	    pnts.insert(pnts.end(), 
			points2.begin()+(kj*par_v.size()+ki)*dimension*par_u.size(),
			points2.begin()+(kj*par_v.size()+ki+1)*dimension*par_u.size());
	    u_interpolator.setBasis(basis_u);
	    u_interpolator.interpolate(par_u, pnts,
				       tg_idx, tg_pnt, coefs);
	    cv_coefs.insert(cv_coefs.end(), coefs.begin(), coefs.end());

// 	    vector<int> type(par_u.size(), 1);
// 	    vector<int> der(par_u.size(), 0);
// 	    int in;
// 	    int kstat = 0;
// 	    vector<double> knots;
// 	    double *coefs3;
// 	    knots.insert(knots.end(), basis_u.begin(), basis_u.end());
// 	    s1891(&par_u[0], 
// 		  &points2[(kj*par_v.size()+ki)*dimension*par_u.size()], 
// 		  dimension, par_u.size(), 1, &der[0], 1, &knots[0], &coefs3, 
// 		  &in, basis_u.order(), 0, 0, &kstat);
// 	    cv_coefs.insert(cv_coefs.end(), coefs3, coefs3+in*dimension);
// 	    free(coefs3);
	      
	  }

	// Interpolate the curves to make a surface
	SplineInterpolator v_interpolator;
	vector<double> coefs2;
	v_interpolator.setBasis(basis_v);
	v_interpolator.interpolate(par_v, cv_coefs, tg_idx, tg_pnt, coefs2);
	sf_coefs.insert(sf_coefs.end(), coefs2.begin(), coefs2.end());
      }

    // Interpolate surfaces to create volume
    SplineInterpolator w_interpolator;
    vector<double> vol_coefs;
    w_interpolator.setBasis(basis_w);
    w_interpolator.interpolate(par_w, sf_coefs, tg_idx, tg_pnt, vol_coefs);
    
    if (rational)
      {
	dimension--;
      }

    // Make surface
    SplineVolume* vol = new SplineVolume(basis_u, basis_v, basis_w,
					 vol_coefs.begin(), dimension,
					 rational);
    return vol;
  }

} // namespace Go
