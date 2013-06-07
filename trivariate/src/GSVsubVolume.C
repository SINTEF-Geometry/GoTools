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
#include "GoTools/trivariate/SplineVolume.h"
#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/trivariate/VolumeTools.h"

using std::vector;

namespace Go {


//===========================================================================
SplineVolume* SplineVolume::subVolume(double from_upar,
					 double from_vpar,
					 double from_wpar,
					 double to_upar,
					 double to_vpar,
					 double to_wpar,
					 double fuzzy) const
//===========================================================================
{
    if (from_upar >= to_upar) {
	THROW("First u-parameter must be strictly less than second.");
    }
    if (from_vpar >= to_vpar) {
	THROW("First v-parameter must be strictly less than second.");
    }
    if (from_wpar >= to_wpar) {
	THROW("First w-parameter must be strictly less than second.");
    }
    if (from_upar < startparam(0) || from_vpar < startparam(1) || from_wpar < startparam(2)) {
	THROW("Subsurface defined outside surface.");
    }

    // Currently turns the volume into curves in higher dimensional spaces,
    // and creates subcurves. Compare to subSurface on SplineSurface, where
    // this method is only used if one of the end parameter lies outside the
    // knot interval

    // First in the u direction.
    shared_ptr<SplineCurve> temp_cv = VolumeTools::representVolumeAsCurve(*this, 0);
    shared_ptr<SplineCurve> temp_sub(temp_cv->subCurve(from_upar, to_upar));
    shared_ptr<SplineVolume> vol
	= VolumeTools::representCurveAsVolume(*temp_sub, 0, basis_v_, basis_w_, rational());
    BsplineBasis new_ubasis = vol->basis(0);

    // Then in the v direction.
    temp_cv = VolumeTools::representVolumeAsCurve(*vol, 1);
    temp_sub.reset(temp_cv->subCurve(from_vpar, to_vpar));
    vol = VolumeTools::representCurveAsVolume(*temp_sub, 1, new_ubasis, basis_w_, rational());
    BsplineBasis new_vbasis = vol->basis(1);

    // Finally in the w direction.
    temp_cv = VolumeTools::representVolumeAsCurve(*vol, 2);
    temp_sub.reset(temp_cv->subCurve(from_wpar, to_wpar));
    vol = VolumeTools::representCurveAsVolume(*temp_sub, 2, new_ubasis, new_vbasis, rational());

    // We have to clone the return value because we cannot return a
    // shared pointer from this function.
    return vol->clone();
}


//===========================================================================
std::vector<shared_ptr<SplineVolume> > 
SplineVolume::split(std::vector<double>& param,
		    int pardir, double fuzzy) const
//===========================================================================
{
  // Represent the volume as a curve in the given parameter direction
  shared_ptr<SplineCurve> temp_cv = VolumeTools::representVolumeAsCurve(*this, pardir);

  // Split curve
  vector<shared_ptr<SplineCurve> > sub_cvs = temp_cv->split(param);

  BsplineBasis basis1 = (pardir > 0) ? basis_u_ : basis_v_;
  BsplineBasis basis2 = (pardir < 2) ? basis_w_ : basis_v_;

  // Represent the sub curves as volumes
  vector<shared_ptr<SplineVolume> > sub_vols;
  for (size_t ki=0; ki<sub_cvs.size(); ++ki)
    {
      shared_ptr<SplineVolume> vol = 
	  VolumeTools::representCurveAsVolume(*(sub_cvs[ki].get()), pardir, basis1, basis2, 
			       rational());
      sub_vols.push_back(vol);
    }

  return sub_vols;
}

} // namespace Go;
