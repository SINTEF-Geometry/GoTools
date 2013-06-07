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

#include "GoTools/trivariate/VolumeSpaceCurve.h"
#include "GoTools/trivariate/ParamVolume.h"
#include "GoTools/geometry/ParamCurve.h"
#include "GoTools/utils/Array.h"
#include "GoTools/utils/Volumes.h"

using namespace Go;
using std::vector;

//===========================================================================
VolumeSpaceCurve::VolumeSpaceCurve(shared_ptr<ParamVolume> vol, 
				   shared_ptr<ParamCurve> crv)
  : vol_(vol), crv_(crv)
//===========================================================================
{
}

//===========================================================================
VolumeSpaceCurve::~VolumeSpaceCurve()
//===========================================================================
{
}
 
//===========================================================================
Point VolumeSpaceCurve::eval(double t) const
//===========================================================================
{
  vector<Point> result(1);
  eval(t, 0, &result[0]);
  return result[0];
}

//===========================================================================
void VolumeSpaceCurve::eval(double t, int n, Point der[]) const
//===========================================================================
{
  if (n > 1)
    n = 1;
  if (n < 0)
    n = 0;

  // Evaluate curve
  vector<Point> cv_der(n+1);
  crv_->point(cv_der, t, n);

  // Evaluate volume
  vector<Point> vol_der(1+n*3);
  vol_->point(vol_der, cv_der[0][0], cv_der[0][1], cv_der[0][2], n);

  der[0] = vol_der[0];
  if (n == 1)
    {
      der[1] = cv_der[1][0]*vol_der[1] + cv_der[1][1]*vol_der[2] +
	cv_der[1][2]*vol_der[3];
    }
}

//===========================================================================
double VolumeSpaceCurve::start() const
//===========================================================================
{
  return crv_->startparam();
}

//===========================================================================
double VolumeSpaceCurve::end() const
//===========================================================================
{
  return crv_->endparam();
}

//===========================================================================
int VolumeSpaceCurve::dim() const
//===========================================================================
{
// Trivariate parameter area
  return 3;
}

//===========================================================================
bool VolumeSpaceCurve::approximationOK(double par, Point approxpos,
					   double tol1, double tol2) const
//===========================================================================
{
  Point pos = eval(par);
  double dist = pos.dist(approxpos);
  return (dist <= tol1);
}

