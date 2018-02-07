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

#include "GoTools/compositemodel/SISLCurveInterface.h"
#include "GoTools/geometry/SISLconversion.h"
#include "GoTools/utils/errormacros.h"
#include "sislP.h"

using std::vector;

using namespace Go;

//===========================================================================
shared_ptr<SplineCurve> SISLCurveInterface::interpolate(vector<Point>& pnts, 
							vector<int>& type,
							int degree)
//===========================================================================
{
  shared_ptr<SplineCurve> crv;
  if (pnts.size() < 2)
    return crv;  // Not enough information to produce a curve

  if (pnts.size() != type.size())
    THROW("Inconsistent interpolation information");

  // Translate input information
  int dim = pnts[0].dimension();
  vector<double> epoint;
  epoint.reserve(dim*pnts.size());

  vector<int> nptype(type.size());
  for (size_t ki=0; ki<pnts.size(); ++ki)
    {
      epoint.insert(epoint.end(), pnts[ki].begin(), pnts[ki].end());
    switch (type[ki])
      {
      case 1:
	nptype[ki] = 1;
	break;
      case 2:
	nptype[ki] = 4;
	break;
      case 3:
	nptype[ki] = 3;
	break;
      default:
	THROW("Unexpected interpolation condition");
      }
    }

  // Interpolate
  SISLCurve *qc = NULL;
  double *parvals = NULL;
  double startpar=0.0, endpar;
  int nmbpar;
  int status;
  s1356(&epoint[0], (int)pnts.size(), dim, &nptype[0], 0, 0, 1, degree+1,
	startpar, &endpar, &qc, &parvals, &nmbpar, &status);

  if (status < 0)
    {
      if (qc != NULL)
	freeCurve(qc);
      if (parvals != NULL)
	delete [](parvals);
      THROW("Error in SISL interpolation");
    }

  crv = shared_ptr<SplineCurve>(SISLCurve2Go(qc));

  if (qc != NULL)
    freeCurve(qc);
  if (parvals != NULL)
    delete [](parvals);
  
  return crv;
}

