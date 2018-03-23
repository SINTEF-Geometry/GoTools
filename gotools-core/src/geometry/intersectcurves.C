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

#include <vector>
using std::vector;
using std::pair;
#include "GoTools/geometry/GeometryTools.h"
#include "GoTools/geometry/SplineCurve.h"
#include "GoTools/geometry/GoIntersections.h"
#include "GoTools/geometry/SISLconversion.h"
#include "sislP.h"
 
//***************************************************************************
//
// Implementation file of the free function intersectcurves defined in
// GoIntersections.h/
//
//***************************************************************************

namespace Go
{

void intersectParamCurves(ParamCurve* cv1, ParamCurve* cv2, double epsge,
		     vector<std::pair<double,double> >& intersections)
  //************************************************************************
  // 
  // Intersect two curves. Collect intersection parameters.
  //
  //***********************************************************************
{
  shared_ptr<SplineCurve> spline_cv1(cv1->geometryCurve());
  shared_ptr<SplineCurve> spline_cv2(cv2->geometryCurve());
  return intersectcurves(spline_cv1.get(), spline_cv2.get(), epsge,
			 intersections);
}

void intersectcurves(SplineCurve* cv1, SplineCurve* cv2, double epsge,
		     vector<std::pair<double,double> >& intersections)
  //************************************************************************
  // 
  // Intersect two spline curves. Collect intersection parameters.
  //
  //***********************************************************************
{

  // Make sisl curves and call sisl.
  SISLCurve *pc1 = Curve2SISL(*cv1, false);
  SISLCurve *pc2 = Curve2SISL(*cv2, false);

  int kntrack = 0;
  int trackflag = 0;  // Do not make tracks.
  SISLTrack **track =0;
  int knpt=0, kncrv=0;
  double *par1=0, *par2=0;
  int *pretop = 0;
  SISLIntcurve **intcrvs = 0;
  int stat = 0;
  sh1857(pc1, pc2, 0.0, epsge, trackflag, &kntrack, &track,
	 &knpt, &par1, &par2, &pretop, &kncrv, &intcrvs, &stat);

  ALWAYS_ERROR_IF(stat<0,"Error in intersection, code: " << stat);


  // Remember intersections points. The intersection curves are
  // skipped.
  int ki;
  for (ki=0; ki<knpt; ki++)
    {
      intersections.push_back(std::make_pair(par1[ki],par2[ki]));
    }

  if (kncrv > 0)
    freeIntcrvlist(intcrvs, kncrv);

  if (par1 != 0) free(par1);
  if (par2 != 0) free(par2);
  if (pretop != 0) free(pretop);
  if (pc1 != 0) freeCurve(pc1);
  if (pc2 != 0) freeCurve(pc2);
}


void intersectParamCurves(ParamCurve* cv1, ParamCurve* cv2, double epsge,
			  vector<pair<double,double> >& intersections,
			  vector<pair<pair<double,double>,
			  pair<double,double> > >& int_cvs)
  //************************************************************************
  // 
  // Intersect two curves. Collect intersection parameters.
  //
  //***********************************************************************
{
  shared_ptr<SplineCurve> spline_cv1(cv1->geometryCurve());
  shared_ptr<SplineCurve> spline_cv2(cv2->geometryCurve());
  return intersectcurves(spline_cv1.get(), spline_cv2.get(), epsge,
			 intersections, int_cvs);
}

void intersectcurves(SplineCurve* cv1, SplineCurve* cv2, double epsge,
		     vector<std::pair<double,double> >& intersections,
		     vector<pair<pair<double,double>,
		     pair<double,double> > >& int_cvs)
  //************************************************************************
  // 
  // Intersect two spline curves. Collect intersection parameters.
  //
  //***********************************************************************
{

  // Make sisl curves and call sisl.
  SISLCurve *pc1 = Curve2SISL(*cv1, false);
  SISLCurve *pc2 = Curve2SISL(*cv2, false);

  int kntrack = 0;
  int trackflag = 0;  // Do not make tracks.
  SISLTrack **track =0;
  int knpt=0, kncrv=0;
  double *par1=0, *par2=0;
  int *pretop = 0;
  SISLIntcurve **intcrvs = 0;
  int stat = 0;
  sh1857(pc1, pc2, 0.0, epsge, trackflag, &kntrack, &track,
	 &knpt, &par1, &par2, &pretop, &kncrv, &intcrvs, &stat);

  ALWAYS_ERROR_IF(stat<0,"Error in intersection, code: " << stat);


  // Remember intersections points. 
  int ki;
  for (ki=0; ki<knpt; ki++)
    {
      intersections.push_back(std::make_pair(par1[ki],par2[ki]));
    }
  
  // Remember endpoints of intersection curves
  for (ki=0; ki<kncrv; ++ki)
    {
      int np = intcrvs[ki]->ipoint;
      int_cvs.push_back(std::make_pair(std::make_pair(intcrvs[ki]->epar1[0],
						      intcrvs[ki]->epar2[0]),
				       std::make_pair(intcrvs[ki]->epar1[np-1],
						      intcrvs[ki]->epar2[np-1])));
    }

  if (kncrv > 0)
    freeIntcrvlist(intcrvs, kncrv);

  if (par1 != 0) free(par1);
  if (par2 != 0) free(par2);
  if (pretop != 0) free(pretop);
  if (pc1 != 0) freeCurve(pc1);
  if (pc2 != 0) freeCurve(pc2);
}



} // namespace Go  
