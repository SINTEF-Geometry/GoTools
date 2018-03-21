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
#include "GoTools/geometry/SISLconversion.h"
#include "sislP.h"

#ifdef __BORLANDC__
using std::free;
#endif

//***************************************************************************
//
// Implementation file of the free function intersect2Dcurves defined in
// GoIntersections.h/
//
//***************************************************************************

using std::pair;
using std::vector;

namespace Go
{

void intersectCurvePoint(const ParamCurve* crv, Point pnt, double epsge,
			 vector<double>& intersections, 
			 vector<pair<double, double> >& int_crvs)
  //************************************************************************
  // 
  // Intersect a curve with a point
  //
  //***********************************************************************
{
  // First make sure that the curve is a spline curve
    ParamCurve* tmpcrv = const_cast<ParamCurve*>(crv);
    SplineCurve* sc = tmpcrv->geometryCurve();
    if (sc == NULL)
        THROW("ParamCurve doesn't have a spline representation.");
    
  // Make sisl curve and call sisl.
  SISLCurve *pc = Curve2SISL(*sc, false);

  int knpt=0, kncrv=0;
  double *par=0;
  SISLIntcurve **vcrv = 0;
  int stat = 0;
  s1871(pc, pnt.begin(), pnt.size(), epsge, &knpt, &par, &kncrv, &vcrv, &stat);
  ALWAYS_ERROR_IF(stat<0,"Error in intersection, code: " << stat);


  // Remember intersections points. 
  if (knpt > 0)
    intersections.insert(intersections.end(), par, par+knpt);

  // Remember intersection curves
  for (int ki=0; ki<kncrv; ++ki)
    int_crvs.push_back(std::make_pair(vcrv[ki]->epar1[0], vcrv[ki]->epar1[vcrv[ki]->ipoint-1]));

  if (kncrv > 0)
    freeIntcrvlist(vcrv, kncrv);

  if (par != 0) free(par);
  if (pc) freeCurve(pc);

  delete sc;
}

void intersect2Dcurves(const ParamCurve* cv1, const ParamCurve* cv2, double epsge,
		       vector<pair<double,double> >& intersections,
		       vector<int>& pretopology,
		       vector<pair<pair<double,double>, pair<double,double> > >& int_crvs)
  //************************************************************************
  // 
  // Intersect two 2D spline curves. Collect intersection parameters
  // and pretopology information.
  //
  //***********************************************************************
{

  // First make sure that the curves are spline curves.
    ParamCurve* tmpcv1 = const_cast<ParamCurve*>(cv1);
    ParamCurve* tmpcv2 = const_cast<ParamCurve*>(cv2);
    SplineCurve* sc1 = tmpcv1->geometryCurve();
    SplineCurve* sc2 = tmpcv2->geometryCurve();
    if (sc1 == NULL || sc2 == NULL)
        THROW("ParamCurves doesn't have a spline representation.");

    MESSAGE_IF(cv1->dimension() != 2,
		"Dimension different from 2, pretopology not reliable.");

  // Make sisl curves and call sisl.
  SISLCurve *pc1 = Curve2SISL(*sc1, false);
  SISLCurve *pc2 = Curve2SISL(*sc2, false);

  int kntrack = 0;
  int trackflag = 0;  // Do not make tracks.
  SISLTrack **track =0;
  int knpt=0, kncrv=0;
  double *par1=0, *par2=0;
  int *pretop = 0;
  SISLIntcurve **vcrv = 0;
  int stat = 0;
  sh1857(pc1, pc2, 0.0, epsge, trackflag, &kntrack, &track,
	 &knpt, &par1, &par2, &pretop, &kncrv, &vcrv, &stat);

  ALWAYS_ERROR_IF(stat<0,"Error in intersection, code: " << stat);


  // Remember intersections points. 
  int ki;
  for (ki=0; ki<knpt; ki++)
    {
      intersections.push_back(std::make_pair(par1[ki],par2[ki]));
      pretopology.insert(pretopology.end(), pretop+4*ki, pretop+4*(ki+1));
    }

  // Remember intersection curves
  for (ki=0; ki<kncrv; ++ki)
    int_crvs.push_back(std::make_pair(std::make_pair(vcrv[ki]->epar1[0], vcrv[ki]->epar2[0]), 
				      std::make_pair(vcrv[ki]->epar1[vcrv[ki]->ipoint-1],
						     vcrv[ki]->epar2[vcrv[ki]->ipoint-1])));

  if (kncrv > 0)
    freeIntcrvlist(vcrv, kncrv);

  if (par1 != 0) free(par1);
  if (par2 != 0) free(par2);
  if (pc1 != 0) freeCurve(pc1);
  if (pc2 != 0) freeCurve(pc2);
  if (pretop != 0) free(pretop);

  delete sc1;
  delete sc2;
}



} // namespace Go  
