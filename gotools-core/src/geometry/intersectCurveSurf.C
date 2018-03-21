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
#include "GoTools/utils/Point.h"
#include "GoTools/geometry/SISLconversion.h"
#include "sislP.h"

#ifdef __BORLANDC__
using std::free;
#endif

//***************************************************************************
//
// Implementation file of the free function intersectCurveSurf defined in
// GoIntersections.h/
//
//***************************************************************************

using std::pair;
using std::vector;

namespace Go
{

 void intersectCurveSurf(const SplineCurve *cv, const SplineSurface *sf,
			 double epsge, 
			 vector<pair<double, Point> >& int_pts,
			 vector<int>& pretopology,
			 vector<pair<pair<double,Point>, 
			 pair<double,Point> > >& int_crvs)
 {
   SISLSurf* sislsf = GoSurf2SISL(*sf, false);
   SISLCurve* sislcv = Curve2SISL(*cv, false);
   int kntrack = 0;
   int trackflag = 0;  // Do not make tracks.
   SISLTrack **track =0;
   int knpt=0, kncrv=0;
   double *par1=0, *par2=0;
   int *pretop = 0;
   SISLIntcurve **vcrv = 0;
   int stat = 0;

   sh1858(sislsf, sislcv, 0.0, epsge, trackflag, &kntrack, &track,
	  &knpt, &par1, &par2, &pretop, &kncrv, &vcrv, &stat);
   ALWAYS_ERROR_IF(stat<0,"Error in intersection, code: " << stat);

   // Remember intersections points. 
   int ki;
   for (ki=0; ki<knpt; ki++)
     {
       int_pts.push_back(std::make_pair(par2[ki],
					Point(par1[2*ki],par1[2*ki+1])));
       pretopology.insert(pretopology.end(), pretop+4*ki+2, pretop+4*(ki+1));
       pretopology.insert(pretopology.end(), pretop+4*ki, pretop+4*ki+2);
     }

   // Remember intersection curves
   for (ki=0; ki<kncrv; ++ki)
     {
       int nmb_pt = vcrv[ki]->ipoint;
       Point par1 = Point(vcrv[ki]->epar1[0],vcrv[ki]->epar1[1]);
       Point par2 = Point(vcrv[ki]->epar1[2*(nmb_pt-1)],
			  vcrv[ki]->epar1[2*nmb_pt-1]);
       int_crvs.push_back(std::make_pair(std::make_pair(vcrv[ki]->epar2[0], 
							par1),
					 std::make_pair(vcrv[ki]->epar2[nmb_pt-1],
							par2)));
     }

   if (kncrv > 0)
     freeIntcrvlist(vcrv, kncrv);

   if (par1 != NULL) free(par1);
   if (par2 != NULL) free(par2);
   if (sislsf != NULL) freeSurf(sislsf);
   if (sislcv != NULL) freeCurve(sislcv);
   if (pretop != NULL) free(pretop);
 }

} // namespace Go  
