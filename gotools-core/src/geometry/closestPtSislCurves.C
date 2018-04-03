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
 
//***************************************************************************
//
// Implementation file of the free function closestPtCurves defined in
// GoIntersections.h/
//
//***************************************************************************

namespace Go
{

void closestPtSislCurves(SplineCurve* cv1, SplineCurve* cv2, double epsge,
		     double& par1, double& par2, double& dist)
  //************************************************************************
  // 
  // Compute the closest point between two curves.
  //
  //***********************************************************************
{

  // Make guess point to the iteration
  // Find position of closest vertices
  std::vector<double>::const_iterator co1 = cv1->coefs_begin();
  std::vector<double>::const_iterator co2 = cv2->coefs_begin();
  std::vector<double>::const_iterator co3;
  std::vector<double>::const_iterator co12 = cv1->coefs_end();
  std::vector<double>::const_iterator co22 = cv2->coefs_end();
  int dim = cv1->dimension();
  double td, tmin=1.0e8;
  int minidx1=0, minidx2=0;
  int ki, k1, k2;
  for (k1=0; co1<co12; co1+=dim, k1++)
    for (k2=0, co3=co2; co3<co22; co3+=dim, k2++)
      {
	for (td=0.0, ki=0; ki<dim; ki++)
	  td += (co1[ki]-co3[ki])*(co1[ki]-co3[ki]);
	if (td < tmin)
	  {
	    tmin = td;
	    minidx1 = k1;
	    minidx2 = k2;
	  }
      }

  // Estimate parameter value of vertices
  std::vector<double>::const_iterator st;
  int kk = cv1->order();

  for (k1=minidx1+1, st=cv1->basis().begin(), par1=0.0;
       k1<minidx1+kk; par1+=st[k1], k1++);
  par1 /=(double)(kk-1);

  kk = cv2->order();
  for (k1=minidx2+1, st=cv2->basis().begin(), par2=0.0;
       k1<minidx2+kk; par2+=st[k1], k1++);
  par2 /=(double)(kk-1);

  // Make sisl curves and call sisl.
  SISLCurve *pc1 = Curve2SISL(*cv1, false);
  SISLCurve *pc2 = Curve2SISL(*cv2, false);

  // Iterate for closest point
  int stat = 0;
  s1770(pc1, pc2, epsge, cv1->startparam(), cv2->startparam(), 
	cv1->endparam(), cv2->endparam(), par1, par2, &par1, &par2, 
	&stat);
  ALWAYS_ERROR_IF(stat<0,"Error in closest point, code: " << stat);


  Point pt1, pt2;
  cv1->point(pt1, par1);
  cv2->point(pt2, par2);
  dist = pt1.dist(pt2);

  if (pc1 != 0) freeCurve(pc1);
  if (pc2 != 0) freeCurve(pc2);
}



} // namespace Go  
