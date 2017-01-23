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

#include "GoTools/creators/SmoothTransition.h"

#include "GoTools/creators/CreatorsOffsetUtils.h"

#include "GoTools/geometry/LineCloud.h"

#include <fstream>

using namespace Go;
using std::vector;
using std::max;
using std::min;


void SmoothTransition::
offsetIntersectionIterate(double arad1, double arad2, std::vector<Point>& epoint,
			  std::vector<Point>& epnt1, std::vector<Point>& epnt2,
			  Point& epar1, Point& epar2,
			  const SplineSurface& psurf1,
			  const SplineSurface& psurf2,
			  double astep, double aepsge, std::vector<Point>& gpnt1,
			  std::vector<Point>& gpnt2, std::vector<Point>& goffpnt1,
			  std::vector<Point>& goffpnt2, Point& gpar1,
			  Point& gpar2)
     /*
*********************************************************************
*                                                                   
* PURPOSE    : To iterate to an intersection point between two surfaces
*              and a plane.
*
*
*
* INPUT      : arad1  - First offset distance.
*              arad2  - Second offset distance.
*              epoint - Array containing parts of plane description.
*                       epoint[0:2] contains a position value.
*                       epoint[3:5] contains the normal to the plane
*                       A point in the plane is defined by
*                       epoint[0:2] + astep*epoint[3:5]
*              epnt1  - 0-2 Derivatives + normal of start point for
*                       iteration in first surface
*              epnt2  - 0-2 Derivatives + normal of start point for
*                       iteration in second surface
*              epar1  - Parameter pair of start point in first surface
*              epar2  - Parameter pair of start point in second surface
*              psurf1 - Description of first surface
*              psurf2 - Description of second surface
*              astep  - Step length
*              aepsge - Absolute tolerance
*
*
* OUTPUT     : gpnt1  - 0-2 Derivatives + normal of result of iteration
*                       in first surface
*              gpnt2  - 0-2 Derivatives + normal of result of iteration
*                       in second surface
*              goffpnt1 - 0-2 Derivatives + normal of result of iteration
*                       in first offset surface.
*              goffpnt2 - 0-2 Derivatives + normal of result of iteration
*                       in second offset surface.
*              gpar1  - Parameter pair of result of iteration in first surface
*              gpar2  - Parameter pair of result of iteration in second
*                       surface
*              jstat  - status messages  
*                       = 2      : Iteration diverged or to many iterations
*                       = 1      : iteration converged, singular point found
*                       = 0      : ok, iteration converged
*                       < 0      : error
*
*
* METHOD     :
*
* USE        : The function is only working i 3-D
*
* REFERENCES :
*
*-
* CALLS      :
*
* WRITTEN BY : Tor Dokken, SI, Oslo, Norway, June-1988
* Revised by : Tor Dokken, SI, OSLO, Norway, 24-Feb-1989
*              Prepared for degenerate points
* Revised by : Tor Dokken, SI, Oslo, Norway, 3-April-1989
*              Correct handling of small determinats
*
*********************************************************************
*/
{  

//     // debug
//     std::ofstream of2("data/debug.g2");
//     vector<double> line_pts1;
//     line_pts1.insert(line_pts1.end(), epoint[0].begin(), epoint[0].end());
//     Point to_pt = epoint[0] + epoint[1];
//     line_pts1.insert(line_pts1.end(), to_pt.begin(), to_pt.end());
//     LineCloud line_cloud1(&line_pts1[0], 1);
//     line_cloud1.writeStandardHeader(of2);
//     line_cloud1.write(of2);
//     // end debug

  int kcont;              /* Indicator telling if iteration is not finished */
  int kder = 2;           /* Derivative indicator                           */
  int klfu=0;             /* Pointer into knot vector                       */
  int klfv=0;             /* Pointer into knot vector                       */
  int klfs=0;             /* Pointer into knot vector                       */
  int klft=0;             /* Pointer into knot vector                       */
  int kstat;              /* Status variable                                */
  int knbit;              /* Counter for number of iterations               */
  int kmaxit = 100;       /* Maximal number of iterations allowed           */
//   int kpos=1;             /* Position indicator ofr errors                  */
  Point spoint(3);       /* Point in intersection plane                    */
  vector<Point>::iterator snorm; /* Pointer to normal vector of intersection plane */
  vector<Point>::iterator sp, spu, spv, spn; /* Pointers into goffpnt1                      */
  vector<Point>::iterator sq, sqs, sqt, sqn; /* Pointers into goffpnt2                         */
  double ta11,ta12,ta21;  /* Variables used in equation systems             */
  double ta22,tb1,tb2;    /* Variables used in equation systems             */
  Point sdiff(3);       /* Difference between two vectors                 */
  double tdum1,tdum2;     /* Dummy variables                                */
  double tdum3,tdum;      /* Dummy variables                                */
  double tdist = 1e100;   /* Distance betweentwo points in iteration        */

  /* Make description of intersection plane */
  spoint = epoint[0] + astep*epoint[1];

//   // debug
//   std::ofstream of3("data/debug.g2");
//   Point normal_pt = epoint[0] + 1.0*epoint[1];
//   SplineCurve plane_normal(epoint[0], normal_pt);
//   plane_normal.writeStandardHeader(of3);
//   plane_normal.write(of3);
//   // end debug

  snorm = epoint.begin() + 1;

  /* Copy input variables to output variables */
  goffpnt1 = epnt1;
  goffpnt2 = epnt2;
  gpar1 = epar1;
  gpar2 = epar2;
  
  /* At the start of the iteration the two point goffpnt1 and goffpnt2 might be
     very close since we in most cases start from a point on the intersection
     curve. */
  
  /* Set a number of local pointers that are used often */
  sp = goffpnt1.begin();
  spu = goffpnt1.begin() + 1;
  spv = goffpnt1.begin() + 2;
  spn = goffpnt1.begin() + 6;
  sq = goffpnt2.begin();
  sqs = goffpnt2.begin() + 1;
  sqt = goffpnt2.begin() + 2;
  sqn = goffpnt2.begin() + 6;
  
  kcont = 1;
  knbit = 0;
  
  while (kcont) {
      
      /* Put a parametric representation of the tangent 
	 plane of surface 1 into
	 the implicit representation of the tangent 
	 plane of surface 2 and also
	 into the implicit representation of 
	 the intersection plane */

      ta11 = (*spu)*(*sqn);
      ta12 = (*spv)*(*sqn);
      ta21 = (*spu)*(*snorm);
      ta22 = (*spv)*(*snorm);

      sdiff = *sq - *sp;
      tb1 = sdiff*(*sqn);
      
      tdum = max(fabs(ta11),fabs(ta12));
      tdum = max(tdum,fabs(tb1));
      //       if (tdum == DNULL) tdum = (double)1.0;
      if (tdum == 0.0)
	  tdum = 1.0;
      ta11 /= tdum;
      ta12 /= tdum;
      tb1  /= tdum;
      
      sdiff = spoint - *sp;
      tb2 = sdiff*(*snorm);
      
      tdum = max(fabs(ta21),fabs(ta22));
      tdum = max(tdum,fabs(tb2));
      if (tdum == 0.0)
	  tdum = 1.0;
      ta21 /= tdum;
      ta22 /= tdum;
      tb2  /= tdum;
      
      /* Calculate determinant of equation system */
      
      tdum1 = ta11*ta22 - ta12*ta21;
      tdum  = max(fabs(ta11),fabs(ta22));
      tdum  = max(fabs(ta12),tdum);
      tdum  = max(fabs(ta21),tdum);
      
      //       if (DEQUAL((tdum+tdum1),tdum)) tdum1 =DNULL;
      if (tdum + tdum1 == tdum)
	  tdum1 = 0.0;

      /* If tdum1 = 0.0, then the equation system is singular, 
	 iteration not possible */

      //       if (DNEQUAL(tdum1,DNULL))
      if (tdum != 0.0) {
	  gpar1[0] += (tb1*ta22-tb2*ta12)/tdum1;
	  gpar1[1] += (ta11*tb2-ta21*tb1)/tdum1;

	  if (gpar1[0] < psurf1.startparam_u())
	      gpar1[0] = psurf1.startparam_u();
	  else if (gpar1[0] > psurf1.endparam_u())
	      gpar1[0] = psurf1.endparam_u();
	  if (gpar1[1] < psurf1.startparam_v())
	      gpar1[1] = psurf1.startparam_v();
	  else if (gpar1[1] > psurf1.endparam_v())
	      gpar1[1] = psurf1.endparam_v();
      }
      
      /* Put a parametric representation of the 
	 tangent plane of surface 2 into
	 the implicit representation of the 
	 tangent plane of surface 1 and also
	 into the implicit representation 
	 of the intersection plane */

      ta11 = (*sqs)*(*spn);
      ta12 = (*sqt)*(*spn);
      ta21 = (*sqs)*(*snorm);
      ta22 = (*sqt)*(*snorm);
      
      sdiff = *sp - *sq;
      tb1 = sdiff*(*spn);
      sdiff = spoint - *sq;
      tb2 = sdiff*(*snorm);
      
      /*Calculate determinant of equation system */

      tdum2 = ta11*ta22 - ta12*ta21;
      
      tdum2 = ta11*ta22 - ta12*ta21;
      tdum  = max(fabs(ta11),fabs(ta22));
      tdum  = max(fabs(ta12),tdum);
      tdum  = max(fabs(ta21),tdum);
      
      //       if (DEQUAL((tdum+tdum2),tdum)) tdum2 =DNULL;
      if (tdum + tdum2 == tdum)
	  tdum2 = 0.0;
      
      /* If tdum2 = 0.0, then the equation system is singular, 
	 iteration not possible */

      //       if (DNEQUAL(tdum2,DNULL))
      if (tdum2 != 0.0) {
	  gpar2[0] += (tb1*ta22-tb2*ta12)/tdum2;
	  gpar2[1] += (ta11*tb2-ta21*tb1)/tdum2;

	  if (gpar2[0] < psurf2.startparam_u())
	      gpar2[0] = psurf2.startparam_u();
	  else if (gpar2[0] > psurf2.endparam_u())
	      gpar2[0] = psurf2.endparam_u();
	  if (gpar2[1] < psurf2.startparam_v())
	      gpar2[1] = psurf2.startparam_v();
	  else if (gpar2[1] > psurf2.endparam_v())
	      gpar2[1] = psurf2.endparam_v();
      }
      
      /* Calculate values of new points */
      
      //       blend_s1421(psurf1,arad1,kder,gpar1,&klfu,&klfv,goffpnt1,goffpnt1+18,
      // 		  gpnt1,gpnt1+18,&kstat);
      // @@sbr Really not using klfu & klfv.
      OffsetUtils::blend_s1421(&psurf1, arad1, kder, gpar1, klfu, klfv,
                               goffpnt1, gpnt1, &kstat);
      ALWAYS_ERROR_IF(kstat < 0,
		  "Method failed.");
      
      /* If the surface normal has zero length no use in continuing */
      if (kstat == 2) {
	  MESSAGE("Surface normal has zero length.");
	  return;
      }
      
      OffsetUtils::blend_s1421(&psurf2, arad2, kder, gpar2, klfs, klft,
                               goffpnt2, gpnt2, &kstat);
      ALWAYS_ERROR_IF(kstat < 0,
		  "Method failed.");

//       // debug
//       std::ofstream of("data/debug.g2");
//       Point from_1 = psurf1.ParamSurface::point(gpar1[0], gpar1[1]);
//       SplineCurve normal_cv1(from_1, goffpnt1[0]);
//       Point from_1_orig = psurf1.ParamSurface::point(epar1[0], epar1[1]);
//       SplineCurve normal_cv1_orig(from_1_orig, goffpnt1[0]);
//       Point from_2 = psurf2.ParamSurface::point(gpar2[0], gpar2[1]);
//       SplineCurve normal_cv2(from_2, goffpnt2[0]);
//       Point from_2_orig = psurf2.ParamSurface::point(epar2[0], epar2[1]);
//       SplineCurve normal_cv2_orig(from_2_orig, goffpnt2[0]);
//       normal_cv1.writeStandardHeader(of);
//       normal_cv1.write(of);
//       normal_cv2.writeStandardHeader(of);
//       normal_cv2.write(of);
//       normal_cv1_orig.writeStandardHeader(of);
//       normal_cv1_orig.write(of);
//       normal_cv2_orig.writeStandardHeader(of);
//       normal_cv2_orig.write(of);
//       // end debug

      
      /* If the surface normal has zero length no use in continuing */
      if (kstat == 2) {
	  MESSAGE("Surface normal has zero length.");
	  return;
      }
      
      /* Make difference between the two points, 
	 and calculate length of difference */
      sdiff = goffpnt1[0] - goffpnt2[0];
      tdum3 = sdiff.length();
      if (tdum3 == 0.0) {
	  /* Length is zero iteration has converged */
	  kcont = 0;
      }
      
      if (knbit==0) {
	  /* First iteration intitate distance variable, if the equation
	     systems were not singular */
	  // 	  if (DEQUAL(tdum1,DNULL) || DEQUAL(tdum2,DNULL)) goto war02;
	  if ((tdum1 == 0.0) || (tdum2 == 0.0)) {
	      MESSAGE("This should not happen.");
	      return;
	  }
	  tdist = tdum3;
	  knbit = 1;
      } else {
	  /* More than one iteration done, stop if distance is not decreasing.
	     Then decide if we converge distance between the points is within
	     the tolerance and the last step had singular or none singular
	     equation systems. */

	  knbit = knbit + 1;
	  if (tdum3>=tdist) {
	      /* Distance is not decreasing */
	      if (tdist <= aepsge) {
		  /* Distance within tolerance */
		  // 		  if (DEQUAL(tdum1,DNULL) || DEQUAL(tdum2,DNULL))
		  if ((tdum1 == 0.0) || (tdum2 == 0.0)) { // @@sbr Tolerance?
		      MESSAGE("Iteration converged, singular point found");
		      return;
		  } else {
		      /* Nonsingular equation system */
		      return;
		  }
	      } else {
		  MESSAGE("Distance is not within tolerance, divergence");
		  return;
	      }
	  }
	  /*      Distance still decreasing */
	  
	  tdist = tdum3;
      }
      
      /*  Make sure that not to many iteration are being done */
      if (knbit > kmaxit) {
	  MESSAGE("Iteration count exceeded!");
	  return;
      }	   
  }

  // We should be done
}
