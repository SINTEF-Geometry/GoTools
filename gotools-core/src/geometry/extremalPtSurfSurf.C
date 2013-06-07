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
using std::max;
using std::min;
#include "GoTools/geometry/GeometryTools.h"
#include "GoTools/geometry/ParamCurve.h"
#include "GoTools/geometry/extremalPtSurfSurf.h"


/* @bsp ------------------------------------
I slutten av extremalPtSurfSurf(..) :
  // Test if the iteration is close to a knot
Hva skal gjoeres her?


Returner flagg =1:Extremum found. =0:Extremum NOT found
I shsing testes paa avstanden mellom punktene i de to flatene.
I shsing_ext testes paa vinkelen mellom flatenormalene i de to punktene.
Kommentar i prog:
 // Unsure about what is right here , angle between normals and difference vector ??


Test for konvergens i shsing og shsing_ext :
		  
          if ((fabs(t1[0]/tdelta[0]) <= REL_COMP_RES) &&
	      (fabs(t1[1]/tdelta[1]) <= REL_COMP_RES) &&
	      (fabs(t1[2]/tdelta[2]) <= REL_COMP_RES) &&	   
	      (fabs(t1[3]/tdelta[3]) <= REL_COMP_RES))
	      {
Endringen i parameterverdiene for flate 2: t1[2] og t1[3] er blitt satt til
null og blir aldri endret. Er det meningen?

Toleranser.
------------------------------- @bsp */


// Anonymous namespace
namespace {
  const double DNULL = (double)0.0;
  const double TOL = 1.0e-16;
  const double REL_PAR_RES =  (double)0.000000000001;     // e-12
  const double REL_COMP_RES = (double)0.000000000000001;  // e-15
}

namespace Go {

// (shsing_ext)
//===========================================================================
int extremalPtSurfSurf(ParamSurface* psurf1, ParamSurface* psurf2,
		      int    constraints[2], double constraints_par[2],
		      double limit[], double enext[],
		      double gpos[], double angle_tol)
//===========================================================================
/*
*
* PURPOSE    : To find one point in each surface where the two normals
*              are parallel to each other and to the difference vector
*              between the two points.
*              The edge info makes it possible to constrain the search
*              to an ISO-curve in one of or both surfaces.
*
* INPUT      : psurf1 - Pointer to first surface
*              psurf2 - Pointer to second surface
*              constraints[0]- Constraint 1. surface
*              constraints[1]- Constraint 2. surface
*                       = -1 No constraints
*                       =  0 Constant in 1. par direction
*                       =  1 constant in 2. par direction
*
*                           Example constrains == 0
*
*                           -----------------
*                           !     |          !
*                           !     |          !
*                           !     |          !
*                           !     |          !
*                           !     |          !
*                           !     |          !
*                           -----------------
*                                 ^
*                                 |
*                                 constraints_par
*              constraints_par[0]- Constraints parameter 1. surface
*              constraints_par[1]- Constraints parameter 2. surface
*
*              limit  - Parameter borders of both surfaces.
*                       limit[0] - limit[1] Parameter interval 1.  surface 1. direction.
*                       limit[2] - limit[3] Parameter interval 1.  surface 2. direction.
*                       limit[4] - limit[5] Parameter interval 2.  surface 1. direction.
*                       limit[6] - limit[7] Parameter interval 2.  surface 2. direction.
*              enext     - Parameter start value for iteration(4 values).
*              angle_tol - The angular tolerance for success.
*
* OUTPUT     : gpos    - Parameter values of the found singularity.(4 values)
*              kstat   - status messages
*                                = 1   : Extremum found.
*                                = 0   : Extremum NOT found.
*
*
* METHOD     :  - Start with a guess value (u,v) in domain of surface 1 (S(u,v))
*           (a) - Find domain value (r,t) of closest point (to S(u,v) in surface 2 (Q(r,t))
*               - If vf1(u,v) = <Su,Normal(Q> and vf2(u,v)= <Sv,Normal(Q> is small enough  stop
*                      (<,> means scalar prod.)
*               - Find du and dv by taylorizing vf1 and vf2.
*                 This include finding the derivatives of the closest point function (r(u,v),t(u,v))
*                 with respect to u and v. (called h(u,v) in article, see comments in nextStep)
*               - u:= u+du v:= v+dv, goto (a)
*
*
* REFERENCES : Solutions of tangential surface and curve intersections.
*              R P Markot and R L Magedson
*              Computer-Aided Design; vol. 21, no 7 sept. 1989, page 421-429
*
*
*********************************************************************
*/
{
  int kstat = 0;            /* Local status variable.                      */
  int ki;                   /* Loop control                                */
  int kp;                   /* Loop control                                */
  const int kder=2;         /* Order of derivatives to be calulated        */
  const int kdim=3;         /* Dimension of space the surface lies in      */
  int knbit;                /* Number of iterations                        */
  double tdelta[4];         /* Length of parameter intervals.              */
  double tdist;             /* The current norm of the cross product       */
                            /* between the two normals                     */
  double tprev;             /* The current norm of the cross product  ?    */
                            /* between the two normals                     */
  double td[2],t1[2],tdn[2];/* Distances between old and new parameter     */
			    /* value in the two parameter directions.     */
  std::vector<Point> sval1(7,Point(kdim)); /* Value ,first and second      */
					   /* derivative of first surface  */
  Point& snorm1=sval1[6];    /* Normal vector of firstrface                */
  std::vector<Point> sval2(7,Point(kdim)); /* Value ,first and second      */
                                           /* derivative of second surface */
  Point& snorm2=sval2[6];   /* Normal vector of second surface             */
  double snext[4];          /* Parameter values                            */
  double start[2];          /* Parameters limit of second surface, used in */
                            /* call to closest point                       */
  double end[2];            /* Parameters limit of second surface, used in */
                            /* call to closest point                       */
  double guess[2];          /* Start point for closest point iteration     */
  Point& ppoint=sval1[0];   /* Contains the current position in first      */
                            /* surface used in closest point iteration     */
  double min_tol = (double)10000.0*REL_COMP_RES;
  int max_iter=20;          /* Maximal number of iteration allowed         */
  double limit_corr[8];     /* Copy of limit with edge constraints         */
  double clo_dist;          /* Distance between point and closest point    */
  Point clo_pt(3);          /* Closest point                               */
  /* --------------------------------------------------------------------- */

  // Test input.
  DEBUG_ERROR_IF(kdim != psurf1->dimension(), "Dimension mismatch.");
  DEBUG_ERROR_IF(kdim != psurf2->dimension(), "Dimension mismatch.");

  // Fetch referance numbers from the serach intervals for the surfaces.
  tdelta[0] = limit[1] - limit[0];
  tdelta[1] = limit[3] - limit[2];
  tdelta[2] = limit[5] - limit[4];
  tdelta[3] = limit[7] - limit[6];

  // Set limit values, used in closest point iteration
  start[0] = limit[4];
  start[1] = limit[6];
  end[0]   = limit[5];
  end[1]   = limit[7];

  // Take into account edge constraints
  for (int i=0;i<8;i++)
    limit_corr[i]=limit[i];    

  for (ki=0,kp=0;ki<2;ki++,kp+=4) {
    switch (constraints[ki])
    {
    case -1:
      break;
    case 0:
      limit_corr[kp]    = limit_corr[kp+1] = constraints_par[ki];
      break;
    case 1:
      limit_corr[kp+2] = limit_corr[kp+3] = constraints_par[ki];
      break;
    default:
      THROW("Wrong constraint: " << constraints[ki]);

      break;
    }
  }

  Vector2D corner1(limit_corr[4],limit_corr[6]);
  Vector2D corner2(limit_corr[5],limit_corr[7]);
  RectDomain rect_dom2(corner1,corner2);

  // Collapsed ?
  for (ki=0;ki<4;ki++) 
    tdelta[ki] = max (tdelta[ki], min_tol);

  // Initiate output variables.
  for (ki=0;ki<4;ki++)
    gpos[ki] = enext[ki];

  for (ki=0;ki<4;ki++) {
    // Snap to the limit box
    gpos[ki] = max(gpos[ki], limit_corr[2*ki]);
    gpos[ki] = min(gpos[ki], limit_corr[2*ki+1]);
  }

  //  Evaluate 0.-2. derivatives of first surface
  psurf1->point(sval1,gpos[0],gpos[1],kder);
  psurf1->normal(snorm1,gpos[0], gpos[1]);

  //  Get closest point in second surface.
  guess[0] = gpos[2];
  guess[1] = gpos[3];
  psurf2->closestPoint(ppoint,gpos[2],gpos[3],clo_pt,clo_dist,
		       REL_PAR_RES, &rect_dom2,guess);
  //		       REL_COMP_RES, &rect_dom2,guess);    @bsp


  //  Evaluate 0.-2. derivatives of second surface
  psurf2->point(sval2,gpos[2],gpos[3],kder);
  psurf2->normal(snorm2,gpos[2], gpos[3]);

  //  Get length of normal cross product
  //s6crss(snorm1,snorm2,temp);
  //tprev = s6length(temp,kdim,&kstat);
  //  temp = SIX_norm(snorm1);
  //  temp = SIX_norm(snorm2);

  // Angle between the normal vectors
  snorm1.normalize();
  snorm2.normalize();
  tprev = snorm1.angle_smallest(snorm2);

  // Compute the Newton stepdistance vector in first surface.
  nextStep(td,sval1,sval2);

  // Adjust if we are not inside the parameter intervall.
  for (ki=0;ki<2;ki++)
    t1[ki] = td[ki];
  insideParamDomain(t1,gpos,limit_corr);


  // Iteratation loop.

  for (knbit = 0; knbit < max_iter; knbit++) {

    snext[0] = gpos[0] + t1[0];
    snext[1] = gpos[1] + t1[1];

    // Evaluate 0.-2. derivatives of first surface
    psurf1->point(sval1,snext[0],snext[1],kder);
    psurf1->normal(snorm1,snext[0], snext[1]);
    sval1[6]=snorm1;

    //   Get closest point in second surface.
    guess[0] = gpos[2];
    guess[1] = gpos[3];
    psurf2->closestPoint(ppoint,snext[2],snext[3],clo_pt,clo_dist,
		         REL_PAR_RES, &rect_dom2,guess);
    //			 REL_COMP_RES, &rect_dom2,guess);   @bsp
    
    // Since the last 2 values in snext has been changed in closestPoint:
    for (ki=2;ki<4;ki++) {
      // Snap to the limit box
      snext[ki] = max(snext[ki], limit_corr[2*ki]);
      snext[ki] = min(snext[ki], limit_corr[2*ki+1]);
    }

    //  Evaluate 0.-2. derivatives of second surface
    psurf2->point(sval2,snext[2],snext[3],kder);
    psurf2->normal(snorm2,snext[2], snext[3]);
    sval2[6]=snorm2;

 
    // Get length of normal cross product
    // s6crss(snorm1,snorm2,temp);
    // tdist = s6length(temp,kdim,&kstat);
    //temp = SIX_norm(snorm1);
    //temp = SIX_norm(snorm2);
    //tdist = s6ang(snorm1,snorm2,3);
    snorm1.normalize();
    snorm2.normalize();
    tdist = snorm1.angle_smallest(snorm2);
    
    // Compute the Newton stepdistance vector.
    nextStep(tdn,sval1,sval2);

    if (tdist <= tprev) {

      // Ordinary converging.
      
      for (ki=0;ki<4;ki++)
	gpos[ki] = snext[ki];
      
      td[0] = t1[0] = tdn[0];
      td[1] = t1[1] = tdn[1];

      
      // Adjust if we are not inside the parameter interval.
      insideParamDomain(t1,gpos,limit_corr);
      
      tprev = tdist;
      
      if ((fabs(t1[0]/tdelta[0]) <= REL_COMP_RES) &&
	  (fabs(t1[1]/tdelta[1]) <= REL_COMP_RES)) {

	gpos[0] += t1[0];
	gpos[1] += t1[1];
	// Evaluate 0.-2. derivatives of first surface
	psurf1->point(sval1,gpos[0],gpos[1],kder);
	psurf1->normal(snorm1,gpos[0], gpos[1]);
	sval1[6]=snorm1;
	
	
	// Get closest point in second surface.
	guess[0] = gpos[2];
	guess[1] = gpos[3];
	psurf2->closestPoint(ppoint,gpos[2],gpos[3],clo_pt,clo_dist,
			     REL_PAR_RES, &rect_dom2,guess);
	//			     REL_COMP_RES, &rect_dom2,guess);
	break;
      }
    }

    else {
      
      // Not converging, half step length and try again.
      
      for (ki=0;ki<2;ki++) 
	t1[ki] *= 0.5;
    }
  }

  // Iteration stopped, test if point is extremum
  // Unsure about what is right here , angle between normals and difference vector ??
  if (tprev <= angle_tol)
    kstat = 1;
  else
    kstat = 0;


  /*
  // Test if the iteration is close to a knot
  if (fabs(gpos[0] - psurf1->et1[kleftt])/tdelta[0] < min_tol)
    gpos[0] = psurf1->et1[kleftt];
  else if (fabs(gpos[0] - psurf1->et1[kleftt+1])/tdelta[0] < min_tol)
    gpos[0] = psurf1->et1[kleftt+1];

  if (fabs(gpos[1] - psurf1->et2[klefts])/tdelta[1] < min_tol)
    gpos[1] = psurf1->et2[klefts];
  else if (fabs(gpos[1] - psurf1->et2[klefts+1])/tdelta[1] < min_tol)
    gpos[1] = psurf1->et2[klefts+1];

  if (fabs(gpos[2] - psurf2->et1[kleftu])/tdelta[2] < min_tol)
    gpos[2] = psurf2->et1[kleftu];
  else if (fabs(gpos[2] - psurf2->et1[kleftu+1])/tdelta[2] < min_tol)
    gpos[2] = psurf2->et1[kleftu+1];

  if (fabs(gpos[3] - psurf2->et2[kleftv])/tdelta[3] < min_tol)
    gpos[3] = psurf2->et2[kleftv];
  else if (fabs(gpos[3] - psurf2->et2[kleftv+1])/tdelta[3] < min_tol)
    gpos[3] = psurf2->et2[kleftv+1];
  */

  // Iteration completed.
  return kstat;
}

// Anonymous namespace
namespace {

//===========================================================================
void insideParamDomain(double gd[], double coef[],double limit[])
//===========================================================================
/*
*
* PURPOSE    : To be sure that we are inside the border of the
*              parameter plane. If we are outside clipping is used
*	       to adjust the step value.
*
*
* INPUT      : coef    - Current position.
*              limit   - Parameter borders of both surfaces.
*
*
*
* INPUT/OUTPUT : gd    - Proposed delta values.
*
*
* METHOD     : Cutting the line towards the parameter box.
*
*
* REFERENCES :
*
*
*********************************************************************
*/
{
  int ki;

  for (ki=0;ki<2;ki++) {
    if (coef[ki] + gd[ki] < limit[2*ki])
      gd[ki] = limit[2*ki]    - coef[ki];
    else if (coef[ki] + gd[ki] > limit[2*ki+1])
      gd[ki] = limit[2*ki +1] - coef[ki];
  }

}

//===========================================================================
void nextStep(double cdiff[],std::vector<Point> evals,
	      std::vector<Point> evalq)
//===========================================================================
/*
*
* PURPOSE    : To calculate the increments in the first surface domain.
*
*
*
* INPUT      : evals - Value and derivatives on first surface.
*              evalq - Value and derivatives on second surface.
*
*
* OUTPUT     : cdiff1  - Parameter increments in two directions.
*
*
*
* METHOD     : See comments in main header.
*              The only thing missing in the article is the derivation of the
*              function h(u,v). Calling this function (r(u,v), t(u,v))
*              we know that
*              For all (u,v) (<,> meaning scalar product)
*              <S(u,v)-Q(r(u,v),t(u,v)),Qt(r(u,v),t(u,v))> = 0
*              <S(u,v)-Q(r(u,v),t(u,v)),Qr(r(u,v),t(u,v))> = 0
*              This means that (derivation by u and v)
*              <Su-[Qt*tu+Qr*ru],Qt> + <S-Q,Qtt*tu + Qtr*ru> = 0
*              <Su-[Qt*tu+Qr*ru],Qr> + <S-Q,Qrt*tu + Qrr*ru> = 0
*              <Sv-[Qt*tv+Qr*rv],Qt> + <S-Q,Qtt*tv + Qtr*rv> = 0
*              <Sv-[Qt*tv+Qr*rv],Qr> + <S-Q,Qrt*tv + Qrr*rv> = 0
*              Solving these four equations gives us ru,rv,tu,tv.
*
* REFERENCES :
*
*
*********************************************************************
*/
{


  int ki;                             // Loop control.                      
  static Point nq_u(3), nq_v(3);      // Derivatives of second surface normal
                                      // (with u and v !)
  static Point help1(3), help2(3);    // Help vectors
  static Point help3(3), help4(3);    // Help vectors
  double matr[4];                     // Matrix in linear equation to be solved
  static Point sq(3);                 // The difference vector S-Q
  double h_u[2];                      // The partial derivative of h() by u
  double h_v[2];                      // The partial derivative of h() by v
  double h[2];                        // Right hand side of equation system
  //  int kstat;                          // Local status
  static Point nq(3);                 // Vector for cross product
  //--------------------------------------------------------------------------

  cdiff[0] = DNULL;
  cdiff[1] = DNULL;

  /* Init, Set references to input values */

  Point& sval = evals[0];    // Reference to first surface value
  Point& qval = evalq[0];    // Reference to second surface value

// References to first surface derivatives
  Point& s_u  = evals[1];
  Point& s_v  = evals[2];
  Point& s_uu = evals[3];
  Point& s_uv = evals[4];
  Point& s_vv = evals[5];
  //  Point& ns   = evals[6];   //References to first surface normal vector

// References to second surface derivatives
  Point& q_t  = evalq[1];
  Point& q_r  = evalq[2];
  Point& q_tt = evalq[3];
  Point& q_tr = evalq[4];
  Point& q_rr = evalq[5];
  //  Point& nq   = evalq[6];

  //  SIX_cross(q_t,q_r,nq);
  nq = q_t.cross(q_r);

  // Get the difference vector S-Q
  sq = sval;
  sq -= qval;

  // Find the derivatives of the h() function by solving 2 2x2 systems
  // (same matrix).  Using Cramer's rule.

  matr[0] = q_tt*sq - q_t*q_t;    matr[1] = q_tr*sq - q_t*q_r;
  matr[2] = matr[1];              matr[3] = q_rr*sq - q_r*q_r;

  double det = matr[0]*matr[3] - matr[1]*matr[2];
  if (fabs(det) < TOL)                 //@bsp
    return;

  // Right hand side
  h[0] = -s_u*q_t;
  h[1] = -s_u*q_r;
  // Solve
  h_u[0] = (h[0]*matr[3]-matr[1]*h[1])/det;
  h_u[1] = (matr[0]*h[1]-h[0]*matr[2])/det;

 // Right hand side
  h[0] = -s_v*q_t;
  h[1] = -s_v*q_r;
  // Solve
  h_v[0] = (h[0]*matr[3]-matr[1]*h[1])/det;
  h_v[1] = (matr[0]*h[1]-h[0]*matr[2])/det;


  //     Construct matrix for finding du and dv

  help1 = q_tt*h_u[0] + q_tr*h_u[1];
  help2 = q_tr*h_u[0] + q_rr*h_u[1];
  help3 = help1.cross(q_r); 
  help4 = q_t.cross(help2);

  nq_u = help3 + help4;

  help1 = q_tt*h_v[0] + q_tr*h_v[1];
  help2 = q_tr*h_v[0] + q_rr*h_v[1];
  help3 = help1.cross(q_r); 
  help4 = q_t.cross(help2); 
  
  nq_v = help3 + help4;

  for (ki=0;ki<4;ki++) 
    matr[ki] = DNULL;

  for (ki=0;ki<3;ki++) {
    matr[0] += s_uu[ki]*nq[ki] + s_u[ki]*nq_u[ki];
    matr[1] += s_uv[ki]*nq[ki] + s_u[ki]*nq_v[ki];
    matr[2] += s_uv[ki]*nq[ki] + s_v[ki]*nq_u[ki];
    matr[3] += s_vv[ki]*nq[ki] + s_v[ki]*nq_v[ki];
  }

  // solve the linear 2x2 system using Cramer's rule.

  det = matr[0]*matr[3] - matr[1]*matr[2];
  if (fabs(det) < TOL) {                 //@bsp
    if     (matr[0] != DNULL) cdiff[0] = - s_u*nq/matr[0];
    else if(matr[1] != DNULL) cdiff[1] = - s_u*nq/matr[1];
    else if(matr[2] != DNULL) cdiff[0] = - s_v*nq/matr[2];
    else if(matr[3] != DNULL) cdiff[1] = - s_v*nq/matr[3];
    
  }
  else {
    h[0] = -s_u*nq;
    h[1] = -s_v*nq;

    cdiff[0] = (h[0]*matr[3]-matr[1]*h[1])/det;
    cdiff[1] = (matr[0]*h[1]-h[0]*matr[2])/det;
  }

}

} // Anonymous namespace

} // namespace Go  
