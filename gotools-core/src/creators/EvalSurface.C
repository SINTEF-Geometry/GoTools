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

#include "GoTools/creators/EvalSurface.h"
#include "GoTools/utils/Point.h"
#include "GoTools/geometry/Utils.h"

using namespace std;

namespace Go
{

EvalSurface::~EvalSurface()
{}

  void EvalSurface::closestPoint(const Point& pt,
				 double&        clo_u,
				 double&        clo_v, 
				 Point&         clo_pt,
				 double&        clo_dist,
				 double         epsilon,
				 int            maxiter,
				 double   *seed) const
  {
    // Iteration 
    double start_par[2], par[2], minpar[2], maxpar[2];
    double dist;
    minpar[0] = start_u();
    minpar[1] = start_v();
    maxpar[0] = end_u();
    maxpar[1] = end_v();
    if (seed)
    {
      start_par[0] = seed[0];
      start_par[1] = seed[1];
    }
    else
    {
      start_par[0] = 0.5*(minpar[0] + maxpar[0]);
      start_par[1] = 0.5*(minpar[1] + maxpar[1]);
    }

    int kstat = 0;
    s1773(pt.begin(), epsilon, minpar, maxpar, start_par, par, 
	  clo_pt, maxiter, &kstat);
    clo_u = par[0];
    clo_v = par[1];
    clo_dist = pt.dist(clo_pt);
   }
 
void 
EvalSurface::s1773(const double ppoint[],double aepsge, 
		   double estart[],double eend[],double enext[],
		   double gpos[], Point& ppoint2,
		   int maxiter, int *jstat) const
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : Newton iteration on the distance function between
*              a surface and a point, to find a closest point or an
*              intersection point.
*              If a bad choice for the guess parameters is given in, the
*              iteration may end at a local, not global closest point.
*
*
* INPUT      : ppoint   - The point in the closest point problem.
*              psurf    - The surface in the closest point problem.
*              aepsge   - Geometry resolution.
*              estart   - Surface parameters giving the start of the search
*                         area (umin, vmin).
*              eend     - Surface parameters giving the end of the search
*                         area (umax, vmax).
*              enext    - Surface guess parameters for the closest point
*                         iteration.
*
*
*
* OUTPUT     : gpos    - Resulting surface parameters from the iteration.
*              ppoint2 - Position
*              jstat   - status messages  
*                                = 2   : A minimum distanse found.
*                                = 1   : Intersection found.
*                                < 0   : error.
*
*
* METHOD     : Newton iteration in two parameter direction.
*
*
* REFERENCES :
*
*
* WRITTEN BY : Arne Laksaa, SI, May 1989
* Revised by : Johannes Kaasa, SINTEF Oslo, August 1995.
*              Introduced a local copy of enext, to avoid changes.
*
*********************************************************************
*/                       
{                        
  // int kstat = 0;            /* Local status variable.                      */
  int kder=1;               /* Order of derivatives to be calculated       */
  int kdim;                 /* Dimension of space the curves lie in        */
  int knbit;                /* Number of iterations                        */
  int kdir;                 /* Changing direction.                         */
  int kdeg;                 /* Degenaracy flag.                            */
  double tdelta[2];         /* Parameter intervals of the surface.         */
  double tdist;             /* Distance between position and origo.        */
  double td[2],t1[2],tdn[2];/* Distances between old and new parameter
			       value in the tree parameter directions.     */
  double tprev;             /* Previous difference between the curves.     */
  vector<double> sdiff;      /* Difference between the point and the surf.  */
  double snext[2];          /* Parameter values                            */
  double guess[2];          /* Local copy of enext.                        */
  double REL_COMP_RES = 1.0e-12; //1.0e-15;
  Point pts[3];
  double fac = 1.5;
  
  guess[0] = enext[0];
  guess[1] = enext[1];
  
  kdim = dim();
  
  
  /* Fetch endpoints and the intervals of parameter interval of curves.  */
  
  tdelta[0] = std::max(start_u(), end_u() - start_u());
  tdelta[1] = std::max(start_v(), end_v() - start_v());
  
  /* Allocate local used memory */
  sdiff.resize(kdim);
  
  /* Initiate variables.  */
  
  tprev = 1.0e10;
  
  /* Evaluate 0-1.st derivatives of surface */
  /* printf("\n lin: \n %#20.20g %#20.20g",
     guess[0],guess[1]); */
  
  eval(guess[0], guess[1], kder, pts);
  
  /* Compute the distanse vector and value and the new step. */
  
  s1773_s9dir(&tdist,td,td+1,&sdiff[0],ppoint,pts,
	      aepsge,kdim,&kdeg);
  
  /* Correct if we are not inside the parameter intervall. */
  
  t1[0] = td[0];
  t1[1] = td[1];
  s1773_s9corr(t1,guess[0],guess[1],estart[0],eend[0],estart[1],eend[1]);
  
  /* Iterate to find the intersection point.  */
  
  tprev = tdist;
  for (knbit = 0; knbit < maxiter; knbit++)
    {
      /* Evaluate 0-1.st derivatives of surface */
      
      snext[0] = guess[0] + t1[0];
      snext[1] = guess[1] + t1[1];
      
      eval(snext[0], snext[1], kder, &pts[0]);
      
      /* Compute the distanse vector and value and the new step. */
      
      s1773_s9dir(&tdist,tdn,tdn+1,&sdiff[0],ppoint,
	    pts,aepsge,kdim,&kdeg);
      
      /* Check if the direction of the step have change. */
      
      kdir = (Utils::inner(td, td+2, tdn) >= 0.0);     /* 0 if changed. */
      
      /* Ordinary converging. */
      
      if (tdist < tprev/(double)2 || (kdir && tdist < fac*tprev))
	{
	   guess[0] += t1[0];
	   guess[1] += t1[1];
  
	  /* printf("\n %#20.20g %#20.20g",
	     guess[0],guess[1]); */
  
	  
          td[0] = t1[0] = tdn[0];
          td[1] = t1[1] = tdn[1];
	  
	  /* Correct if we are not inside the parameter intervall. */
	  
	  s1773_s9corr(t1,guess[0],guess[1],estart[0],eend[0],estart[1],eend[1]);
          tprev = tdist;

	  if ( (fabs(t1[0]/tdelta[0]) <= REL_COMP_RES) &&
	      (fabs(t1[1]/tdelta[1]) <= REL_COMP_RES)) break;
	}
      
      /* Not converging, adjust and try again. */
      
      else if (knbit < maxiter-1)
	{
          t1[0] /= (double)2;
          t1[1] /= (double)2;
          /* knbit--;  */
	}
      if (guess[0]==guess[0]+t1[0] &&
	  guess[1]==guess[1]+t1[1]) break;
    }
  
  /* Iteration stopped, test if point founds found is within resolution */
  
  if (tdist <= aepsge)
  {
     *jstat = 1;
     /* printf("\n SUCCESS!!"); */
     
  }
  else if(kdeg)
     *jstat = 9;
  else
     *jstat = 2;
  
  if (tdist < tprev)
    {
      gpos[0] = snext[0];
      gpos[1] = snext[1];
      ppoint2 = pts[0];
    }
  else
    {
      gpos[0] = guess[0];
      gpos[1] = guess[1];
      ppoint2 = eval(gpos[0], gpos[1]);
    }

  return;

}

void
EvalSurface::s1773_s9corr(double gd[],double acoef1,double acoef2,
			    double astart1,double aend1,double astart2,double aend2) const
/*
*********************************************************************
*
*********************************************************************
*                                                                   
* PURPOSE    : To be sure that we are inside the boorder of the
*              parameter plan. If we are outside clipping is used
*	       to corrigate the step value.
*
*
* INPUT      : acoef1  - Coeffisient in the first direction.
*              acoef2  - Coeffisient in the second direction.
*              astart1 - The lower boorder in first direction.
*              aend1   - The higher boorder in first direction.
*              estart2 - The lower boorder in second direction.
*              eend2   - The higher boorder in second direction.
*
*
*
* INPUT/OUTPUT : gd    - Old and new step value.
*
*
* METHOD     : We are cutting a line inside a rectangle.
*	       In this case we always know that the startpoint of
*	       the line is inside the rectangel, and we may therfor
*	       use a simple kind of clipping.
*
*
* REFERENCES :
*
*
* WRITTEN BY : Arne Laksaa, SI, Feb 1989
*
*********************************************************************
*/                       
{
  if (acoef1 + gd[0] < astart1)  gd[0] = astart1 - acoef1;
  else if (acoef1 + gd[0] > aend1) gd[0] = aend1 - acoef1;
  
  if (acoef2 + gd[1] < astart2)  gd[1] = astart2 - acoef2;
  else if (acoef2 + gd[1] > aend2) gd[1] = aend2 - acoef2;
}

void
EvalSurface::s1773_s9dir(double *cdist,double *cdiff1,double *cdiff2,
			 double PS[],const double *eval1, Point eval2[],
			 double aepsge, int idim,int *jstat) const
/*
*********************************************************************
*
*********************************************************************
*                                                                   
* PURPOSE    : To compute the distance vector and value beetween
*	       a point and a point on a surface.
*	       And to compute a next step on both parameter direction
*	       This is equivalent to the nearest way to the
*	       parameter plan in the tangent plan from a point in the
*	       distance surface between a point and a surface.
*
*
* INPUT      : eval1 - Value in point.
*              eval2 - Value +1 and 2. der in surface.
*	       aepsge- Geometry tolerance.
*	       idim  - Dimension of space the surface lie in.
*
*
* OUTPUT     : PS   - Array to use when computing the differens vector.
*	       cdiff1  - Relative parameter value in intersection 
*                        point in first direction.
*              cdiff2  - Relative parameter value in intersection 
*                        point in second direction.
*              cdist   - The value to the point in the distance surface.
*              jstat   - 0 OK, new No degeneracy.
*                        1 Degeneracy.
*
*
* METHOD     : The method is to compute the parameter distance to the points
*	       on both tangents which is closest to the point.
*	       The difference vector beetween these points are orthogonal
*	       to both tangents. If the distance vector beetween the point and
*	       point on the surface is "diff" and the two derivativ vectors
*	       are "der1" and "der2", and the two wanted parameter distance
*	       are "dt1" and "dt2", then we get the following system of 
*	       equations:
*		 <-dist+dt1*der1+dt2*der2,der1> = 0
*		 <-dist+dt1*der1+dt2*der2,der2> = 0
*	       This is futher:
*
*		 | <der1,der1>   <der1,der2> |  | dt1 |   | <dist,der1> |
*		 |                           |  |     | = |     	|
*		 | <der1,der2>   <der2,der2> |  | dt2 |   | <dist,der2> |
*
*	       The solution of this matrix equation is the
*	       following function.
*
*
* REFERENCES :
*
*
* WRITTEN BY : Arne Laksaa, SI, Feb 1989
*
*********************************************************************
*/                       
{                        
  // int kstat=0;		          /* Local status variable.       */
  register double tdet;		  /* Determinant                  */
  register double t1,t2,t3,t4,t5; /* Variables in equation system */
  register double *S, *Su, *Sv;   /* Pointers to surf values      */
  register double ref, ang;       /* Referance value, angle       */
  register double l1, l2;         /* Vector norm                  */
  register double tcos;
  register double min_ang=10e-11; /* Min angle                    */
  register double ptol = 1.0e-12; /* Replace DEQUAL               */
  /* ____________________________________________________________ */
  
  /* Init */
  *jstat = 0;
  *cdiff1 = 0.0;
  *cdiff2 = 0.0;
  
  /* Set pointers */
  S   = eval2[0].begin();
  Su  = eval2[1].begin();
  Sv  = eval2[2].begin();

  /* Degenerate if Su=0 v Sv=0 v Su||Sv */
  l1 = sqrt(Utils::inner(Su, Su+idim, Su));
  l2 = sqrt(Utils::inner(Sv, Sv+idim, Sv));
  tcos = Utils::inner(Su, Su+idim, Sv)/(l1*l2);
  ang = acos(tcos);
  if (std::min(l1,l2) < aepsge || ang < min_ang) 
    *jstat = 1;

  /* Computing difference vector and lenght */
  for (int ki=0; ki<idim; ki++)
      PS[ki] = eval1[ki] - S[ki];
  *cdist = sqrt(Utils::inner(PS, PS+idim, PS));
  
  if (*jstat == 1)
  {
     if (l1 < aepsge)
     {
	if (l2 > aepsge)
	   /* Su = 0 */
	   *cdiff2 = Utils::inner(PS, PS+idim, Sv)/l2*l2;
     }
     else if (l2 < aepsge)
	   /* Sv = 0 */
	   *cdiff1 = Utils::inner(PS, PS+idim, Su)/(l1*l1);
     else /* Su,Sv || */
     {
	/* Best strategy? */
	*cdiff1 = Utils::inner(PS, PS+idim, Su)/(l1*l1);
      }
	
  }
  else /* *jstat == 0 */
     
  {
     
     t1 =  Utils::inner(Su, Su+idim, Su);
     t2 =  Utils::inner(Su, Su+idim, Sv);
     t3 =  Utils::inner(Sv, Sv+idim, Sv);
     t4 =  Utils::inner(PS, PS+idim, Su);
     t5 =  Utils::inner(PS, PS+idim, Sv);
     
     ref = std::max(fabs(t1),fabs(t2));
     ref = std::max(ref,fabs(t3));
     /* Computing the determinant. */
     
     tdet = t1*t3 - t2*t2;
     
     if (fabs(tdet) < ptol)
     {
	*jstat = 1;
     }
     else 
     {
	/* Using Cramer's rule to find the solution of the system. */
	
	*cdiff1 =  (t4*t3-t5*t2)/tdet;
	*cdiff2 =  (t1*t5-t2*t4)/tdet;
     }
  }
}

void EvalSurface::write(std::ostream& out) const
{
  return;
}

} // namespace Go

