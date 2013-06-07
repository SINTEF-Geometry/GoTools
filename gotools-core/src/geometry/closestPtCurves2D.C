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
#include "GoTools/geometry/GeometryTools.h"
#include "GoTools/geometry/ClosestPoint.h"
#include "GoTools/utils/Values.h"          // MAXDOUBLE

//***************************************************************************
//
// Implementation file of the free function ClosestPoint::closestPtCurves2D in namespace
// ClosestPoint, declared in ClosestPoint.h
//
//***************************************************************************

// Anonymous namespace
namespace {
  const double TOL = 1.0e-16;
  const double REL_COMP_RES = 0.000000000000001;
}


// Anonymous namespace for helper functions declaration.
namespace
{
using namespace Go;
void insideParamDomain2D(Point& gd, const Point& acoef, double astart1, 
			 double aend1, double astart2,double aend2, int& corr);

void nextStep2D(double& dist, Point& diff, Point& delta, double& det,
		int& kstat, std::vector<Point>& eval1, 
		std::vector<Point>& eval2, int method);

void makeEqSys1(double A[], double b[], const std::vector<Point>& eval1,
		const std::vector<Point>&eval2, const Point& diff);

void makeEqSys2(double A[], double b[], const std::vector<Point>& eval1,
		const std::vector<Point>&eval2, const Point& diff);

void makeEqSys3(double A[], double b[], const std::vector<Point>& eval1,
		const std::vector<Point>&eval2, const Point& diff);

int localPretop(double dist,const Point& diff,
		const std::vector<Point>& eval1,
		const std::vector<Point>& eval2);

void singular(ParamCurve* cv1, ParamCurve* cv2,Point& par_val,double& dist,
	      int quick,double aepsge,double delta,const Point& diff,
	      const std::vector<Point>&eval1,const std::vector<Point>&eval2,
	      double astart1,double astart2,double aend1,double aend2);

void secant2D(ParamCurve *pcurve1,ParamCurve *pcurve2,Point& par_val,
	      double& dist,int& jstat, double delta,double aepsge,
	      double astart1,double astart2,double aend1,double aend2);

void setReturnValues(const Point& par_val,ParamCurve* cv1,ParamCurve* cv2,
		     double aepsge,  int& jstat,double& par1,double& par2,
		     double& dist,Point& ptc1,Point& ptc2);

void newPointEval(ParamCurve *pcurve1,double par1,ParamCurve *pcurve2,
		  double astart2,double aend2,double aepsge,
		  double& par2,double& y,double& dist,int& jstat);

} // End of anonymous namespace for helper functions declaration.


namespace Go {

// (s1770_2D)
//===========================================================================
void ClosestPoint::closestPtCurves2D(ParamCurve* cv1, ParamCurve* cv2, double aepsge,
		     double tmin1, double tmax1, double tmin2, double tmax2,
		     double seed1, double seed2, int method, bool quick, 
		     double& par1, double& par2, double& dist, 
		     Point& ptc1, Point& ptc2, int& istat)
//===========================================================================
/*
* PURPOSE    : Newton iteration on the distance function between
*              two curves in 2D to find a closest point or an 
*              intersection point.
*
*
* INPUT      : cv1     - Pointer to the first curve in the intersection.
*              cv2     - Pointer to the second curve in the intersection.
*              aepsge  - Geometry resolution.
*              tmin1   - Start parameter value of the first curve.
*              tmax1   - End parameter value of the first curve.
*              tmin2   - Start parameter value of the second curve.
*              tmax2   - End parameter value of the second curve.
*              seed1   - Start parameter value of the iteration on
*                        the first curve.
*              seed2   - Start parameter value of the iteration on
*                        the second curve.
*
*
*
* OUTPUT     : par1    - Parameter value of of first curve in intersection
*                        point.
*              par2    - Parameter value of of second curve in intersection
*                        point.
*              dist    - Distance in space between the points.
*              ptc1    - Point on curve number one.
*              ptc2    - Point on curve number two.
*              istat   - status messages.
*                                = 1   : Intersection found.
*                                = 2   : A minimum distance found.
*                                = 3   : Nothing found.
*
*
* METHOD     : Newton iteration in two parameter directions.
*
*********************************************************************
*/
{
  const double DZERO = (double)0.0;
  const int dim = cv1->dimension();
  DEBUG_ERROR_IF(dim != cv2->dimension(), "Dimension mismatch.");
  DEBUG_ERROR_IF(dim != 2, "Only implemented for 2D");


  // Estimated start values  
  Point par_val(seed1, seed2);

  // Parameter end points
  double aend1 = cv1->endparam();
  double aend2 = cv2->endparam();

  double delta[2];          // Parameter interval of the curves. 
  delta[0] = aend1 - cv1->startparam();
  delta[1] = aend2 - cv2->startparam();

  //  double tprev = MAXDOUBLE;
  bool from_right;
  int alt_method;
  if (method==2)
    alt_method=1;
  else
    alt_method=method;

  int order;
  if (method==2)
    order=2;
  else
    order=1;

  // Evaluate 0-2.st derivatives of both curves.
  std::vector<Point> sval1(3,Point(0.0,0.0)), sval2(3,Point(0.0,0.0));
  from_right = (par_val[0] == aend1) ? false : true;
  cv1->point(sval1, par_val[0], 2, from_right);
  from_right = (par_val[1] == aend2) ? false : true;
  cv2->point(sval2, par_val[1], 2, from_right);

  Point d(2), c_d(2), nc_d(2);  // Distances between old and new parameter-
                                // value in the two parameter directions.
  Point normal_cv1(dim), normal_cv2(dim);  // Curve normals 
  Point diff(dim), prev_diff(dim);  //Difference vector between the curves.

 // Compute the distance vector, value and the new step.
  double det;
  int kstat;
  nextStep2D(dist, diff, c_d, det, kstat,  sval1, sval2, method);

  if (kstat == 1) {
    // Singular matrix. Try second order.
    if (order!=2) {
      order=2;
      method=2;
      nextStep2D(dist, diff, c_d, det, kstat,  sval1, sval2, method);
    }
    if (kstat == 1) {
      // Singular matrix.
      if (dist > aepsge)  
	// Singularity (closest point) found.
	singular(cv1,cv2,par_val,dist,quick,aepsge,c_d[0],diff,sval1,sval2,
		 tmin1,tmin2,tmax1,tmax2);
      setReturnValues(par_val,cv1,cv2,aepsge,istat,par1,par2,dist,ptc1,ptc2);
      return;
    }
  }

  if (method == 4) {
    secant2D(cv1,cv2,par_val,dist,kstat,c_d[0],aepsge,
	     tmin1,tmin2,tmax1,tmax2);
    setReturnValues(par_val,cv1,cv2,aepsge,istat,par1,par2,dist,ptc1,ptc2);
    return;
  }

  
  double prev_dist;         // Previous difference between the curves.
  int max_it = 20;          // Maximum number of iterations.
  int sing = 0;		    // Mark that singularity has ocured.
  int p_dir;                // Changing direction in par-space.
  int keep_order = 0;
  int corr = 0;
  int div2 = 0;
  int g_up,ng_up,g_dir;     // Changing direction in geometric space.	

  if(quick)
    max_it=10;

  // Adjust the step size if we are not inside the parameter interval
  d = c_d;
  normal_cv1[0] = -sval1[1][1];  normal_cv1[1] = sval1[1][0];  // Only 2D
  normal_cv2[0] = -sval2[1][1];  normal_cv2[1] = sval2[1][0];
  g_up =  (diff*normal_cv2 >= 0.0) ? 1 : -1;
  g_up += (diff*normal_cv1 >= 0.0) ? 10 : -10;
  insideParamDomain2D(d,par_val,tmin1,tmax1,tmin2,tmax2,corr);

  prev_dist = dist;
  prev_diff = diff;


  // Start of iteration to find the intersection point
  int knbit;
  for (knbit = 0; knbit < max_it; knbit++) {
    par_val+=d;

    // Test if the current method is OK, or should we change method?
    while (1) {
      // Evaluate derivatives of curve
      from_right = (par_val[0] == aend1) ? false : true;
      cv1->point(sval1, par_val[0], order, from_right);
      from_right = (par_val[1] == aend2) ? false: true;
      cv2->point(sval2, par_val[1], order, from_right);
      
      // Compute the distanse vector and value and the new step.
      nextStep2D(dist, diff, nc_d, det, kstat,  sval1, sval2, method);
      if (kstat == 1) {
	// Singular matrix.
	sing++;
	if (order==2) {
	  if (dist>aepsge)
	    // Singularity (closest point) found.
	    singular(cv1,cv2,par_val,dist,quick,aepsge,nc_d[0],diff,
		     sval1,sval2,tmin1,tmin2,tmax1,tmax2);
	  setReturnValues(par_val,cv1,cv2,aepsge,istat,par1,par2,dist,ptc1,ptc2);
	  return;
	} 
	else {
	  order=2;  // Use more terms in the series expansion
	  method=2;
	}
      }
      else {
	// Normal to curve 1
	normal_cv1[0] = -sval1[1][1];  normal_cv1[1] = sval1[1][0];
	// Normal to curve 2
	normal_cv2[0] = -sval2[1][1];  normal_cv2[1] = sval2[1][0];

	// Have normal and difference vectors the same direction ?
	ng_up = (diff*normal_cv2 >= 0.0) ? 1 : -1;
	ng_up += (diff*normal_cv1 >= 0.0) ? 10 : -10;

	// g_dir=1 if we have not changed position to the other side
	// of the intersection point. 0 if changed.
	g_dir = (ng_up+g_up != 0);

	// p_dir=1 if the steps in the parameter interval continues
	// in the same direction. 0 if direction has changed.
	p_dir = (c_d[0]*nc_d[0] >= DZERO &&
		 c_d[1]*nc_d[1] >= DZERO);
	
	if (order!=2 && g_dir && (!p_dir || dist > 0.4*prev_dist)
	    && !keep_order) {           // Few terms in the series expansion.
	  if (!quick && div2) div2 = 0; // Not good enough convergence. 
	  order=2;                      // Change method.
	  method=2;
	}
	else if (order==2 && !g_dir) { // Max terms in the series expansion.
	  if (sing) {        // Found closest point?
	    if (dist>aepsge)
	      singular(cv1,cv2,par_val,dist,quick,aepsge,nc_d[0],diff,
		       sval1,sval2,tmin1,tmin2,tmax1,tmax2);

	    setReturnValues(par_val,cv1,cv2,aepsge,istat,par1,par2,dist,
			    ptc1,ptc2);
	    return;
	  }
	  if (div2) div2 = 0; // Closest point not found. Change method
	  order=1;             
	  method=alt_method;
	}
	else {                   // Decreasing distance. 
	  keep_order = 0;        // Continue with current method.
	  if (sing) sing = 0;    
	  break;
	}
      }
    } // end while(1)


    // We have decided method and computed a new step.
    // Test if we are inside the intervals, and if the iteration is
    // converging fast enough.
    if (corr)
      if (!(p_dir && g_dir)) corr = 0;
    
    if (dist < prev_dist || p_dir) {
      // Adjust if we are not inside the parameter interval.

	g_up = ng_up;
	d = c_d = nc_d;
	insideParamDomain2D(d,par_val,tmin1,tmax1,tmin2,tmax2,corr);

	prev_dist = dist;
	prev_diff = diff;

	// if (corr > 3) break;

	if (corr > 2 ||
	    ((fabs(d[0]/std::max(par_val[0],delta[0])) <= REL_COMP_RES) &&
	     (fabs(d[1]/std::max(par_val[1],delta[1])) <= REL_COMP_RES)))
	  break;
	if (div2) div2 = 0;

	if (corr > 1 && order==2) {
	  keep_order = 1;
	  order=1;
	  method=alt_method;
	}
     }
     else if (corr > 2 ||
	      ((fabs(d[0]/std::max(par_val[0],delta[0])) <= REL_COMP_RES) &&
	       (fabs(d[1]/std::max(par_val[1],delta[1])) <= REL_COMP_RES)))
       break;
     else {
       // Not converging, adjust and try again.
       if (dist > prev_dist && div2 > 5)
	 break;
       if (quick && dist > prev_dist && div2 > 3)
	 break;
       div2++;                   // Try half of the step length
       par_val -= d;
       d[0] *= 0.5;  d[1] *=0.5;
     }
  }


  // Iteration stopped, test if point found is within resolution

  if (fabs(det)<0.1) {
    if (order < 2) {      // Make sure we have enough derivatives for
      order=2;            // further computations
      method=2;
      
      from_right = (par_val[0] == aend1) ? false : true;
      cv1->point(sval1, par_val[0], order, from_right);
      from_right = (par_val[1] == aend2) ? false : true;
      cv2->point(sval2, par_val[1], order, from_right);
    }
    singular(cv1,cv2,par_val,dist,quick,aepsge,d[0],diff,sval1,sval2,
	     tmin1,tmin2,tmax1,tmax2);
  }

  setReturnValues(par_val,cv1,cv2,aepsge, istat,par1,par2,dist,ptc1,ptc2);

  // Iteration completed.

}
}  // End of namespace Go



// Anonymous namespace for helper functions definition.
namespace {
using namespace Go;

// (s1770_2D_s9corr)
//===========================================================================  
void insideParamDomain2D(Point& gd, const Point& acoef, double astart1, 
			 double aend1, double astart2,double aend2, int& corr)
//===========================================================================
/*
* PURPOSE    : To be sure we are inside the intervals [astart1,aend1] and
*              [astart2,aend2]. If we are outside clipping is used to adjust
*              the step value.
*              Ported from the sisl function s1770_2D_s9corr.
*
*
* INPUT      : acoef[] - Parameter values.
*              astart1 - The lower border in first direction.
*              aend1   - The higher border in first direction.
*              astart2 - The lower border in second direction.
*              aend2   - The higher border in second direction.
*
*
* INPUT/OUTPUT : gd   - Old and new step value.
*
* OUTPUT       : corr - If the step value has been changed, corr is increased
*                       by one. If not, corr is set to zero.
*
*********************************************************************
*/
{
  int lcorr = 0;
  if (acoef[0] + gd[0] < astart1) {
    gd[0] = astart1 - acoef[0];
    lcorr=1;
  }
  else if (acoef[0] + gd[0] > aend1) {
    gd[0] = aend1 - acoef[0];
    lcorr=1;
  }


  if (acoef[1] + gd[1] < astart2) {
    gd[1] = astart2 - acoef[1];
    lcorr=1;
  }
  else if (acoef[1] + gd[1] > aend2) {
    gd[1] = aend2 - acoef[1];
    lcorr=1;
  }

  if (lcorr)
    corr++;
  else
    corr = 0;
}


// (s1770_2D_s9dir)
//===========================================================================
void nextStep2D(double& dist,Point& diff,Point& delta,double& det,int& kstat,
		std::vector<Point>& eval1,std::vector<Point>& eval2, 
		int method)
//===========================================================================
/*
* PURPOSE    : Make and solve an equation system to compute the distance
*              vector and value beetween a point on the first curve and a 
*	       point on the second curve. And to compute a next step on both
*              curves.
*	       This is equivalent to the nearest way to the
*	       parameter plane in the tangent plane from a point in the
*	       distance surface between two curves.
*
*
* INPUT      : eval1  - Value and derivatives of the first curve.
*              eval2  - Value and derivatives of the second curve.
*              method -  
*                      = 1 : The method in s1770_2D, starting with order one.
*                      = 2 : The method in s1770_2D, starting with order two.
*                      = 3 : The method from s1770 (order one).
*                      = 4 : The secant method from s1770_2D (order one). 
*                            For test purpose.
*
* OUTPUT     : dist    - Length of the distance vector between the points.
*	       diff    - Distance vector between the points.
*	       delta   - Relative step parameter values towards intersection on
*                        the curve 1 delta[0] and the curve 2 delta[1].
*              det     - Value of determinant of the equation system
*              kstat   - Status messages
*                      = 0 : Solution found.
*                      = 1 : Singular system.
*
*********************************************************************
*/
{

  double A[4];   // Equation system matrix
  double b[2];   // Equation system right hand side

  diff = eval1[0]-eval2[0]; // Difference vector between the points
  dist = diff.length();

  // Make equation system. Fills matrix A and right hand side b
  switch (method) {
    case 1:
      { 
	makeEqSys1(A, b, eval1, eval2, diff);
	break;
      }
    case 2:
      {
	makeEqSys2(A, b, eval1, eval2, diff);
	break;
      }
    case 3:
      {
	makeEqSys3(A, b, eval1, eval2, diff);
	break;
      }
    case 4:
      {
	makeEqSys1(A, b, eval1, eval2, diff);
	break;
      }
    default:
      {
	  THROW("Wrong method: " << method);
	break;
      }
      
  }

  
  // Solve the 2x2 equation system Ax=b
  det = A[0]*A[3] - A[2]*A[1];   // Determinant
  if (fabs(det) < TOL) {
    kstat = 1;
    delta[0] = 0.0;
    delta[1] = 0.0;
  }
  else {   // Using Cramer's rule to find the solution of the system
    kstat = 0;
    delta[0] = (b[0]*A[3] - b[1]*A[1])/det;
    delta[1] = (A[0]*b[1] - A[2]*b[0])/det;
  }

}

//===========================================================================
void makeEqSys1(double A[], double b[], const std::vector<Point>& eval1,
		const std::vector<Point>&eval2, const Point& diff)
//===========================================================================
/*
* PURPOSE    : Fill the 2x2 matrix A and right hand side b of the equation 
*              system Ax=b. This is an equation for computing the next step
*              of the parameters s and t for the 2D curves f and g in an
*              iteration for computing their intersection point or closest
*              point. The equation is the one in the sisl function 
*              s1770_2D_s9dir with order equal to one.
*
*
*                           | ds |   |        |
*         | -g'(s)  f'(t) | |    | = | d(s,t) |
*                           | dt |   |        |  
*
*
* INPUT   : eval1 Curve g evaluated at s
*           eval2 Curve f evaluated at t
*           diff  Distance vector d(s,t) between the space points g(s) and f(t)
*
* OUTPUT  : A Matrix
*           b Right hand side
*
*********************************************************************
*/
{
  A[0]=-eval1[1][0];   A[1]=eval2[1][0]; 
  A[2]=-eval1[1][1];   A[3]=eval2[1][1];

  b[0]=diff[0];
  b[1]=diff[1];
}

//===========================================================================
void makeEqSys2(double A[], double b[], const std::vector<Point>& eval1,
		const std::vector<Point>&eval2, const Point& diff)
//===========================================================================
/*
* PURPOSE    : Fill the 2x2 matrix A and right hand side b of the equation 
*              system Ax=b. This is an equation for computing the next step
*              of the parameters s and t for the 2D curves f and g in an
*              iteration for computing their intersection point or closest
*              point. The equation is the one in the sisl function 
*              s1770_2D_s9dir with order equal to two.
*
* 
*  | -<g'(s),g'(s)>-<d(s,t),g"(s)  <g'(s),f'(t)> | | ds |   | <d(s,t),g'(s)> |
*  |                                             | |    | = |                |
*  | -<g'(s),f'(t)>  <f'(s),f'(s)>-<d(s,t),f"(s)>| | dt |   | <d(s,t),f'(t)> |
* 
*
* INPUT    : eval1 Curve g evaluated at s
*            eval2 Curve f evaluated at t
*            diff Distance vector d(s,t) between the space points g(s) and f(t)
*
* OUTPUT     : A Matrix
*              b Right hand side
*
*********************************************************************
*/
{
  A[0]=-eval1[1]*eval1[1]-diff*eval1[2];  A[1]=eval1[1]*eval2[1]; 
  A[2]=-A[1];                             A[3]=eval2[1]*eval2[1]-diff*eval2[2];

  b[0]=diff*eval1[1]; 
  b[1]=diff*eval2[1]; 
}

//===========================================================================
void makeEqSys3(double A[], double b[], const std::vector<Point>& eval1,
		const std::vector<Point>&eval2, const Point& diff)
//===========================================================================
/*
* PURPOSE  : Fill the 2x2 matrix A and right hand side b of the equation 
*            system Ax=b. This is an equation for computing the next step
*            of the parameters s and t for the 2D curves f and g in an
*            iteration for computing their intersection point or closest
*            point. The equation is the one in the sisl function s1770_s9dir.
*
*       | -<g'(s),f'(t)>  <f'(t),f'(t)> | | ds |   | <d(s,t),f'(t)> |
*       |                               | |    | = |                |
*       | -<g'(s),g'(s)>  <g'(s),f'(t)> | | dt |   | <d(s,t),g'(s)> |
*
*
* INPUT    : eval1 Curve g evaluated at s
*            eval2 Curve f evaluated at t
*            diff Distance vector d(s,t) between the space points g(s) and f(t)
*
* OUTPUT   : A Matrix
*            b Right hand side
*
*********************************************************************
*/
{
  A[0]=-eval1[1]*eval2[1];   A[1]=eval2[1]*eval2[1]; 
  A[2]=-eval1[1]*eval1[1];   A[3]=-A[0];

  b[0]=diff*eval2[1]; 
  b[1]=diff*eval1[1];


 
}

// (s1770_2D_s6local_pretop)
//===========================================================================
int localPretop(double dist,const Point& diff,
		const std::vector<Point>& eval1,
		const std::vector<Point>& eval2)
//===========================================================================
/*
*   PURPOSE : To find if we have a minimum or a maximum or a point of
*             inflection situation. This function assumes that it is a
*	      singular situation. 
*             Ported from the sisl function s1770_2D_s6local_pretop.
*
*
*
*   INPUT   : dist     - The length of the distance vector between the points.
*	      diff     - The distance vector between the points.
*	      eval1    - Value and derivatives on the first curve.
*             eval2    - Value and derivatives on the second curve.
*
*
*   OUTPUT  : return value - 	= 0: Maximum position or inflection point.
*				= 1: Minimum position. 
*
*
*   METHOD  : Testing size and direction of the second derivatives
*             of the curve.
*
************************************************************************
*/
{
  DEBUG_ERROR_IF(diff.dimension() != 2, "Only implemented for 2D");

  int return_val;	 // For return value.
  double l_1,l_2;
  double v_1,v_2;
  double k_1,k_2;
  double r_1,r_2;

  //  const Point& c1=eval1[0];
  const Point& c1_t=eval1[1];
  const Point& c1_tt=eval1[2]; 
  //  const Point& c2=eval2[0];
  const Point& c2_t=eval2[1];
  const Point& c2_tt=eval2[2]; 

  l_1 = c1_tt*diff;
  l_2 = c2_tt*diff;

  if (( l_1 < 0.0 && l_2 > 0.0) || (l_1 > 0.0 && l_2 < 0.0))
    return 1;   // Minimum position


  v_1 = c1_t*c1_t;
  v_1 = v_1*sqrt(v_1);
  k_1 = fabs(c1_t[0]*c1_tt[1] - c1_tt[0]*c1_t[1]);
  if (k_1 < REL_COMP_RES)	r_1 = 0.0;
  else				r_1 = v_1/k_1;

  v_2 = c2_t*c2_t;
  v_2 = v_2*sqrt(v_2);
  k_2 = fabs(c2_t[0]*c2_tt[1] - c2_tt[0]*c2_t[1]);
  if (k_2 < REL_COMP_RES)	r_2 = 0.0;
  else				r_2 = v_2/k_2;


  if (( l_1 < 0.0 || l_2 < 0.0) && (r_1 > r_2 + dist))
    return_val = 1;
  else if (( l_1 > 0.0 || l_2 > 0.0) && (r_2 > r_1 + dist))
    return_val = 1;
  else
    return_val = 0;  // Maximum position or inflection point.

  return return_val;
}

//===========================================================================
void singular(ParamCurve* cv1,ParamCurve* cv2,Point& par_val,double& dist,
	      int quick,double aepsge,double delta,const Point& diff,
	      const std::vector<Point>&eval1,const std::vector<Point>&eval2,
	      double astart1,double astart2,double aend1,double aend2)
//===========================================================================
/*
*   PURPOSE : Singularity (closest point) found. If we have a point of 
*             inflection or a maximum position situation call the
*             secant method.
*
*
*   INPUT   : cv1     - Curve number one.
*	      cv2     - Curve number two.
*             quick   - Reduce requirement on exactness.
*             aepsge  - Geometry resolution.
*             delta   - Parameter distance on first curve beetveen start values
*             diff    - Distance vector between the points.
*	      eval1   - Value and derivatives of the first curve.
*             eval2   - Value and derivatives of the second curve.
*
*
*   OUTPUT  : par_val - Parameter values for the two curves.
*	      dist    - Distance between the points.
*
*
*   METHOD  : Uses function "localPretop" to test if we have a point
*             of inflection or a maximum position situation.
*
************************************************************************
*/
{
  if (!quick && dist > aepsge)
  {
    int ki=localPretop(dist,diff,eval1,eval2);
    if (ki == 0) {
      int kstat;
      secant2D(cv1,cv2,par_val,dist,kstat, delta,aepsge,
	       astart1,astart2,aend1,aend2);
    }
  }
}

// (s1770_2D_s6sekant1)
//===========================================================================
void secant2D(ParamCurve *pcurve1,ParamCurve *pcurve2,Point& par_val,
	      double& dist,int& jstat, double delta,double aepsge,
	      double astart1,double astart2,double aend1,double aend2)
//===========================================================================
/*
* PURPOSE  : Secant method iteration on the distance function between
*            two curves to find a closest point or an intersection point.
*            Ported from the sisl function s1770_2D_s6sekant1.
*
*
* INPUT    : pcurve1   - Pointer to the first curve in the intersection.
*            pcurve2   - Pointer to the second curve in the intersection.
*            delta     - Parameter distance on first curve betveen start values
*            aepsge    - Geometry resolution.
*
*
* IN/OUT   : par_val   - Parameter value of the curves in intersection point.
*
*
* OUTPUT   : dist      - Distance in space between the points.
*          : jstat     - status messages
*                                = 0   : OK
*                                = 2   : Can't find a reasonable start point
*
*
* METHOD     : Secant method in two parameter directions.
*
*********************************************************************
*/
{
  int ki;		    /* Counter.					   */
  int kstat = 0;            /* Local status variable.                      */
  int knbit;                /* Number of iterations                        */
  Point cu_val(2); 	    /* Parameter values on curve.		   */
  double new_cu_val;	    /* New parameter value on curve.		   */
  double y[2],new_y,delta_y;/* Signed distance.				   */
  //  int cu1_left = 0;	    // Keep left knot information for evaluator.
  //  int cu2_left = 0;	    // Keep left knot information for evaluator.
  int shift = 0;	    /* Mark that the direction has been changed.  */

  jstat = 0;

  // Test input.
  const int dim = pcurve1->dimension();
  DEBUG_ERROR_IF(dim != pcurve2->dimension(), "Dimension mismatch.");
  DEBUG_ERROR_IF(dim != 2, "Only implemented for 2D");

  std::vector<Point> c1(2,Point(0.0,0.0));  // Derivatives of the first curve
  std::vector<Point> c2(2,Point(0.0,0.0));  // Derivatives of the second curve
  Point norm(dim);       // Normal vector
  Point diff(dim);       // Difference vector between point/curve.
  Point pt(dim);         // Point for use in closest point point/curve.
  Point clo_pt(dim);     // Point for use in closest point point/curve.
  //  double pt_dist;        // Distance between points;

  if (delta == 0.0) delta =1e-15;

  if (par_val[0] < astart1) 
    par_val[0] = astart1;
  else if (par_val[0] > aend1) 
    par_val[0] = aend1;

  if ((par_val[0] == astart1 && delta < 0.0) ||
      (par_val[0] == aend1   && delta > 0.0)) {
    delta = -delta;       // Change direction of search along curve 1
    shift++;
  }


  // Make a "segmentation" of the parameter interval on curve 1 to get
  // a set of test values for finding start point number two.

  if (fabs(delta) < (aend1 -astart1)/100.0) {
    if (delta < 0.0)
	delta = (astart1 - aend1)/100.0;
     else
	delta = (aend1 - astart1)/100.0;
  }
  else if (fabs(delta) > (aend1 -astart1)/10.0)
  {
     if (delta < 0.0)
	delta = (astart1 - aend1)/10.0;
     else
	delta = (aend1 - astart1)/10.0;
  }



  cu_val[0] = par_val[0];
  newPointEval(pcurve1,cu_val[0],pcurve2,astart2,aend2,aepsge,
		 par_val[1],y[0],dist,kstat);
  if (kstat < 0)
    return;

  // Set parameter value to sample point  
  cu_val[1] = cu_val[0] + delta;
  if (cu_val[1] < astart1)
    cu_val[1] = astart1;
  else if (cu_val[1] > aend1)
    cu_val[1] = aend1;


  // Search for the second start point on the other side of the
  // intersection point.
  for (ki=0; ki<20; ki++) {
    newPointEval(pcurve1,cu_val[1],pcurve2,astart2,aend2,aepsge,
		   par_val[1],y[1],dist,kstat);
    if (kstat < 0) {
      par_val[0]=cu_val[1];
      return;
    }
 
    new_y = y[1]/y[0];
    if (new_y > 1.0000000000001) {
      // We are closer to a singularity. (Closest point, not intersection.)
      if (shift) {
	// We have already changed the direction. Give up.
 	par_val[0] = cu_val[1];
	return;
      }
      // Search in the other direction
      delta = -delta;
      cu_val[1] = cu_val[0] + delta;
      if (cu_val[1] < astart1) cu_val[1] = astart1;
      else if (cu_val[1] > aend1) cu_val[1] = aend1;
      shift++;
    }

    // Start points on different side of the intersection point or far from
    // a singularity?
    else if (y[0]*y[1] <= 0.0 || fabs(new_y) < 0.6)
      break;
    else {
      // try another point as start point number two
      if (cu_val[1]+delta <= aend1 && cu_val[1]+delta >= astart1)
	cu_val[1] += delta;
      else if (cu_val[1] < aend1)
  	cu_val[1] = aend1;
      else if (cu_val[1] > astart1)
	cu_val[1] = astart1;
      else {
	par_val[0] = cu_val[1];
	return;
      }
    }
   }

  // Can't find a reasonable start point
  if (ki == 20) {
    jstat = 2;
    return;
  }

  new_cu_val=par_val[0];

  // Here starts the real search for the intersection point

  for (knbit=0; knbit < 25; knbit++) {
    delta_y = y[0]-y[1];
    if (fabs(delta_y) < REL_COMP_RES)  
      break;   // Start points are too close

    // Are the start points on different side of the intersection point and
    // is the distance from one start point to the intersection point much
    // smaller than the other?
    if (y[0]*y[1] < 0.0 &&
	(fabs(y[0]) < 6*fabs(y[1]) || fabs(y[1]) < 6*fabs(y[0])))
      new_cu_val = 0.5*(cu_val[1]+cu_val[0]);  // Yes, try the midpoint
    else
      new_cu_val = cu_val[1] + y[1]*(cu_val[1]-cu_val[0])/delta_y;  // secant

    // Make sure we are inside the interval
    if (new_cu_val >= aend1) {
      new_cu_val = aend1;
      if (cu_val[0] == aend1 || cu_val[1] == aend1) {
	par_val[0] = new_cu_val;
	return;
      }
    }
    else if (new_cu_val <= astart1) {
      new_cu_val = astart1;
      if (cu_val[0] == astart1 || cu_val[1] == astart1) {
	par_val[0] = new_cu_val;
	return;
      }
    }


    // Evaluation and closest point calls in the new test point
    newPointEval(pcurve1,new_cu_val,pcurve2,astart2,aend2,aepsge,
		   par_val[1],new_y,dist,kstat);
    if (kstat < 0) {
      par_val[0]=new_cu_val;
      return;
    }

    // y[] tells if the difference vector and the curvee's normal vector
    // points in the same direction. If y changes sign, we move to the
    // other side of the intersection point.

    // How the points are distributed around the intersection point
    // decides which start point tis replaced by the new point.

    if ((y[0] < 0.0 && y[1] > 0.0) ||
	(y[0] > 0.0 && y[1] < 0.0)) {
 
      if ((new_y > 0.0 && y[0] > 0.0) ||
	  (new_y < 0.0 && y[0] < 0.0)) {
	cu_val[0] = new_cu_val;
	y[0] = new_y;
      }
      else {
	cu_val[1] = new_cu_val;
	y[1] = new_y;
      }
    }
    else {
      if ( y[0] < 0.0 && new_y > 0.0) {
	if (y[0] < y[1]) {
	  cu_val[0] = new_cu_val;
	  y[0] = new_y;
	}
	else {
	  cu_val[1] = new_cu_val;
	  y[1] = new_y;
	}
      }
      else if ( y[0] > 0.0 && new_y < 0.0) {
	if (y[0] > y[1]) {
	  cu_val[0] = new_cu_val;
	  y[0] = new_y;
	}
	else {
	  cu_val[1] = new_cu_val;
	  y[1] = new_y;
	}
      }
      else if (y[0] > 0.0) {
	if (y[0] > y[1]) {
	  if (new_y >=  y[0]) break;
	  cu_val[0] = new_cu_val;
	  y[0] = new_y;
	}
	else {
	  if (new_y >=  y[1]) break;
	  cu_val[1] = new_cu_val;
	  y[1] = new_y;
	}
	
      }
      else if (y[0] < 0.0) {
	if (y[0] < y[1]) {
	  if (new_y <=  y[0]) break;
	  cu_val[0] = new_cu_val;
	  y[0] = new_y;
	}
	else {
	  if (new_y <=  y[1]) break;
	  cu_val[1] = new_cu_val;
	  y[1] = new_y;
	}
      }
    }
  }
  par_val[0] = new_cu_val;

}


//===========================================================================
void setReturnValues(const Point& par_val,ParamCurve* cv1,ParamCurve* cv2,
		     double aepsge,  int& istat,double& par1,double& par2,
		     double& dist,Point& ptc1,Point& ptc2)
//===========================================================================
/*
* PURPOSE    :  To set the values to be returned from the function 
*               ClosestPoint::closestPtCurves2D.
*
*
* INPUT      : par_val  - Parameter values for the two curves.
*              pcurve1  - Curve number one.
*              pcurve2  - Curve number two.
*              aepsge   - Geometry resolution.
*
*
* OUTPUT     : istat    - Status variable.
*                       	= 1 : Intersection.
*                       	= 2 : Minimum value.
*                       	= 3 : Nothing found.
*              par1     - Parameter value of the first curve's point.
*              par2     - Parameter value of the second curve's point.
*              dist     - Distance in space between the points.
*              ptc1     - Point on curve number one.
*              ptc2     - Point on curve number two.
*
************************************************************************
*/
{
  const double ANGULAR_TOLERANCE = 0.01;
  const double PIHALF = 1.57079632679489661923;

  par1 = par_val[0];
  par2 = par_val[1];

  // Compute the points and the distance between them.
  std::vector<Point> eval1(2), eval2(2);
  cv1->point(eval1,par1,1); 
  cv2->point(eval2,par2,1);

  ptc1=eval1[0];
  ptc2=eval2[0];
  dist=ptc1.dist(ptc2);

  if (dist <= aepsge)
    istat=1;        // Intersection
  else {
    Point diff(ptc1.dimension());
    diff = ptc1 - ptc2;
    if ((PIHALF-eval1[1].angle(diff)) < ANGULAR_TOLERANCE &&
	(PIHALF-eval2[1].angle(diff)) < ANGULAR_TOLERANCE)
      istat = 2;    // A minimum distance found
    else
      istat = 3;    // Nothing found
  }
    
}

//===========================================================================
void newPointEval(ParamCurve *pcurve1,double par1,ParamCurve *pcurve2,
		    double astart2,double aend2,double aepsge,
		    double& par2,double& y,double& dist,int& jstat)
//===========================================================================
/*
* PURPOSE     : To evaluate curve 1 in a start point, find the closest point
*               on curve 2, compute the distance vector between the points,
*               the length of the distance vector, and compute the scalar 
*               product of the distance vector and the normal vector on
*               curve 2 in the closest point.
*
* INPUT       : pcurve1  - Curve number one.
*               par1     - Parameter value of start point on curve 1.
*               astart2  - Start of curve 2's parameter interval.
*               aend2    - End of curve 2's parameter interval. 
*               pcurve2  - Curve number two.
*               aepsge   - Geometry resolution
*
*
* OUTPUT      : par2     - Parameter value of closest point on curve 2.
*               y        - Scalar product of the distance vector and the 
*                          normal vector on curve 2 in the closest point
*               dist     - Length of the distance vector.
*               jstat - Status message
*                        0 - OK.
*                       -1 - The curve has a singularity.(Should not happen.)
*                       -2 - We have found an intersection point.
*/
{
  double clo_dist;
  static Point clo_pt(2);
  static Point pt1(2);
  static Point diff(2);
  static Point normal2(2);
  static std::vector<Point> eval2(2);

  jstat=0;

  // Evaluate curve 1 in the start point
  pcurve1->point(pt1,par1);

 // Find closest point on curve 2
  pcurve2->closestPoint(pt1,astart2,aend2,par2,clo_pt,clo_dist);

  // Evaluate curve 2 in the closest point
  pcurve2->point(eval2,par2,1);

  // Difference vector between the points
  diff = eval2[0]-pt1;

  // Normal vector, curve 2  
  normal2[0] = -eval2[1][1];
  normal2[1] =  eval2[1][0];

  // Normalize normal vector to length 1.
  // Return if the curve has a singularity. (Should not happen.)
  if (normal2.normalize_checked() == 0.0) {
    dist=diff.length();
    jstat = -1;
    return;
  }

// Length of difference vector
  dist=diff.length();

// Return if we have found an intersection point.
  if (dist < aepsge) {
    jstat = -2;
    return;
  }

  y = normal2*diff;  
}

}  // End of anonymous namespace for helper functions definition.


