#include <vector>
using std::vector;
#include "GoTools/geometry/GeometryTools.h"
#include "GoTools/utils/MatrixXD.h"
#include "GoTools/geometry/ClosestPoint.h"
#include "GoTools/geometry/RectDomain.h"
#include "GoTools/utils/GeneralFunctionMinimizer.h"
#include "GoTools/geometry/GeometryTools.h"


//***************************************************************************
//
// Implementation file of the free function ClosestPoint::closestPtCurveSurf defined in
// ClosestPoint::closestPtCurveSurf.h/
//
//***************************************************************************

using namespace Go;

// Anonymous namespace
namespace {
  const double DZERO = (double)0.0;
  const double TOL = 1.0e-17; //1.0e-16;
  const double REL_COMP_RES = 0.000000000000001;
  const double ANGULAR_TOLERANCE = 0.01;
  const double SINGULAR = 1.0e-16;
}

// Anonymous namespace for definition of class CrvSrfDistFun
namespace {

// Distance function between a curve and a surface  Used by the minimization algorithm
// initiated by ClosestPoint::closestPtCurveSurf
class CrvSrfDistFun {
public:
    CrvSrfDistFun(const ParamSurface* sf, 
		  const ParamCurve* cv,
		  const double* const minpar = 0,
		  const double* const maxpar = 0);
    
    inline double operator()(const double* arg) const;
    inline double grad(const double* arg, double* res) const;
    inline double minPar(int pardir) const;
    inline double maxPar(int pardir) const;

private:
    double minpar_[6];
    double maxpar_[6];
    const ParamSurface * const sf_;
    const ParamCurve * const cv_;
    mutable Point p1_, p2_, d_;
    mutable vector<Point> p1vec_, p2vec_;
};
} // End of anonymous namespace for declaration of class CrvSurfDistFun


// Helper functions declaration
// Anonymous namespace for helper functions declaration.
namespace
{
void insideParamDomain(Point& gd, const Point& acoef, double astart1, 
		       double aend1, double astart2[],double aend2[],
		       int& corr);

void nextStep(double& dist,Point& diff,Point& delta,int& kstat,
	      std::vector<Point>& eval_cv,std::vector<Point>& eval_su, 
	      int order);

int localPretop(double dist,const Point& diff,const Point& normal,
		const std::vector<Point>& eval_cv,
		const std::vector<Point>& eval_su);

void singular(ParamCurve* pcurve,ParamSurface* psurf,Point& par_val,
	      double& dist, double aepsge,double delta,const Point& diff,
	      const Point& norm_vec,
	      const std::vector<Point>&eval_cv,
	      const std::vector<Point>&eval_su,
	      double astart1,double estart2[],double aend1,double eend2[]);

void secant(ParamCurve *pcurve,ParamSurface *psurf,Point& par_val,
	    double& dist,int& jstat, double delta,double aepsge,
	    double astart1,double estart2[],double aend1,double eend2[]);

void setReturnValues(const Point& par_val,ParamCurve* pcurve,
		     ParamSurface* psurf,double aepsge,
		     int& jstat,double& par_cv,double par_su[],
		     double& dist,Point& pt_cv,Point& pt_su);

void newPointEval(ParamCurve *pcurve,double par_cv,ParamSurface *psurf,
		  double estart2[],double eend2[],double aepsge,
		  Point& par_su,double& y,double& dist,int& jstat);

} // End of anonymous namespace for helper functions declaration.


namespace Go {

//===========================================================================
void ClosestPoint::closestPtCurveSurf(ParamCurve* pcurve, ParamSurface* psurf, double aepsge,
			double astart1, double aend1, RectDomain *domain,
		       double anext1, double enext2[],
		       double& cpos1, double gpos2[],
		       double& dist, Point& pt_cv, Point& pt_su,
			bool second_order)
//===========================================================================
{
  int status = 0;
  double rel_comp_res = 1.0e-15;
  double sstart[2], send[2];
  sstart[0] = domain->umin();
  sstart[1] = domain->vmin();
  send[0] = domain->umax();
  send[1] = domain->vmax();
  ClosestPoint::closestPtCurveSurf(pcurve, psurf, aepsge, astart1, sstart, aend1, 
		     send, anext1, enext2, cpos1, gpos2, dist,
		     pt_cv, pt_su, status, second_order);
  DEBUG_ERROR_IF(status < 0,"Error in closest point");

  if (status >= 10 && dist > rel_comp_res && !second_order)
  {
      // Iteration encountered a singularity. Try fallback iteration
    const double TOL = 1.0e-8;
    double seed[3], minpar[3], maxpar[3];
    double dist2;
    minpar[0] = sstart[0];
    minpar[1] = sstart[1];
    minpar[2] = astart1;
    maxpar[0] = send[0];
    maxpar[1] = send[1];
    maxpar[2] = aend1;
    seed[0] = gpos2[0];
    seed[1] = gpos2[1];
    seed[2] = cpos1;

    CrvSrfDistFun distfun(psurf, pcurve, minpar, maxpar);
    FunctionMinimizer<CrvSrfDistFun> funmin(3, distfun, seed, TOL);
    minimise_conjugated_gradient(funmin);//, 3); // number of iterations in each cycle

    dist2 = sqrt(funmin.fval());
    if (dist2 < dist)
    {
	dist = dist2;
	cpos1 = funmin.getPar(2);
	gpos2[0] = funmin.getPar(0);
	gpos2[1] = funmin.getPar(1);

	pcurve->point(pt_cv, cpos1);
	psurf->point(pt_su, gpos2[0], gpos2[1]);

	status = (dist <= aepsge) ? 1 : 2;
    }
    else 
	status = status % 10;
  }

}

// (s1772)
//===========================================================================
void ClosestPoint::closestPtCurveSurf(ParamCurve* pcurve, ParamSurface* psurf, double aepsge,
		       double astart1, double estart2[], double aend1,
		       double eend2[], double anext1, double enext2[],
		       double& cpos1, double gpos2[],
		       double& dist, Point& pt_cv, Point& pt_su, int& istat,
			bool second_order)
//===========================================================================
/*
* PURPOSE    : Newton iteration on the distance function between
*              a curve and a surface to find a closest point
*              or an intersection point.
*              Ported from the sisl-function s1772.
*
*
* INPUT      : pcurve    - Pointer to the curve in the intersection.
*              psurf     - Pointer to the surface in the intersection.
*              aepsge    - Geometry resolution.
*              astart1   - Start parameter value of the curve.
*              estart2[] - Start parameter values of surface.
*              aend1     - End parameter value of the curve.
*              eend2[]   - End parameter values of the surface.
*              anext1    - Start parameter value of the iteration on
*                          the curve.
*              enext2[]  - Start parameter values of the iteration on
*                          the surface.
*
*
* OUTPUT     : cpos1   - Parameter value of the curve in the intersection
*                        point.
*              gpos2[] - Parameter values of the surface in the
*                        intersection point.
*              dist    - Distance between the points.
*              pt_cv   - The point on the curve.
*              pt_su   - The point on the surface.
*              istat   - status messages 
*                                = 1   : Intersection found.
*                                = 2   : A minimum distance found.
*                                = 3   : Nothing found.
*
*
* METHOD     : Newton iteration in tree parameter directions.
*
*********************************************************************
*/ 
{

  const int dim = pcurve->dimension();
  DEBUG_ERROR_IF(dim != psurf->dimension(), "Dimension mismatch.");
  DEBUG_ERROR_IF(dim != 3, "Only implemented for 3D");


 // Parameter interval of surface and curve 
  Point delta(3);
  delta[0]=psurf->containingDomain().umax() - psurf->containingDomain().umin();
  delta[1]=psurf->containingDomain().vmax() - psurf->containingDomain().vmin();
  delta[2]=pcurve->endparam() - pcurve->startparam();


  int knbit;                // Number of iterations.
  int p_dir;                // Changing direction in par-space.
  int g_up,ng_up,g_dir;     // Changing direction in geometric space.
  int order;		    // Order of method.
  int sing = 0;		    // Mark that singularity has occured.
  int left[3];              // Variables used in the evaluator.
  Point d(3);		    // Clipped distances between old and new par.
                            // value in the tree parameter directions.
  Point c_d(3);	            // Computed distances ....
  Point nc_d(3);	    // New computed distances ....
  double prev_dist;         // Previous difference between curve and surface.
  Point par_val(3);         // Parameter values
  Point norm_vec(dim);      // Normal vector to the surface
  Point diff(dim);          // Difference vector between points.
  int corr = 0, div2 = 0;
  std::vector<Point> c0(3,Point(dim)), s0(6,Point(dim));

  par_val[0] = enext2[0];
  par_val[1] = enext2[1];
  par_val[2] = anext1;

  left[0]=left[1]=left[2]=0;

  bool from_right;
  int ki, kstat;


  for (ki=1; ki<3; ki++) {
      if (second_order && ki==1)
	  continue;
    order=ki;

    // Evaluate 0-order derivatives of curve
    from_right = (par_val[2] == aend1) ? false : true;
    pcurve->point(c0, par_val[2], order, from_right);

   // Evaluate 0-order derivatives of surface
    psurf->point(s0, par_val[0], par_val[1], order);

    // Compute the distance vector, value and the new step.
    nextStep(dist,diff,c_d,kstat,c0,s0,order);

    // Surface's normal vector
    norm_vec=s0[1].cross(s0[2]);

    if (kstat == 1) { // Singular matrix
      if (order == 2) {
	singular(pcurve,psurf,par_val,dist,aepsge,c_d[0],diff,norm_vec,c0,
		 s0, astart1,estart2,aend1,eend2);
	setReturnValues(par_val,pcurve,psurf,aepsge,istat,cpos1,gpos2,dist,
			pt_cv,pt_su);
	istat += 10;
	return;
      }
    }
    else break;
  }


  // Adjust the step size if we are not inside the parameter interval

  g_up =  (diff*norm_vec >= DZERO) ? 1 : -1;
  d = c_d;
  insideParamDomain(d,par_val, astart1,aend1,estart2,eend2,corr);

  prev_dist = dist;

  // Start of iteration to find the intersection point
  const int max_it=30;
  for (knbit = 0; knbit < max_it; knbit++) {
    par_val+=d;

    // Test if the current method is OK, or should we change method?
    while (1) {
      // Evaluate derivatives of the curve
      from_right = (par_val[2] == aend1) ? false : true;
      pcurve->point(c0, par_val[2], order, from_right);

      // Evaluate derivatives of the surface
      psurf->point(s0, par_val[0], par_val[1], order);      

      // Compute the distanse vector and value and the new step.
      nextStep(dist,diff,nc_d,kstat,c0,s0,order);
      if (kstat == 1) {
	// Singular matrix.
	sing++;
	if (order==2) {
	  // Singularity (closest point) found.
	  singular(pcurve,psurf,par_val,dist,aepsge,nc_d[0],diff,
		   norm_vec,c0,s0, astart1,estart2,aend1,eend2);
	  setReturnValues(par_val,pcurve,psurf,aepsge,istat,cpos1,gpos2,dist,
			  pt_cv,pt_su);
	  istat += 10;
	  return;
	} 
	else
	  order=2;  // Use more terms in the series expansion
      }
      else {
	norm_vec=s0[1].cross(s0[2]);
	
	// Have normal and difference vectors the same direction ?
	ng_up = (diff*norm_vec >= DZERO) ? 1 : -1;
	
	// g_dir=1 if we have not changed position to the other side
	// of the intersection point. 0 if changed.
	g_dir = (ng_up+g_up != 0);
	
	// p_dir=1 if the steps in the parameter interval continues
	// in the same direction. 0 if direction has changed.
	p_dir = (c_d*nc_d >= DZERO);
	
	if (order!=2 && g_dir && (!p_dir || dist > 0.3*prev_dist)) {
	  if (div2)            // Few terms in the series expansion.
	    div2 = 0;          // Not good enough convergence. 
	  order=2;             // Change method.
	}
	else if (order==2 && !g_dir) { // Max terms in the series expansion.
	  if (sing) {        // Found closest point?
	    singular(pcurve,psurf,par_val,dist,aepsge,nc_d[0],diff,
		     norm_vec,c0,s0, astart1,estart2,aend1,eend2);
	    setReturnValues(par_val,pcurve,psurf,aepsge,istat,cpos1,gpos2,dist,
			    pt_cv,pt_su);
	    istat += 10;
	    return;
	  }
	  if (div2)
	    div2 = 0; // Closest point not found. Change method
	  order=1;             
	}
	else {           // Decreasing distance. 
	  if (sing)      // Continue with current method. 
	    sing = 0;
	  break;
	}
      }
    } // end while(1)
    
 
    // We have decided method and computed a new step.
    // Test if we are inside the intervals, and if the iteration is
    // converging fast enough.
    if (corr)
      if (!(p_dir && g_dir))
	corr = 0;
    
    if (dist < prev_dist) {
      // Adjust if we are not inside the parameter interval.
      if (div2)
	div2 = 0;
      g_up = ng_up;
      d = c_d = nc_d;
      insideParamDomain(d,par_val, astart1,aend1,estart2,eend2,corr);

      prev_dist = dist;


      if (corr > 3) 
	break;
    }
    else if ( corr > 3 ||
	      ((fabs(d[0]/delta[0]) <= REL_COMP_RES) &&
	       (fabs(d[1]/delta[1]) <= REL_COMP_RES) &&
	       (fabs(d[2]/delta[2]) <= REL_COMP_RES)))
      break;
    else {
      // Not converging, adjust and try again.
      if (corr)
	corr++;

      if (dist > prev_dist && div2)
	break;
      div2++;                   // Try half of the step length
      par_val -= d;
      d[0] *= 0.5;  d[1] *=0.5; d[2] *=0.5;
    }

  }

  // end of iteration.


  setReturnValues(par_val,pcurve,psurf,aepsge, istat,cpos1,gpos2,dist,
		  pt_cv,pt_su);
}

} // End of namespace Go


// Anonymous namespace for helper functions definition.
namespace {
// (s1772_s9corr)
//===========================================================================
void insideParamDomain(Point& gd, const Point& acoef, double astart1, 
		       double aend1, double astart2[],double aend2[],
		       int& corr)
//===========================================================================
/*
*                                                                   
* PURPOSE    : To be sure that we are inside the border of the
*              parameter plane. If we are outside clipping is used
*	       to adjust the step value.
*              Ported from the sisl-function s1772_s9corr.
*
*
* INPUT      : acoef     - Parameter values.
*              astart1   - The lower border in curve.
*              aend1     - The higher border in curve.
*              astart2[] - The lower border in surface.
*              aend2[]   - The higher border in surface.
*
*
*
* INPUT/OUTPUT : gd   - Old and new step values.
*                corr - If the step value has been changed, corr is increased
*                       by one. If not, corr is set to zero.
*
*/
{
  int lcorr = 0;

  // Surface. First parameter direction.
  if (acoef[0] + gd[0] < astart2[0]) { 
    gd[0] = astart2[0] - acoef[0]; 
    lcorr=1;
  }
  else if (acoef[0] + gd[0] > aend2[0]) {
    gd[0] = aend2[0] - acoef[0]; 
    lcorr=1;
  }

  // Surface. Second parameter direction.  
  if (acoef[1] + gd[1] < astart2[1]) {
    gd[1] = astart2[1] - acoef[1]; 
    lcorr=1;
  }
  else if (acoef[1] + gd[1] > aend2[1]){
    gd[1] = aend2[1] - acoef[1]; 
    lcorr=1;
  }

  // Curve  
  if (acoef[2] + gd[2] < astart1) {
    gd[2] = astart1 - acoef[2]; 
    lcorr=1;
  }
  else if (acoef[2] + gd[2] > aend1){
    gd[2] = aend1 - acoef[2]; 
    lcorr=1;
  }
  
  if (lcorr) 
     corr++;
  else
    corr = 0;
}


// (s1772_s9dir)
//===========================================================================
void nextStep(double& dist,Point& diff,Point& delta,int& kstat,
	      std::vector<Point>& eval_cv,std::vector<Point>& eval_su, 
	      int order)
//===========================================================================
/*                                                                   
* PURPOSE    : To compute the distance vector and the length of the
*              distance vector, and to compute a next step on all three 
*              parameter directions towards an intersection or a closest point.
*              Ported from the sisl-function s1772_s9dir.
*
*
* INPUT      : eval_cv - Value and derivatives on the curve.
*              eval_su - Value and derivatives on the surface.

*              order - 2 if we have to use second order method 
*		       1 if only first order method.
*
* OUTPUT     : dist  - The length of the distance vector between the points.
*	       diff  - Distance vector between the points.
*	       delta - Relative step parameter values towards intersection on
*                      the curve delta[2] and the surface (delta[0],delta[1]).
*              kstat - kstat=0:Solution found.  kstat=1:Singular system. 
*
*
* METHOD     : We have a point on the curve and a point on the surface.
*	       The distance vector between these two points are going towards
*	       either zero length or to be orthogonal to the derivative on
*	       the curve and the derivatives on the surface.
*	       The method is to use Newton itarations.
*	       The dot product between the distance vector and the derivatives
*	       are going to be zero.
*	       We then use Taylor expansion both on the position and
*	       the derivatives. Then we just remove the second degree parts
*	       of the equations. The matrix is then splitted into a first order
*	       part A:
*		     | -g_u |	
*		 K = | -g_v |    ,and       A = K*K(transpost)
*		     |  f_t | 	
*
*	         where f_t is the first derivative of the curve,
*	         and  g_u is the first derivative of the surface in first dir,
*	         and  g_v is the first derivative of the surface in second dir,
*
*	       and a second order part B:
*		     | -<d,g_uu>   -<d,g_uv> 		 |
*		 B = | -<d,g_uv>   -<d,g_vv>		 |
*		     |  			<d,f_tt> |
*
*	         where d is the distance vector,
*	         and f_tt is second derivative of the curve,
*	         and  g_uu is second derivative of the surface in first dir,
*	         and  g_vv is second derivative of the surface in second dir,
*	         and  g_uv is cross derivative of the surface.
*
*	       We then have the following possible equations:
*
*		             A*delta = -K*d, or    (1)
*	                 (A+B)*delta = -K*d, or    (2)
*	         K(transposed)*delta = -d.         (3)         
*
*              Here we use (1) and (2). In s1772 (1) and (3) were used.
*                            
*	       The solutions of these matrix equations are the
*	       following function.
*
*********************************************************************
*/
{

  MatrixXD<double, 3> A;     // Equation system matrix
  MatrixXD<double, 3> mat;   // Equation system matrix
  double h[3];               // Equation system right hand side

  const Point& f=   eval_cv[0];  // Value in point on curve.
  const Point& f_t= eval_cv[1];  // First derivative of the curve.
  const Point& g=   eval_su[0];  // Value in point on surface.
  const Point& g_u= eval_su[1];  // Derivative of the surface in first par-dir.
  const Point& g_v= eval_su[2];  // Derivative of the surface in second par-dir.

  // Difference vector between the points
  diff = f - g;  // Difference vector between the points

  // Length of the difference vector 
  dist = diff.length();

  // Make equation system. Fills matrix  and right hand side h

  double a1,a2,a3,a4,a5,a6;	// The A matrix, diagonal and A12 A13 A23.
  double b1,b2,b3,b4;		// The B matrix, diagonal and B23.

  // Scalar products
  a1 = f_t*f_t;
  a2 = g_u*g_u;
  a3 = g_v*g_v;
  a4 = f_t*g_u;
  a5 = f_t*g_v;
  a6 = g_u*g_v;
     
  if (order==2) {
    const Point& f_tt=eval_cv[2];  // Second derivative of the curve.
    const Point& g_uu=eval_su[3];  // 2nd deriv of the surface in first par-dir
    const Point& g_uv=eval_su[4];  // Cross derivative of the surface.
    const Point& g_vv=eval_su[5];  // 2nd deriv of the surface in second par-dir
    b1 = diff*f_tt;
    b2 = diff*g_uu;
    b3 = diff*g_vv;
    b4 = diff*g_uv;
  }
  else
    b1=b2=b3=b4=0.0;

  A(0,0) = a2-b2;	A(0,1) = a6-b4;		A(0,2) = -a4;
  A(1,0) = a6-b4;	A(1,1) = a3-b3;		A(1,2) = -a5;
  A(2,0) = -a4;	        A(2,1) = -a5;		A(2,2) = a1+b1;
  
  h[0] =  diff*g_u;
  h[1] =  diff*g_v;
  h[2] = -diff*f_t;


  double det = A.det();  // Determinant

  if (fabs(det) < TOL) {   // @bsp
    kstat = 1;
    delta[0] = 0.0;
    delta[1] = 0.0;
    delta[2] = 0.0;
  }
  else {   // Using Cramer's rule to find the solution of the 3x3 system
    int i;
    kstat = 0;
    for(i=0;i<3;i++) {
      mat(i,0)=h[i];
      mat(i,1)=A(i,1);
      mat(i,2)=A(i,2);
    }
    delta[0] = mat.det()/det;

    for(i=0;i<3;i++) {
      mat(i,0)=A(i,0);
      mat(i,1)=h[i];
    }
    delta[1] = mat.det()/det;

    for(i=0;i<3;i++) {
      mat(i,1)=A(i,1);
      mat(i,2)=h[i];
    }
    delta[2] = mat.det()/det;
  }

}



// (s1772_s6local_pretop)
//===========================================================================
int localPretop(double dist,const Point& diff,const Point& normal,
		const std::vector<Point>& eval_cv,
		const std::vector<Point>& eval_su)
//===========================================================================
/*
*   PURPOSE : To find if we have a minimum or a maximum or a point of 
*             inflection situation. This function assumes that it is a 
*             singular situation.
*             Ported from the sisl function s1772_s6local_pretop.
*
*
* INPUT      : dist    - The length of the difference vector.
*	       diff    - The difference vector.
*	       normal  - The normal vector on the surface.
*              eval_cv - Value and derivatives in the point on the curve.
*              eval_su - Value and derivatives in the point on the surface
*   
*   OUTPUT  :  return value - 	= -1: degenerated system.
*				=  0: Maximum position or inflection point.
*				=  1: Minimum position.
*
*
*   METHOD  : Computing and interpretation of curvatures.
*
************************************************************************
*/
{
  const int dim = diff.dimension();
  DEBUG_ERROR_IF(dim != 3, "Only implemented for 3D");

  //  int kstat = 0;	// Status variable.
  double a1,a2,a3,a4;   // Matrix.

  Point S_u(dim);       // Normalized s_u.
  Point S_v(dim);	// Normalized s_v.
  Point S_uxS_v(dim);	// Cross between S_u and S_v.
  Point s_d(dim);	// Second derivative in direction f_t.
  Point N(dim);		// Normalized normal.
  Point d_uv(2);	// Normalized direction vector in par-plane.

  //  const Point& f=   eval_cv[0];  // Value in point on curve.
  const Point& f_t= eval_cv[1];  // First derivative of the curve.
  const Point& f_tt=eval_cv[2];  // Second derivative of the curve.

  //  const Point& s=   eval_su[0];  // Value in point on surface.
  const Point& s_u= eval_su[1];  // Derivative of the surface in first par-dir.
  const Point& s_v= eval_su[2];  // Derivative of the surface in second par-dir.
  const Point& s_uu=eval_su[3];  // 2nd deriv. of the surface in first par-dir.
  const Point& s_uv=eval_su[4];  // Cross derivative of the surface.
  const Point& s_vv=eval_su[5];  // 2nd deriv. of the surface in second par-dir.


  if (diff.angle(normal) > ANGULAR_TOLERANCE)
    return -1;

  S_u=s_u;
  if (S_u.normalize_checked() == 0.0)
    return -1;

  S_v=s_v;
  if (S_v.normalize_checked() == 0.0)
    return -1;

  S_uxS_v= S_u.cross(S_v);
  a1 = S_u*S_v;
  a2 = f_t*S_u;
  a3 = f_t*S_v;
  a4 = S_uxS_v*S_uxS_v;
  if (a4 < SINGULAR)
    return -1;
  
  d_uv[0] = (a2 - a1*a3)/a4;
  d_uv[1] = (a3 - a1*a2)/a4;
  if (d_uv.normalize_checked() == 0.0)
    return -1;

  a1 = d_uv[0]*d_uv[0];
  a2 = d_uv[1]*d_uv[1];
  a3 = 2*d_uv[0]*d_uv[1];

  int ki;
  for (ki=0; ki<dim; ki++)
     s_d[ki] = a1*s_uu[ki] + a3*s_uv[ki] + a2*s_vv[ki];
  
  for (ki=0; ki<dim; ki++)
     N[ki] = diff[ki]/dist;  
  
  a1 = N*f_tt - N*s_d;

  return  (a1 > 1.0e-10);

}

//===========================================================================
void singular(ParamCurve* pcurve,ParamSurface* psurf,Point& par_val,
	      double& dist, double aepsge,double delta,const Point& diff,
	      const Point& norm_vec,
	      const std::vector<Point>&eval_cv,
	      const std::vector<Point>&eval_su,
	      double astart1,double estart2[],double aend1,double eend2[])
//===========================================================================
/*
*   PURPOSE : Singularity (closest point) found. If we have a point of 
*             inflection or a maximum position situation call the
*             secant method.
*
*
*   INPUT   : pcurve   - Pointer to the curve.
*	      psurf    - Pointer to the surface.
*             aepsge   - Geometry resolution.
*             delta    - Parameter distance on the curve between start values
*             diff     - Distance vector between the points.
*             norm_vec - The normal vector on the surface.
*	      eval_cv  - Value and derivatives of the curve.
*             eval_su  - Value and derivatives of the surface.
*
*   INPUT/OUTPUT : par_val - Parameter values for the curve and the surface.
*	           dist    - Distance between the points.
*
*
*   METHOD  : Uses function "localPretop" to test if we have a point
*             of inflection or a maximum position situation.
*
************************************************************************
*/
{
  if (dist > aepsge)
  {
    int ki=localPretop(dist,diff,norm_vec,eval_cv,eval_su);
    if (ki == 0) {
      int kstat;
      secant(pcurve,psurf,par_val,dist,kstat, delta,aepsge,
	     astart1,estart2,aend1,eend2);
    }
  }
}

// (s1772_s6sekant1)
//===========================================================================
void secant(ParamCurve *pcurve,ParamSurface *psurf,Point& par_val,
	    double& dist,int& jstat, double delta,double aepsge,
	    double astart1,double estart2[],double aend1,double eend2[])
//===========================================================================
/*
* PURPOSE    : Secant method iteration on the distance function between
*              a curve and a surface to find a closest point
*              or an intersection point.
*              Ported from the sisl function s1772_s6sekant1
*
*
* INPUT      : pcurve    - Pointer to the curve in the intersection.
*              psurf     - Pointer to the surface in the intersection.
*              delta     - Parameter distance on the curve between start values
*              aepsge    - Geometry resolution.
*              astart1   - Start parameter value of the curve.
*              estart2[] - Start parameter values of surface.
*              aend1     - End parameter value of the curve.
*              eend2[]   - End parameter values of the surface.
*
*
*
* INPUT/
* OUTPUT     : par_val   - Parameter value of the surface in
*                          intersection point.
*              dist      - Distance in space between the points.
* OUTPUT     : jstat     - status messages  
*                                = 0   : OK
*                                = 2   : Can't find a reasonable start point
*
*
* METHOD     : Secant method in three parameter directions.
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
  const int dim = pcurve->dimension();
  DEBUG_ERROR_IF(dim != psurf->dimension(), "Dimension mismatch.");
  DEBUG_ERROR_IF(dim != 3, "Only implemented for 3D");


  std::vector<Point> s0(2,Point(dim));  // Derivatives of the curve
  std::vector<Point> c0(3,Point(dim));  // Derivatives of the surface
  Point norm(dim);       // Normal vector
  Point diff(dim);       // Difference vector between point/surface.
  Point pt(dim);         // Point for use in closest point point/surface.
  Point clo_pt(dim);     // Point for use in closest point point/surface.
  //  double pt_dist;        // Distance between points;

  if (delta == 0.0) 
    delta =1e-15;

  if ((par_val[2] == astart1 && delta < 0.0) ||
      (par_val[2] == aend1   && delta > 0.0)) {
    delta = -delta;       // Change direction of search along curve
    shift++;
  }



  // Make a "segmentation" of the parameter interval on the curve to get
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



  cu_val[0] = par_val[2];
  newPointEval(pcurve,cu_val[0],psurf,estart2,eend2,aepsge,
		 par_val,y[0],dist,kstat);
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

    newPointEval(pcurve,cu_val[1],psurf,estart2,eend2,aepsge,
		   par_val,y[1],dist,kstat);
    if (kstat < 0) {
      par_val[2]=cu_val[1];
      return;
    }
 
    new_y = y[1]/y[0];
    if (new_y > 1.0000000000001) {
      // We are closer to a singularity. (Closest point, not intersection.)
      if (shift) {
	// We have already changed the direction. Give up.
 	par_val[2] = cu_val[1];
	return;
      }

      // Search in the other direction
      delta = -delta;
      cu_val[1] = cu_val[0] + delta;
      if (cu_val[1] < astart1)
	cu_val[1] = astart1;
      else if (cu_val[1] > aend1)
	cu_val[1] = aend1;
      shift++;
    }

    // Start points on different side of the intersection point or far from
    // a singularity?
    else if (y[0]*y[1] <= 0.0 || fabs(new_y) < 0.5)
      break;
    else {
      // try another point as start point number two
      if (cu_val[1]+delta <= aend1 &&
	  cu_val[1]+delta >= astart1)
	cu_val[1] += delta;
      else if (cu_val[1] < aend1)
  	cu_val[1] = aend1;
      else if (cu_val[1] > astart1)
	cu_val[1] = astart1;
      else {
	par_val[2] = cu_val[1];
	return;
      }
    }
  }


  // Can't find a reasonable start point
  if (ki == 20) {
    jstat = 2;
    return;
  }

  new_cu_val=par_val[2];


  //   Here starts the real search for the intersection point

  for (knbit=0; knbit < 50; knbit++) {
    delta_y = y[0]-y[1];
    if (fabs(delta_y) < REL_COMP_RES)  
      break;   // Start points are too close

    // Fra 2D :
    // Are the start points on different side of the intersection point and
    // is the distance from one start point to the intersection point much
    // smaller than the other?
    //    if (y[0]*y[1] < 0.0 &&
    //	(fabs(y[0]) < 6*fabs(y[1]) || fabs(y[1]) < 6*fabs(y[0])))
    //  new_cu_val = 0.5*(cu_val[1]+cu_val[0]);  // Yes, try the midpoint
    //else
      
    new_cu_val = cu_val[1] + y[1]*(cu_val[1]-cu_val[0])/delta_y;  // secant

    // Make sure we are inside the interval
    if (new_cu_val >= aend1) {
      new_cu_val = aend1;
      if (cu_val[0] == aend1 || cu_val[1] == aend1) {
	par_val[2] = new_cu_val;
	return;
      }
    }
    else if (new_cu_val <= astart1) {
      new_cu_val = astart1;
      if (cu_val[0] == astart1 || cu_val[1] == astart1) {
	par_val[2] = new_cu_val;
	return;
      }
    }


    // Evaluation and closest point calls in the new test point
    newPointEval(pcurve,new_cu_val,psurf,estart2,eend2,aepsge,
		   par_val,new_y,dist,kstat);
    if (kstat < 0) {
      par_val[2]=new_cu_val;
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
	  if (new_y >=  y[0])
	    break;
	  cu_val[0] = new_cu_val;
	  y[0] = new_y;
	}
	else {
	  if (new_y >=  y[1]) 
	    break;
	  cu_val[1] = new_cu_val;
	  y[1] = new_y;
	}
	
      }
      else if (y[0] < 0.0) {
	if (y[0] < y[1]) {
	  if (new_y <=  y[0]) 
	    break;
	  cu_val[0] = new_cu_val;
	  y[0] = new_y;
	}
	else {
	  if (new_y <=  y[1]) 
	    break;
	  cu_val[1] = new_cu_val;
	  y[1] = new_y;
	}
      }
    }
  }
  par_val[2] = new_cu_val;

}


//===========================================================================
void setReturnValues(const Point& par_val,ParamCurve* pcurve,
		     ParamSurface* psurf,double aepsge,
		     int& istat,double& par_cv,double par_su[],
		     double& dist,Point& pt_cv,Point& pt_su)
//===========================================================================
/*
* PURPOSE    :  To set the values to be returned from the function 
*               ClosestPoint::closestPtCurveSurf.
*
*
* INPUT      : par_val  - Parameter values for the surface and the curve.
*              pcurve   - Pointer to the curve.
*              psurf    - Pointer to the surface
*              aepsge   - Geometry resolution.
*
*
* OUTPUT     : istat    - Status variable.
*                       	= 1 : Intersection.
*                       	= 2 : Minimum value.
*                       	= 3 : Nothing found.
*              par_cv   - Parameter value of the curve's point.
*              par_su   - Parameter value of the surfaces's point.
*              dist     - Distance in space between the points.
*              pt_cv    - The point on the curve.
*              pt_su    - The point on the surface.
*
************************************************************************
*/
{
  const double PIHALF = 1.57079632679489661923;

  par_su[0] = par_val[0];
  par_su[1] = par_val[1];
  par_cv    = par_val[2];

  // Compute the points and the distance between them.
  pcurve->point(pt_cv,par_cv); 
  
  psurf->point(pt_su,par_su[0],par_su[1]);

  dist=pt_cv.dist(pt_su);
 
  if (dist <= aepsge)
    istat=1;        // Intersection
  else {
    Point tang_vec(3);     // Curve's tangent vector
    Point norm_vec(3);     // Surface's normal vector

    std::vector<Point> evalc(2), evals(3);

    pcurve->point(evalc,par_cv,1);
    tang_vec=evalc[1];

    psurf->point(evals,par_su[0],par_su[1],1);
    norm_vec=evals[1].cross(evals[2]);

    if ((PIHALF-tang_vec.angle(norm_vec)) < ANGULAR_TOLERANCE)
      istat = 2;    // A minimum distance found
    else
      istat = 3;    // Nothing found
  }
    
}


//===========================================================================
void newPointEval(ParamCurve *pcurve,double par_cv,ParamSurface *psurf,
		    double estart2[],double eend2[],double aepsge,
		    Point& par_su,double& y,double& dist,int& jstat)
//===========================================================================
/*
* PURPOSE  : To evaluate the curve in a start point, find the closest point on
*            the surface, compute the distance vector between the points,
*            the length of the distance vector, and compute the scalar product
*            of the distance vector and the normal vector on the surface
*            in the closest point.
*
* INPUT    : pcurve   - Pointer to the curve.
*            par_cv   - Parameter value of start point on the curve.
*            psurf    - Pointer to the surface.
*            estart2  - Start of the surface's parameter interval.
*            eend2    - End of the surface's parameter interval. 
*            aepsge   - Geometry resolution
*
*
* OUTPUT   : par_su   - Parameter value of closest point on the surface.
*            y        - Scalar product of the distance vector and the 
*                       normal vector on the surface in the closest point.
*            dist     - Length of the distance vector.
*            jstat - Status message
*                     0 - OK.
*                    -1 - The curve has a singularity.(Should not happen.)
*                    -2 - We have found an intersection point.
*/
{
  double clo_u, clo_v, clo_dist;
  static double seed_par[2];
  static Point clo_pt(3);
  static Point pt1(3);
  static Point diff(3);
  static Point normal2(3);
  static std::vector<Point> eval_su(5);

  Vector2D corner1(estart2[0],estart2[1]);
  Vector2D corner2(eend2[0],eend2[1]);
  static RectDomain rect_dom(corner1,corner2);

  jstat=0;

  // Evaluate curve in the start point
  pcurve->point(pt1,par_cv);

 // Find closest point on surface
  seed_par[0]=par_su[0];
  seed_par[1]=par_su[1];
  psurf->closestPoint(pt1,clo_u,clo_v,clo_pt,clo_dist,aepsge,&rect_dom,seed_par);
  par_su[0] = clo_u;
  par_su[1] = clo_v;

  // Evaluate in the closest point
  psurf->point(eval_su,clo_u,clo_v,1);

  // Difference vector between the points
  diff = eval_su[0]-pt1;

  // Surface normal vector   
  psurf->normal(normal2,clo_u,clo_v);

  // Normalize normal vector to length 1.
  // Return if the surface has a singularity. (Should not happen.)
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


} // Anonymous namespace



namespace {

//===========================================================================
CrvSrfDistFun::CrvSrfDistFun(const ParamSurface* sf, 
		       const ParamCurve* cv,
		       const double* const minpar,
		       const double* const maxpar)
//===========================================================================
    : sf_(sf), cv_(cv), p1vec_(3), p2vec_(2)
{
    RectDomain dom = sf->containingDomain();
    if (!minpar) {
	minpar_[0] = dom.umin();
	minpar_[1] = dom.vmin();
	minpar_[2] = cv_->startparam();
    } else {
	minpar_[0] = minpar[0];
	minpar_[1] = minpar[1];
	minpar_[2] = minpar[2];
    }
    if (!maxpar) {
	maxpar_[0] = dom.umax();
	maxpar_[1] = dom.vmax();
	maxpar_[2] = cv_->endparam();
    } else {
	maxpar_[0] = maxpar[0];
	maxpar_[1] = maxpar[1];
	maxpar_[2] = maxpar[2];
    }
}

//===========================================================================    
double CrvSrfDistFun::operator()(const double* arg) const
//===========================================================================
{
    sf_->point(p1_, arg[0], arg[1]);
    cv_->point(p2_, arg[2]);
    return p1_.dist2(p2_);
}

//===========================================================================
double CrvSrfDistFun::grad(const double* arg, double* res) const
//===========================================================================
{
    sf_->point(p1vec_, arg[0], arg[1], 1);
    cv_->point(p2vec_, arg[2], 1);
    d_ = p1vec_[0] - p2vec_[0];
    
    res[0] = 2 * d_ * p1vec_[1];
    res[1] = 2 * d_ * p1vec_[2];
    res[2] = -2 * d_ * p2vec_[1];
    
    return d_.length2();
}

//===========================================================================
double CrvSrfDistFun::minPar(int pardir) const
//===========================================================================
{
    //ASSERT(pardir == 0 || pardir == 1);
    return minpar_[pardir];
}

//===========================================================================
double CrvSrfDistFun:: maxPar(int pardir) const
//===========================================================================
{
    //ASSERT(pardir == 0 || pardir == 1);
    return maxpar_[pardir];
}

} // End of anonymous namespace for helper functions definition.
