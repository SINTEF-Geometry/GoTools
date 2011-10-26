//===========================================================================
//                                                                           
// File: CurvatureUtils.C                                                    
//                                                                           
// Created: 060603                                        
//                                                                           
// Author: Vibeke Skytt
//                                                                           
// Revision: 
//                                                                           
// Description:
//                                                                           
//===========================================================================


#define  _USE_MATH_DEFINES 
#include "GoTools/utils/CurvatureUtils.h"
#include <math.h>

using namespace std;

namespace Go
{
  // Tolerance global to the file to avoid dividing by zero
  double ltol = 1.0e-11;

//===========================================================================
  double curvatureRadius(const std::vector<Point>& der, 
			 std::vector<Point>& unitder)
  // Given position, first and second derivative
  // of a curve passing through a point, compute
  // the unit tangent, curvature vector and curvature 
  // radius of this curve
//===========================================================================
{
  
//     Let c = c(w) be a parameterized curve.
//     The curvature vector is defined as the derivative of the unit tangent
//     vector with respect to the arc length a. If we don't have an arclength
//     parametrization then this parametrization can be written as a function
//     of the arc length w = w(a). By using the kernel rule for differentiation
//     we get:
   
//            d            d       dw   d    c'(w)    dw   d    c'(w)      da
//     k(a) = -- T(w(a)) = -- T(w) -- = -- ---------- -- = -- ---------- / --
//            da           dw      da   dw sqrt(c'c') da   dw sqrt(c'c')   dw
      
//            d       c'(w)                c"        c' (c'c'')
//            -- ----------------- =   ---------- - ------------- 
//            dw sqrt(c'(w) c'(w))     sqrt(c'c')   sqrt(c'c')**3
     
//            da
//            -- = sqrt(c'c')
//            dw 

  // Assuming 2D or 3D
  ASSERT(der.size() >= 3);
  ASSERT(der[0].dimension() == 2 || der[0].dimension() == 3);

  // Initiate output vector
  int ki;
  unitder.resize(3);
  for (ki=0; ki<3; ki++)
    unitder[ki] = der[ki];

  double length = unitder[1].length();
  if (length < ltol)
      return -1;

  unitder[1].normalize();
  
  // Check tangent length
  if (length < ltol)
    THROW("Tangent of length zero");
   
  // Make curvature vector
  double dum1 = (der[2] * unitder[1])/length;
  unitder[2] = (der[2]/length - unitder[1]*dum1)/length;

  // Make curvature radius
  double dum2 = unitder[2].length();
  if (dum2 < ltol)
    return -1;
  else
    return 1.0/dum2;
}

//===========================================================================
  double stepLenFromRadius(double radius, double aepsge)

// Computes the step length along a curve based on radius of curvature
// at a point on the curve, and an absolute tolerance.
//===========================================================================
{
  double tstep;

  if (radius > ltol) {
    // Estimate the opening angle of the segments based on the error formula.
    double talpha = M_PI/0.4879*pow(aepsge/radius,1.0/6.0);

    // Estimate step length equal to curve length of this circular arc,
    // We limit the step length to half the radius of curvature.
    tstep = std::min(talpha,0.5)*radius;
  }

  else if (radius >= 0.0)  //  Radius of curvature is zero 
    tstep = 100.0*aepsge;

  else       //  Infinite radius of curvature
    tstep = 1.e4*aepsge;   // @@ ????  tstep = amax in s1311

  return tstep;
}

//===========================================================================
  double tanLenFromRadius(double radius, double angle)

// To create the tangent length for interpolating a
// circular arc with an almost equi-oscillating Hermit qubic
//===========================================================================
{
  double tcos,tsin;          /* Dummy variables                     */
  double ta,tb,tc,tl;        /* Dummy variables                     */
  double tconst = (double)1.85530139760811990992528773586425;
                             /* Constant used in the calculation    */
  
  
  
  tcos = cos(angle);
  tsin = sin(angle);
  
  /*  Calculate length of tangents
   *   tconst = (3-2sqrt(2))**1/3 + (3+2sqrt(2))**1/3 - 0.5 */
  
  ta     = (double)0.6*tconst - (double)0.9*tcos;
  tb     = ((double)0.4*tconst+(double)1.8)*tsin;
  tc     = ((double)0.4*tconst+(double)1.0)
           * tcos - (double)0.4*tconst - (double)1.0;
  tl     = radius*(-tb+sqrt(tb*tb-4*ta*tc))/((double)2.0*ta);
  
  return(tl);
}


//===========================================================================
void getHermiteData(const vector<Point>& der1, const vector<Point>& der2, 
		    double& parint, double& len1, double& len2)

// Given position, first and second derivative in both ends of
// an Hermite segment, compute parameter interval and tangent lengths
// in order to stay close to a circular segment
//===========================================================================
{
  // First compute unit tangent, curvature and curvature radius
  vector<Point> unitder1, unitder2;
  double crad1, crad2;
  crad1 = curvatureRadius(der1, unitder1);
  crad2 = curvatureRadius(der2, unitder2);
  
  // Compute the angle between the unit tangents
  double angle = unitder1[1].angle(unitder2[1]);
  angle = fabs(angle);

  // Compute distance between endpoints
  parint = der1[0].dist(der2[0]);

  // Compute tangent lengths
  if (angle < ltol || crad1 < 0.0)
    len1 = parint/3.0;
  else
    len1 = tanLenFromRadius(crad1, angle);

  if (angle < ltol || crad2 < 0.0)
    len2 = parint/3.0;
  else
    len2 = tanLenFromRadius(crad2, angle);

  // Make sure that the tangent does not explode due to numeric errors,
  // and make a controlled tangent when the radius is zero or almost zero
  double max_dist;
  if (angle < 0.1)		
    max_dist = (double)0.35*parint;
  else if (angle < 0.35)	
    max_dist = (double)0.40*parint;
  else if (angle < 0.75)	
    max_dist = (double)0.50*parint;
  else 			
    max_dist = (double)0.70*parint;

  if ( len1 > max_dist) 
    len1 = max_dist;
  if ( len2 > max_dist) 
    len2 = max_dist;

}

} // end namespace Go


