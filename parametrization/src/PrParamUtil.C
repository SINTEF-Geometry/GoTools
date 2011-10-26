#include "GoTools/parametrization/PrParamUtil.h"
//#include "GoTools/utils/Values.h"
#ifdef _WIN32
#define M_PI 3.14159265358979
#endif



#include <math.h>


//-----------------------------------------------------------------------------
void baryCoords(double x, double y, double x0, double y0, double x1, double y1,
                double x2, double y2, double& tau0, double& tau1, double& tau2)
//-----------------------------------------------------------------------------
// Find the three barycentric coordinates of the point (x,y) with respect
// to the triangle with vertices (x0,y0), (x1,y1), (x2,y2).
// The three vertices should not be collinear!
{
  double det0 = det(x1-x,y1-y,x2-x,y2-y);
  double det1 = det(x2-x,y2-y,x0-x,y0-y);
  double det2 = det(x0-x,y0-y,x1-x,y1-y);
  double det3 = det(x1-x0,y1-y0,x2-x0,y2-y0);

  if(det3 <= 0.0) det3 = 1.0;

  tau0 = det0 / det3;
  tau1 = det1 / det3;
  tau2 = det2 / det3;
}

//-----------------------------------------------------------------------------
void baryCoords0(double& u0, double& v0, double& u1, double& v1,
                 double& u2, double& v2, double& tau0, double& tau1, double& tau2)
//-----------------------------------------------------------------------------
//   Find the barycentric coordinates of the origin (0,0) with respect
//   to the triangle formed by three vectors (u0,v0),(u0,v1),(u0,v2).
//   This is bit more efficient than calling  baryCoords.
//   M.F. Feb 97.
{
  double det0 = det(u1,v1,u2,v2);
  double det1 = det(u2,v2,u0,v0);
  double det2 = det(u0,v0,u1,v1);
  double det3 = det(u1-u0, v1-v0, u2-u0, v2-v0);

  if(det3 <= 0.0)
  {
    //printf("Zero determinant in baryTriang.\n");
    //printf("Setting it to 1.0.\n");
    det3 = 1.0;
  }

  tau0 = det0 / det3;
  tau1 = det1 / det3;
  tau2 = det2 / det3;
}

//-----------------------------------------------------------------------------
double det(const double& u1, const double& v1,
         const double& u2, const double& v2)
//-----------------------------------------------------------------------------
{
  return u1 * v2 - v1 * u2;
}

//-----------------------------------------------------------------------------
double area(const double& x0, const double& y0,
          const double& x1, const double& y1,
          const double& x2, const double& y2)
//-----------------------------------------------------------------------------
// Find the signed area of triangle [(x0,y0),(x1,y1),(x2,y2)],
// where (x0,y0),(x1,y1),(x2,y2) are assumed to be ordered anticlockwise.
{
  return 0.5 * det(x1-x0,y1-y0,x2-x0,y2-y0);
}

//-----------------------------------------------------------------------------
double area(const Vector3D& a,
	    const Vector3D& b,
	    const Vector3D& c)
//-----------------------------------------------------------------------------
// Find the positive area of a 3D triangle.
{
    return 0.5 * ((b-a).cross(c-a)).length();
}

//-----------------------------------------------------------------------------
void polarCoords(Vector2D v, double& r, double& theta)
//-----------------------------------------------------------------------------
// Represent a vector v in polar coordinates r(cos(theta),sin(theta))
// where 0 <= r < 2 pi.
{

  r = v.length();
  if(r > 0.0)
  {
    theta = acos(v.x() / r);
#if ((_MSC_VER > 0) && (_MSC_VER < 1300))
    const double PI = 3.14159265358979323846264338327950288419716939937510;
    if(v.y() < 0.0) theta = 2.0 * PI - theta;
#else
    if(v.y() < 0.0) theta = 2.0 * M_PI - theta;
#endif
  }
  else
  {
    theta = 0.0;
  }
}

//-----------------------------------------------------------------------------
void polarCoords(double u, double v, double& r, double& theta)
//-----------------------------------------------------------------------------
// Represent a vector (u,v) in polar coordinates r(cos(theta),sin(theta))
// where 0 <= r < 2 pi.
{
  Vector2D vec(u,v);
  polarCoords(vec,r,theta);
}


//-----------------------------------------------------------------------------
double tanThetaOverTwo(const Vector3D& a,
		       const Vector3D& b,
		       const Vector3D& c)
//-----------------------------------------------------------------------------
// Return tangent of half the angle between vectors b-a and c-a
// without using trig functions.
// Use fact that tan(alpha/2) = (1-cos(alpha)) / sin(alpha).
// and use scalar and dot products to get cos(alpha) and sin(alpha).
//  M.F. Apr. 2002.
{
  Vector3D bb = b - a;
  Vector3D cc = c - a;

  double cp = (bb.cross(cc)).length();
       // length of cross product of bb and cc
  double dp = bb * cc; // scalar product of bb and cc
  double bc = bb.length() * cc.length();
  return (bc - dp) / cp;

  return true;
}

//-----------------------------------------------------------------------------
double cotangent(const Vector3D& a, const Vector3D& b, const Vector3D& c)
//-----------------------------------------------------------------------------
{
    Vector3D ba = a - b;
    Vector3D bc = c - b;
    return (ba*bc)/(ba%bc).length();
}
