#include <vector>
using std::vector;
#include "GoTools/utils/Point.h"
#include "GoTools/geometry/GeometryTools.h"
#include "GoTools/geometry/ParamCurve.h"
#include "GoTools/geometry/ClosestPoint.h"

//***************************************************************************
//
// Implementation file of the free function ClosestPoint::closestPtSurfSurfPlaneGeometrical
// declared in ClosestPoint.h
//
//***************************************************************************

// Anonymous namespace
namespace {
  const double DZERO = (double)0.0;
  const double TOL = 1.0e-16;
  const double REL_PAR_RES = (double)0.000000000001;
}

// Anonymous namespace for helper function declaration.
namespace
{
    void nextStep(const std::vector<Go::Point>& fpnt,const std::vector<Go::Point>& gpnt,
		  const Go::Point& snorm, const Go::Point& spoint,
		  Go::Point& delta, int& kstat);
}

namespace Go {

// (s9iterate)
//===========================================================================
void ClosestPoint::closestPtSurfSurfPlaneGeometrical(const std::vector<Point>& epoint,
				       const std::vector<Point>& epnt1,
				       const std::vector<Point>& epnt2,
				       const Point& epar1,
				       const Point& epar2,
				       const ParamSurface* psurf1,
				       const ParamSurface* psurf2,
				       double aepsge,
				       std::vector<Point>& gpnt1,
				       std::vector<Point>& gpnt2,
				       Point& gpar1, Point& gpar2, int& jstat)
//===========================================================================
/*
* PURPOSE    : To iterate to an intersection point between two surfaces
*              and a plane.
*              Ported from the sisl function s9iterate.
*
*
*
* INPUT      : epoint - Vector of Points containing parts of plane description.
*                       epoint[0] contains a point in the plane.
*                       epoint[1] contains the normal vector to the plane
*              epnt1  - 0-2 Derivatives + normal of start point for
*                       iteration in first surface
*              epnt2  - 0-2 Derivatives + normal of start point for
*                       iteration in second surface
*              epar1  - Parameter pair of start point in first surface
*              epar2  - Parameter pair of start point in second surface
*              psurf1 - Pointer to the first surface
*              psurf2 - Pointer to the second surface
*              aepsge - Absolute tolerance
*
*
* OUTPUT     : gpnt1  - 0-2 Derivatives + normal of result of iteration
*                       in first surface
*              gpnt2  - 0-2 Derivatives + normal of result of iteration
*                       in second surface
*              gpar1  - Parameter pair of result of iteration in first surface
*              gpar2  - Parameter pair of result of iteration in second
*                       surface
*              jstat  - status messages  
*                       = 3      : Iteration diverged or too many iterations
*                       = 2      : iteration converged, singular point found
*                       = 1      : ok, iteration converged
*
*
* METHOD     :
*
* USE        : The function is only working in 3-D
*

*
*********************************************************************
*/
{
    const int kdim = psurf1->dimension();
    DEBUG_ERROR_IF(kdim != psurf2->dimension(), "Dimension mismatch.");
    DEBUG_ERROR_IF(kdim != 3, "The function is only working i 3-D");



    int kcont;              /* Indicator telling if iteration is not finished */
    int kder = 2;           /* Derivative indicator                           */
    int kstat1;             /* Status variable                                */
    int kstat2;             /* Status variable                                */
    int knbit;              /* Counter for number of iterations               */
    int kmaxit = 100;       /* Maximal number of iterations allowed           */
    Point sdiff(kdim);      /* Difference between two vectors                 */
    double tdum3;           /* Dummy variables                                */
    double tdist=0;         /* Distance between two points in iteration       */
    Point normal_vec(3);    /* Surface normal vector                          */
    Point delta(2);         /* Parameter step values in the iteration         */
  
    // Make description of intersection plane

    const Point& spoint = epoint[0];    // Point in the intersection plane.
    const Point& snorm  = epoint[1];    // Normal vector of intersection plane

    // determining parameter domain (only rectangular domain - trimmed surfaces not supported
    const RectDomain pdom1 = psurf1->containingDomain();
    const RectDomain pdom2 = psurf2->containingDomain();
    const Vector2D lower_left_1 = pdom1.lowerLeft();
    const Vector2D lower_left_2 = pdom2.lowerLeft();
    const Vector2D upper_right_1 = pdom1.upperRight();
    const Vector2D upper_right_2 = pdom2.upperRight();

    // Copy input variables to output variables

    for (int i=0; i<7; i++) {
	gpnt1[i] = epnt1[i];
	gpnt2[i] = epnt2[i];
    }

    gpar1 = epar1;
    gpar2 = epar2;

    // new values of gpnt1 and gpnt2
    vector<Point> gpnt1_new(7);
    vector<Point> gpnt2_new(7);
    Point gpar1_new, gpar2_new;

    // At the start of the iteration the two points gpnt1 and gpnt2 might be very
    // close since we in most cases start from a point on the intersection curve.
  
    kcont = 1;
    knbit = 0;
  
    while (kcont) {
    
	// Put a parametric representation of the tangent plane of surface 1
	// into the implicit representation of the tangent plane of surface 2 
	// and also into the implicit representation of the intersection plane.
	// The solution of the equation system is the next step along
	// the two parameter directions of surface 1.    

	nextStep(gpnt1,gpnt2,snorm,spoint, delta,kstat1);
	gpar1_new = gpar1 + delta;
	//gpar1 += delta;	
	
	// Put a parametric representation of the tangent plane of surface 2
	// into the implicit representation of the tangent plane of surface 1 
	// and also into the implicit representation of the intersection plane.
	// The solution of the equation system is the next step along
	// the two parameter directions of surface 2.

	nextStep(gpnt2,gpnt1,snorm,spoint, delta,kstat2);
	gpar2_new = gpar2 + delta;
	//gpar2 += delta;

	// domain checks, kick parameters back into domain!  @@ does this always work?
	for (int i = 0; i < 2; ++i) {
	    gpar1_new[i] = std::max(gpar1_new[i], lower_left_1[i]);
	    gpar1_new[i] = std::min(gpar1_new[i], upper_right_1[i]);
	    gpar2_new[i] = std::max(gpar2_new[i], lower_left_2[i]);
	    gpar2_new[i] = std::min(gpar2_new[i], upper_right_2[i]);
	}       

	//  Calculate values of new points and normal vector on surface 1.

	psurf1->point(gpnt1_new, gpar1_new[0], gpar1_new[1], kder);
	psurf1->normal(normal_vec, gpar1_new[0], gpar1_new[1]);
	gpnt1_new[6]=normal_vec;

	//   If the surface normal has zero length no use in continuing
	if (normal_vec.length() == 0.0) {
	    jstat=3;
	    break;
	}      

	//  Calculate values of new points and normal vector on surface 2.

	psurf2->point(gpnt2_new, gpar2_new[0], gpar2_new[1], kder);
	psurf2->normal(normal_vec, gpar2_new[0], gpar2_new[1]);
	gpnt2_new[6]=normal_vec;      

	//  If the surface normal has zero length no use in continuing
	if (normal_vec.length() == 0.0) {
	    jstat=3;
	    break;
	}


	// Make difference between the two points, 
	//   and calculate length of difference
	sdiff = gpnt1_new[0];
	sdiff -= gpnt2_new[0];
    
	tdum3 = sdiff.length();      
	if (tdum3 < TOL) {
	    // Length is zero. Iteration has converged.
	    kcont = 0;
	    jstat = 1;
	    break;
	}
      
	if (knbit==0) {
	    // First iteration inititate distance variable, if the equation
	    // systems were not singular
      
	    if (kstat1 || kstat2) {
		jstat=3;
		break;
	    }
	    tdist = tdum3;
	    knbit = 1;
	} else {
	    // More than one iteration done, stop if distance is not decreasing.
	    // Then decide if we converge distance between the points is within
	    // the tolerance and the last step had singular or none singular
	    // equation systems.
      
	    knbit = knbit + 1;
	    if (tdum3 >= tdist) {
		// Distance is not decreasing
		if (tdist <= aepsge) {
		    // Distance within tolerance
		    if (kstat1 || kstat2) {
			// Singular equation system
			jstat=2;
			return; // return without swapping new results to old
			//break;
		    }
		    else {
			// Nonsingular equation system
			jstat=1;
			return; // return without swapping new results to old
			//break;
		    }
		}
		else {
		    // Distance is not within tolerance, divergence
		    jstat=3;
		    return; // return without swapping new results to old
		    //break;
		}
	    }
      	    //      Distance still decreasing
	    tdist = tdum3;
	}
      
	//  Make sure that not too many iterations are being done
	if (knbit > kmaxit) {
	    jstat=3;
	    break;
	}
	gpnt1_new.swap(gpnt1);
	gpnt2_new.swap(gpnt2);
	// Hmm, seems to need an old-fashioned swap here in order to
	// get rid of and internal compiler error on gcc 4.2.4 in
	// optimized mode. @jbt
	Point tmp = gpar1_new;
	gpar1_new = gpar1;
	gpar1 = tmp;
	tmp = gpar2_new;
	gpar2_new = gpar2;
	gpar2 = tmp;
// 	gpar1_new.swap(gpar1);
// 	gpar2_new.swap(gpar2);
    }  // while (kcont)

    // If we got here, we want to keep the newly calculated values before returning.
    // Effectuating swap.
    gpnt1_new.swap(gpnt1);
    gpnt2_new.swap(gpnt2);
    gpar1_new.swap(gpar1);
    gpar2_new.swap(gpar2);

}
} // namespace Go  

	
// Anonymous namespace for helper functions definition.
namespace {
using namespace Go;
//===========================================================================
void nextStep(const std::vector<Point>& fpnt,const std::vector<Point>& gpnt,
	      const Point& snorm, const Point& spoint,
	      Point& delta, int& kstat)
//===========================================================================
/*
*                                                                   
* PURPOSE    : To compute the next step values along the parameters of a 
*              surface in an iteration to an intersection point between 
*              two surfaces and a plane.
*
*
* INPUT        fpnt   - 0-2 Derivatives + normal of start point for
*                       iteration in first surface
*              gpnt   - 0-2 Derivatives + normal of start point for
*                       iteration in second surface
*              snorm  - Normal vector of the intersection plane
*              spoint - Point in the intersection plane
*
*
* OUTPUT     : delta  - New step values for first surface.
*              kstat  - status messages  
*                       = 1      : The equation system is singular.
*                       = 0      : The equation system is not singular.
*
*
* METHOD     : On the first surface we have a point f and the derivatives
*              in that point f_u and f_v in the two parameter directions
*              u and v. On the second surface we have a point g with
*              derivatives g_s and g_t and normal vector g_n.
*              In the intersecting plane we have a point p and the plane's
*              normal vector p_n. With the distance vector between the two
*              surface points d(g,f) and the distance vector between the
*              point in the plane and the point on the first surface d(p,f)
*              we have the following equaton system for computing the
*              next step values along the two parameter directions delta_u
*              and delta_v.
*
*
*              | <f_u,g_n>  <f_v,g_n> |   |delta_u|   | <d(g,f),g_n> |
*              |                      | * |       | = |              |
*              | <f_u,p_n>  <f_v,p_n> |   |delta_v|   | <d(p,f),p_n> |
*
*   
*********************************************************************
*/
{
  double tdum;
  double A[4];   // Equation system matrix
  double b[2];   // Equation system right hand side

  static Point sdiff(3);

          //  First row

  A[0] = fpnt[1]*gpnt[6];
  A[1] = fpnt[2]*gpnt[6];

  sdiff = gpnt[0];
  sdiff -= fpnt[0];
  b[0] = sdiff*gpnt[6];

  // Row scaling
  tdum = std::max(fabs(A[0]),fabs(A[1]));
  tdum = std::max(fabs(tdum),fabs(b[0]));
  if (tdum != DZERO) {
    A[0] /= tdum;
    A[1] /= tdum;
    b[0] /= tdum;
  }


            //  Second row

  A[2]=fpnt[1]*snorm;
  A[3]=fpnt[2]*snorm; 

  sdiff = spoint;
  sdiff -= fpnt[0];
  b[1] = sdiff*snorm;

  // Row scaling
  tdum = std::max(fabs(A[2]),fabs(A[3]));
  tdum = std::max(fabs(tdum),fabs(b[1]));
  if (tdum != DZERO) {
    A[2] /= tdum;
    A[3] /= tdum;
    b[1] /= tdum;
  }

  //    Calculate determinant of equation system
  double det = A[0]*A[3] - A[1]*A[2];

  tdum = std::max(fabs(A[0]),fabs(A[1]));
  tdum = std::max(fabs(A[2]),tdum);
  tdum = std::max(fabs(A[3]),tdum);

  //      if (DEQUAL((tdum+det),tdum)) det =DZERO;
  if (fabs((tdum+det)-tdum) <= REL_PAR_RES*std::max((tdum+det),tdum))
    det  = DZERO;
      
  // If det = 0.0, then the equation system is singular, 
  //	 iteration not possible
  if (det == DZERO) {
    kstat = 1;
    delta[0] = 0.0;
    delta[1] = 0.0;
  }
  else {
    // Using Cramer's rule to find the solution of the system
    kstat = 0;
    delta[0] = (b[0]*A[3] - b[1]*A[1])/det;
    delta[1] = (A[0]*b[1] - A[2]*b[0])/det;
  }

}

} // End of anonymous namespace for helper functions definition.



