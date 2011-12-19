//===========================================================================
//                                                                           
// File: HahnsSurfaceGen.C                                                      
//                                                                           
// Created: Wed Dec 12 12:48:04 2001                                         
//                                                                           
// Author: Sverre Briseid <Sverre.Briseid@math.sintef.no>
//                                                                           
// Revision: $Id: HahnsSurfaceGen.C,v 1.3 2005-11-11 14:06:24 oan Exp $
//                                                                           
// Description: Source files for funcitons in namespace HahnsSurfaceGen.
//              In addiiton several utility functions are included.
//                                                                           
//===========================================================================

#include <fstream> // For debugging
#include "GoTools/creators/HahnsSurfaceGen.h"
#include "GoTools/creators/CoonsPatchGen.h"
#include "GoTools/geometry/GeometryTools.h"
#include "GoTools/geometry/SplineCurve.h"
#include "GoTools/utils/Point.h"
#include "GoTools/geometry/SplineInterpolator.h"
//#include "newmat.h"
#include "GoTools/utils/LUDecomp.h"

using namespace Go;
using std::vector;
using std::max;
using std::min;
using std::swap;

// The following functions are utility functions included in last part of this file.
//namespace /*local_file*/ {


//}


namespace
{

    /// Hermite interpolation of position and icond-3 derivatives
    /// in one endpoint and position and derivative in the other
    /// endpoint, represented as a Bezier curve on the interval [0,1].
    void s9hermit(std::vector<double>& econd, int icond, int idim,
		  double int_length = 1.0);

    /// Compute factors of expression to find derivatives of 
    /// patches in the midpoint of the vertex region.
    void s9coef(Go::Vector3D& evec1, Go::Vector3D& evec2, Go::Vector3D& evec3,
		Go::Vector3D& evec4, int idim, double *cro00,double *cro01,
		double *cro10, double *cro11);

    /// Compute that cross tangent curve belonging to one of the 
    /// blending patches, that is a mapping of the derivative
    /// curves of the previous patch.
    void s9chcoor(std::vector<double>& eblend, int iordblend, 
		  std::vector<double>& epos, int iordpos, 
		  std::vector<double>& ecrt1, int iordcrt1, 
		  double length, Go::Vector3D& evec1, 
		  Go::Vector3D& evec2, Go::Vector3D& evec3,
		  int idim, std::vector<double>& ecrt2,int *jordcrt2);
    /// Compute the product of two Bezier curves of order 4,
    /// one of which is a blending function of dimension 1.
    void s9mult(std::vector<double>& eblend, std::vector<double>& ecoef, int iord,
		int idim, std::vector<double>& ecoefnew);

    /// Compute the angle (in radians) between two vectors
    double s9ang(Go::Vector3D& evec1, Go::Vector3D& evec2, int idim);

    /// Compute given 2. derivative of a blending surface at the 
    /// midpoint of the vertex region when the 2. derivatives of the 
    /// first blending patch are known.
    void s9comder(int ider1, int ider2, std::vector<double>::iterator ederprev, int idim,
		  double aro00, double aro01, double aro10,
		  double aro11, std::vector<double>::iterator eder);


    /// Utility function for computing midpoint of 3-sided region.
    void
    midpoint3(std::vector<shared_ptr<Go::SplineCurve> >& curves, int icurv,
	      std::vector<double>& etwist, std::vector<double>& etang, std::vector<double>& eder);

    /// Utility function for computing midpoint of 4-sided region.
    /// As function is to generate a Coons patch inside, we include two parameters.
    void
    midpoint4(std::vector<shared_ptr<Go::SplineCurve> >& curves, int icurv,
	      std::vector<double>& etwist, std::vector<double>& etang, std::vector<double>& eder,
	      double neighbour_tol, double kink_tol);

    /// Utility function for computing midpoint of 5-sided region.
    void
    midpoint5(std::vector<shared_ptr<Go::SplineCurve> >& curves, int icurv,
	      std::vector<double>& etwist, std::vector<double>& etang, std::vector<double>& eder);

    /// Utility function for computing midpoint of 6-sided region.
    void
    midpoint6(std::vector<shared_ptr<Go::SplineCurve> >& curves, int icurv,
	      std::vector<double>& etwist, std::vector<double>& etang, std::vector<double>& eder);

    // Utility functions used for computing midpoint

    /// Compute second order derivatives of first surface in the midpoint.
    void
    midpoint6_s9der2(std::vector<Go::Point>& ebound, Go::Point& epoint,
		     std::vector<Go::Point>& etang, Go::Point& enorm,
		     std::vector<Go::Point>& evec,
		     int icurv, int idim, std::vector<double>& eder2);

    /// Given both endpoints of a curve segment and the tangent in
    /// one of the endpoints, estimate the curvature radius of the 
    /// curve segment in the endpoint where the tangent is given.
    void curvrad(const Go::Point& epnt1, const Go::Point& epnt2,
		 const Go::Point& etang, int idim, double *crad);

    /// Creates the tangent length for interpolating a
    /// circular arc with an almost equi-oscillating Hermit qubic
    double local_s1325(double aradiu,double angle);

    /// Given position, first and second derivative of a curve at
    /// a point, compute curvature vector.
    void curvature(const std::vector<Go::Point>& eder, int idim, Go::Point& ecurv);


    /// Given the barycentric coordinates of a point in a 3-sided
    /// vertex region, evaluate the value of the ideal blending 
    /// surface of the vertex region in this point.
    void gregoryCharrotFunction3(std::vector<shared_ptr<Go::SplineCurve> >& ecurve,
				 const std::vector<double>& etwist, int ider,
				 const std::vector<double>& ebar, std::vector<double>& eval);

    /// Given the barycentric coordinates of a point in a 5-sided
    /// vertex region, evaluate the value of the ideal blending 
    /// surface of the vertex region in this point.
    void gregoryCharrotFunction5(std::vector<shared_ptr<Go::SplineCurve> >& curves,
				 std::vector<double>& etwist, int ider,
				 std::vector<double>& ebar, std::vector<double>& eval);

    // Utility functions used for computing Gregory Charrot function.

    /// Compute value and 1. derivative of product of barycentric coordinates.
    void sh1467_s9fac1(std::vector<double>& ebar, std::vector<double>& ed1bar,
		       std::vector<double>& ed2bar, int i1, int i2, int i3,
		       double *cfac, double *cfac1, double *cfac2);

    /// Compute value and 1. and 2. derivatives of product of barycentric
    /// coordinates.
    void sh1467_s9fac2(std::vector<double>& ebar, std::vector<double>& ed1bar,
		       std::vector<double>& ed2bar, int i1, int i2, int i3,
		       double *cfac, double *cfac1, double *cfac2,
		       double *cfac11, double *cfac12, double *cfac22);

} // end anonymous namespace


//===========================================================================
vector<shared_ptr<ParamSurface> >
HahnsSurfaceGen::constructPolygonialSurface(vector<shared_ptr<ParamCurve> >&
					      bnd_curves,
					      vector<shared_ptr<ParamCurve> >&
					      cross_curves,
					      double neighbour_tol,
					      double kink_tol,
					      double knot_diff_tol)
//===========================================================================
{
    const int nmb_crvs = (int)bnd_curves.size();
    double startpar = 0.0, endpar=0.0;

    ALWAYS_ERROR_IF((nmb_crvs < 3) || (nmb_crvs > 6),
		"Routine handles input of 3-6 curves.");
    // We expect that all edges are assigned a cross curve, possibly a 0 pointer.
    ALWAYS_ERROR_IF(int(cross_curves.size()) != nmb_crvs,
		"Every boundary curve must be given a cross curve!");

    // If # curves equals 4, CoonsPatch seems more appropriate...
    MESSAGE_IF(nmb_crvs == 4,
		  "As number of edges equals 4, CoonsPatch seems more appropriate.");

    vector<shared_ptr<SplineCurve> > curves;
    for (size_t i = 0; i < bnd_curves.size(); ++i) {
	shared_ptr<SplineCurve> cv =
	    dynamic_pointer_cast<SplineCurve, ParamCurve>(bnd_curves[i]);
	ALWAYS_ERROR_IF(cv.get() == 0,
		    "Input curves must be of type SplineCurve.");
	curves.push_back(cv);
	curves.push_back(dynamic_pointer_cast<SplineCurve, ParamCurve>(cross_curves[i]));
    }

    // Make sure all curves live in 3 dimensions.
    for (size_t i = 0; i < curves.size(); ++i)
	if (curves[i].get() != 0)
	    ALWAYS_ERROR_IF(curves[i]->dimension() != 3,
			"Curves must be 3 dimensional!");

    // We make sure corresponding (bnd and cross) curves share spline space.
    vector<shared_ptr<SplineCurve> > dummy_vec(2);
    for (int i = 0; i < nmb_crvs; ++i) {
	if (curves[2*1+1].get() == 0)
	    continue; // Cross curve may be 0, we do not unify space for 1 curve.
	double tmin = 0.5*(curves[2*i]->startparam() + curves[2*i+1]->startparam());
	double tmax = 0.5*(curves[2*i]->endparam() + curves[2*i+1]->endparam());
 	curves[2*i]->setParameterInterval(tmin, tmax);
 	curves[2*i+1]->setParameterInterval(tmin, tmax);
	copy(curves.begin() + 2*i, curves.begin() + 2*(i + 1), dummy_vec.begin());
	const double knot_diff_tol = 1e-05; // Output basis with knotdiff < 1e-05.
	unifyCurveSplineSpace(dummy_vec, knot_diff_tol);

	// We then seize the opportunity to reparam bnd_curve (and cross tangent).
 	const double tangent_ratio = 1/4.0;
 	CoonsPatchGen::reparamBoundaryCurve(dummy_vec, tangent_ratio);
 	curves[2*i] = dummy_vec[0];
 	curves[2*i+1] = dummy_vec[1];
	
	double ta = curves[2*i]->startparam();
	double tb = curves[2*i]->endparam();
	if (tb - ta > endpar - startpar)
	  {
	    startpar = ta;
	    endpar = tb;
	  }
    }

    // Choose the longest parameter interval.
    for (int i=0; i < nmb_crvs; ++i) {
	curves[2*i]->setParameterInterval(0.0, endpar-startpar);
	curves[2*i+1]->setParameterInterval(0.0, endpar-startpar);
    }

//     // For debugging
//     std::ofstream dump("data/input_dump.g2");
//     for (i = 0; i < curves.size(); ++i) {
// 	curves[i]->writeStandardHeader(dump);
// 	curves[i]->write(dump);
//     }
//     // end of debugging

    // We modify input cross curves to make them consistent wrt tangents and twists.
    vector<shared_ptr<SplineCurve> > mod_cross_curves;
    CoonsPatchGen::getCrossTangs(curves, mod_cross_curves,
				   neighbour_tol, kink_tol);

    // We next add missing cross curves.
    vector<shared_ptr<SplineCurve> > bd_curves(nmb_crvs);
    for (int i = 0; i < nmb_crvs; ++i)
	bd_curves[i] = curves[2*i];
    CoonsPatchGen::addMissingCrossCurves(bd_curves, mod_cross_curves);

    // We use Hahns method to compute the resulting nmb_crvs G1-continuous surfaces.
    vector<shared_ptr<ParamSurface> > surfaces =
	constructHahnsSurface(bd_curves, mod_cross_curves,
			       neighbour_tol, kink_tol, knot_diff_tol);

    return surfaces;
}


//===========================================================================
vector<shared_ptr<ParamSurface> >
HahnsSurfaceGen::constructHahnsSurface(vector<shared_ptr<SplineCurve> >&
					 bnd_curves,
					 vector<shared_ptr<SplineCurve> >&
					 mod_cross_curves,
					 double neighbour_tol,
					 double kink_tol,
					 double knot_diff_tol)
//===========================================================================
{
    // * USE        : 3D geometry only.
    // *
    // *-
    // *              s1401 - Rectangular blending.
    // *              newCurve   - Create new curve.
    // *              freeCurve  - Free space occupied by a curve.
    // *              sh1461_s9coef - Compute coefficients used for computing 
    // *                              derivatives in the midpoint.  
    // *              sh1461_s9hermit - Compute vertices of Bezier curve.  
    // *              sh1461_s9chcoor - Perform coordinate change on derivative curve
    // *              sh1461_s9mult   - Multiply two 4-order Bezier curves.  
    // *              sh1461_s9comder - Compute given 2. deriv. of patch in midpoint
    // *			         of region.
    // *              sh1461_s9ang    - The angle between two vectors.  
    // *              

    int icurv = (int)bnd_curves.size();
    ALWAYS_ERROR_IF(icurv < 3,
		"Expecting input of at least 3 curves.");
    ALWAYS_ERROR_IF(int(mod_cross_curves.size()) != icurv,
		    "Number of bnd curves must equal number of cross tangent curves.");


    vector<shared_ptr<ParamSurface> > surfaces(icurv);
    vector<shared_ptr<SplineCurve> > curves(2*icurv);
    for (size_t i = 0; i < bnd_curves.size(); ++i) {
	curves[2*i] = bnd_curves[i];
	curves[2*i+1] = mod_cross_curves[i];
    }


    int kdim = 3;                 /* Dimension of geometry space.    */
    int kder = 2;                 /* Number of derivatives to
				     compute while evaluating curve. */
    int knmb = 0;                 /* Number of derivatives computed
				     in midpoint of patch.           */
    int kkcrt2;                   /* Order of cross tangent curve.   */
    vector<int> lder(4);
    double tpar;                  /* Parameter value at which to 
                                   evaluate curve.                 */
    double tro00,tro01,tro10,tro11;  /* Coefficients used to compute
					derivatives at midpoint of region. */
//     vector<double>::iterator sder; // = SISL_NULL;/* Array of derivatives at
//                                    // midpoint    */
    vector<double> stwist(icurv*kdim); // = SISL_NULL;/* Twist vectors of
// 	                             //corners of region. */
    vector<double> stang(icurv*kdim); // =SISL_NULL;/*Tangent vectors at midpoint.*/
    vector<double> spos(15);              /* Conditions and vertices of position
                                   curve along inner edge of region.   */
    vector<double> scrt1(21);             /* Conditions and vertices of cross tangent
				     curve along inner edge of region.   */
    vector<double> scrt2(12);             /* Conditions and vertices of cross tangent
				     curve along inner edge of region.   */
    vector<double> stpos(10);             /* Knot vector corresponding to spos.  */
    vector<double> stcrt1(14);            /* Knot vector corresponding to scrt1. */
    vector<double> stcrt2(8);             /* Knot vector corresponding to scrt2. */
    vector<double> sblend(4);             /* Vertices of blending function.      */
    vector<double> sdum(6);               /* Value and 1. derivative of position
				     boundary curve.                     */
    vector<double> sdum2(6);              /* Value and 1. derivative of cross
				     tangent boundary curve.             */
    vector<double> stwist2(4*icurv*kdim); // = SISL_NULL; /* Twist vectors of 4-sided
                                      // region.    */
    vector<double>::iterator st,st2;        /* Pointers into knot vectors. Used to
                                   change parametrization of curve.    */
    shared_ptr<SplineCurve> qc;    /* Pointer to curve.                   */
    vector<shared_ptr<SplineCurve> > qbound(8*icurv); // = SISL_NULL;
                                                        /* Boundary conditions of
							// 4-sided regions. */
    for (size_t ki = 0; ki < qbound.size(); ++ki)
	qbound[ki] = shared_ptr<SplineCurve>(new SplineCurve());

  
    for (int ki=0; ki<4; ki++) lder[ki] = 2;
  
    for (int ki=0; ki<=kder; ki++) knmb += ki + 1;
  
//     /* Initiate knot vectors of position curve and cross tangent
//      curves of inner edges.  */

//     for (ki=0; ki<5; ki++)
// 	{
// 	    stpos[ki] = 0.0;
// 	    stpos[5+ki] = 1.0;
// 	}
//     for (ki=0; ki<4; ki++)
// 	{
// 	    stcrt2[ki] = 0.0;
// 	    stcrt2[4+ki] = 1.0;
// 	}
//     for (ki=0; ki<7; ki++)
// 	{
// 	    stcrt1[ki] = 0.0;
// 	    stcrt1[7+ki] = 1.0;
// 	}


//     /* Allocate scratch for local arrays.  */

//     if ((sder = new0array(icurv*kdim*knmb,DOUBLE)) == SISL_NULL) goto err101;
//     if ((stwist = newarray(icurv*kdim,DOUBLE)) == SISL_NULL) goto err101;
//     if ((stwist2 = newarray(4*icurv*kdim,DOUBLE)) == SISL_NULL) goto err101;
//     if ((stang = newarray(icurv*kdim,DOUBLE)) == SISL_NULL) goto err101;
//     if ((qbound = newarray(8*icurv,SplineCurve*)) == SISL_NULL) goto err101;
    vector<double> sder(icurv*kdim*knmb);

    double length = 0.0;
    for (int ki=0; ki<icurv; ki++)
	{
	  /* Compute parameter interval of boundary curve. */
	  length += curves[2*ki]->endparam() - curves[2*ki]->startparam();

	    /* Initiate twist vector of n-sided region.  */

	    tpar = 0.0;
	    vector<Point> t_point(2); // Position and first derivative.
	    curves[2*ki+1]->point(t_point, tpar, 1);
	    // We transfer from t_point to sdum.
	    for (size_t i = 0; i < t_point.size(); ++i)
		for (int j = 0; j < kdim; ++j)
		    sdum[i*kdim + j] = t_point[i][j];

	    //	    memcopy(stwist+ki*kdim,sdum+kdim,kdim,double);
	    copy(sdum.begin() + kdim, sdum.begin() + 2*kdim,
		 stwist.begin() + ki*kdim);
	}
    length /= 2.0*icurv;
 
    /* Initiate those twist vectors of the rectangular blending surfaces
     that coincide with the twist vectors of the n-sided region.  */

    for (int ki=0; ki<kdim; ki++)
	for (int kj=0; kj<icurv; kj++)
	    {
		int kk = (kj + 1) % icurv;
		stwist2[(kj*4+2)*kdim+ki] = stwist[kk*kdim+ki]/4.0;
	    }
  
      
    /* Set up blending function.  */

    sblend[0] = 1.0;
    sblend[1] = sblend[2] = sblend[3] = 0.0;
    s9hermit(sblend,4,1);
  
    /* Set up value and derivative of the 1. blending surface in the
     midpoint of the vertex region.  */

    if (icurv == 3)
	// *              sh1462 - Evaluate midpoint of 3-sided region.
	midpoint3(curves, icurv, stwist, stang, sder);
    else if (icurv == 4)
	// *              sh1463 - Evaluate midpoint of 4-sided region.
	midpoint4(curves, icurv, stwist, stang, sder, neighbour_tol, kink_tol);
    else if (icurv == 5)
	// *              sh1464 - Evaluate midpoint of 5-sided region.
	midpoint5(curves, icurv, stwist, stang, sder);
    else if (icurv == 6)
	// *              sh1465 - Evaluate midpoint of 6-sided region.
	midpoint6(curves, icurv, stwist, stang, sder);

    for (int ki=0; ki<icurv; ki++)
      for (int kj=0; kj<kdim; kj++)
	stang[ki*kdim+kj] /= length;

    for (int kj=0; kj<kdim; kj++)
      {
	sder[1*kdim+kj] /= length;
	sder[2*kdim+kj] /= length;
	sder[3*kdim+kj] /= (length*length);
	sder[4*kdim+kj] /= (length*length);
	sder[5*kdim+kj] /= (length*length);
      }
  
//     std::ofstream midpoint("debugmidpoint.g2");

    Point p01(sder[0],sder[1],sder[2]);
    Point p02(sder[0]+stang[0],sder[1]+stang[1],
		     sder[2]+stang[2]);
//     SplineCurve midder0(p01, p02);
//     midder0.writeStandardHeader(midpoint);
//     midder0.write(midpoint);
    for (int ki=1; ki<icurv; ki++)
	{
// 	  // @@ VSK, DEBUG OUTPUT, 0203
// 	  Point p1(sder[0],sder[1],sder[2]);
// 	  Point p2(sder[0]+stang[ki*3],sder[1]+stang[ki*3+1],
// 		     sder[2]+stang[ki*3+2]);
// 	  SplineCurve midder(p1, p2);
// 	  midder.writeStandardHeader(midpoint);
// 	  midder.write(midpoint);
// 	  // END DEBUG
	  
	  
	    /* Compute value and derivatives of the other blending surfaces
	       in the midpoint of the vertex region.  */

	    /* Copy the value in the midpoint.  */

	    //	    memcopy(sder+ki*knmb*kdim,sder,kdim,DOUBLE);
	    copy(sder.begin(), sder.begin() + kdim, sder.begin() + ki*knmb*kdim);

	    /* Copy the first derivatives of the current patch.  */

	    //	    memcopy(sder+(ki*knmb+1)*kdim,stang+(ki-1)*kdim,kdim,DOUBLE);
	    copy(stang.begin() + (ki-1)*kdim, stang.begin() + ki*kdim,
		 sder.begin() + (ki*knmb+1)*kdim);
	    //	    memcopy(sder+(ki*knmb+2)*kdim,stang+ki*kdim,kdim,DOUBLE);
	    copy(stang.begin() + ki*kdim, stang.begin() + (ki+1)*kdim,
		 sder.begin() + (ki*knmb+2)*kdim);
      
	    /* The second derivative in the first parameter direction 
	       is equal to the first derivative in the second parameter
	       direction of the previous surface.    */

	    //  memcopy(sder+(ki*knmb+3)*kdim,sder+((ki-1)*knmb+5)*kdim,kdim,DOUBLE);
      	    copy(sder.begin() + ((ki-1)*knmb + 5)*kdim,
		 sder.begin() + ((ki-1)*knmb + 6)*kdim,
		 sder.begin() + (ki*knmb+3)*kdim);

	    /* Compute help variables used to define the rest of the derivatives. */
	    // s9coef(stang + (icurv-1)*kdim, stang, stang+(ki-1)*kdim,
	    //   stang+ki*kdim,kdim, &tro00, &tro01, &tro10, &tro11);
	    Vector3D evec1(&(*(stang.begin() + (icurv-1)*kdim)));
	    Vector3D evec2(&(*(stang.begin())));
	    Vector3D evec3(&(*(stang.begin()+(ki-1)*kdim)));
	    Vector3D evec4(&(*(stang.begin()+ki*kdim)));
	    s9coef(evec1, evec2, evec3, evec4, kdim, &tro00, &tro01, &tro10, &tro11);
      
	    /* Compute mixed derivative of the patch.  */

	    s9comder(1,1,sder.begin()+3*kdim,kdim,tro00,tro01,tro10,tro11,
		     sder.begin()+(ki*knmb+4)*kdim);
	  
	    /* Compute the second derivative in the second parameter direction. */

	    s9comder(0, 2, sder.begin()+3*kdim, kdim, tro00, tro01, tro10, tro11,
		     sder.begin()+(ki*knmb+5)*kdim);
	}  


    /* Set up conditions for 4-edged blending.  */
    for (int ki=0; ki<icurv; ki++)
	{
	    /* Copy value and derivatives of the current patch of the midpoint of
	       the vertex region into arrays containing interpolation conditions 
	       and twist vector of current blending surface. */

	    //	    memcopy(spos,sder+ki*knmb*kdim,kdim,DOUBLE);
	    copy(sder.begin() + ki*knmb*kdim, sder.begin() + (ki*knmb + 1)*kdim,
		 spos.begin());
	    //	    memcopy(spos+kdim,sder+(ki*knmb+2)*kdim,kdim,DOUBLE);
	    copy(sder.begin() + (ki*knmb+2)*kdim, sder.begin() + (ki*knmb+3)*kdim,
		 spos.begin() + kdim);
	    //	    memcopy(spos+2*kdim,sder+(ki*knmb+5)*kdim,kdim,DOUBLE);
	    copy(sder.begin() + (ki*knmb+5)*kdim, sder.begin() + (ki*knmb+6)*kdim,
		 spos.begin() + 2*kdim);
	    //	    memcopy(scrt2,sder+(ki*knmb+1)*kdim,kdim,DOUBLE);
	    copy(sder.begin() + (ki*knmb+1)*kdim, sder.begin() + (ki*knmb+2)*kdim,
		 scrt2.begin());
	    //	    memcopy(scrt2+kdim,sder+(ki*knmb+4)*kdim,kdim,DOUBLE);
	    copy(sder.begin() + (ki*knmb+4)*kdim, sder.begin() + (ki*knmb+5)*kdim,
		 scrt2.begin() + kdim);
	    //	    memcopy(stwist2+4*ki*kdim,sder+(ki*knmb+4)*kdim,kdim,DOUBLE);
      	    copy(sder.begin() + (ki*knmb+4)*kdim, sder.begin() + (ki*knmb+5)*kdim,
		 stwist2.begin() + 4*ki*kdim);

	    /* Compute value and derivatives of the current inner edge curve at 
	       the boundary of the vertex region.  */

	    int kj = (ki < icurv-1) ? ki+1 : 0;
      
	    qc = curves[2*kj];
	    tpar = (*(qc->basis().begin() + qc->order() - 1) +
		    *(qc->basis().begin() + qc->numCoefs()))/2.0;

	    vector<Point> t_point(2); // Position and first derivative.
	    qc->point(t_point, tpar, 1);
	    // We transfer from t_point to sdum.
	    for (size_t i = 0; i < t_point.size(); ++i)
		for (int j = 0; j < kdim; ++j)
		    sdum[i*kdim + j] = t_point[i][j];

	    curves[2*kj+1]->point(t_point, tpar, 1);
	    // We transfer from t_point to sdum.
	    for (size_t i = 0; i < t_point.size(); ++i)
		for (int j = 0; j < kdim; ++j)
		    sdum2[i*kdim + j] = t_point[i][j];
      
	    /* Copy value and derivatives at the current edge into arrays containing 
	       interpolation conditions and twist vector of current blending surface.
	    */

	    //	    memcopy(spos+4*kdim,sdum,kdim,DOUBLE);
	    copy(sdum.begin(), sdum.begin() + kdim, spos.begin() + 4*kdim);

	    for (int kh=0; kh<kdim; kh++)
		{
//  		    spos[3*kdim+kh] = sdum2[kh]/2.0;
//  		    scrt2[3*kdim+kh] = -sdum[kdim+kh]/2.0;
//  		    stwist2[(ki*4+3)*kdim+kh] =
//  			scrt2[2*kdim+kh] = -sdum2[kdim+kh]/4.0;
//  		    stwist2[(kj*4+1)*kdim+kh] = -stwist2[(ki*4+3)*kdim+kh];
		    spos[3*kdim+kh] = sdum2[kh];
		    scrt2[3*kdim+kh] = -sdum[kdim+kh];
		    stwist2[(ki*4+3)*kdim+kh] =
			scrt2[2*kdim+kh] = -sdum2[kdim+kh];
		    stwist2[(kj*4+1)*kdim+kh] = -stwist2[(ki*4+3)*kdim+kh];
		}

	    /* Initiate knot vectors of position curve and cross tangent
	       curves of inner edges.  */
	    //	    length = 1.0; // for debugging
	    for (int kh=0; kh<5; kh++)
		{
		    stpos[kh] = 0.0;
		    //		    stpos[5+kh] = 1.0;
		    stpos[5+kh] = length;
		}
	    for (int kh=0; kh<4; kh++)
		{
		    stcrt2[kh] = 0.0;
		    //		    stcrt2[4+kh] = 1.0;
		    stcrt2[4+kh] = length;
		}
	    for (int kh=0; kh<7; kh++)
		{
		    stcrt1[kh] = 0.0;
		    //		    stcrt1[7+kh] = 1.0;
		    stcrt1[7+kh] = length;
		}

      
	    /* Perform Hermit interpolation of position and first cross tangent
	       curve. */
	    s9hermit(spos,5,kdim, length);
      
	    s9hermit(scrt2,4,kdim, length);
      
	    /* Construct second cross tangent curve. */
	    Vector3D evec1(&(*(sder.begin()+(knmb*ki+1)*kdim)));
	    Vector3D evec2(&(*(sder.begin()+(knmb*ki+2)*kdim)));
	    Vector3D evec3(&(*(sder.begin()+(kj*knmb+2)*kdim)));
	    s9chcoor(sblend, 4, spos, 5, scrt2, 4, length, evec1, evec2, 
		     evec3, kdim, scrt1, &kkcrt2);

      
	    /* Represent inner boundary conditions as curves. */

 	    qbound[8*kj] =
		shared_ptr<SplineCurve>(new SplineCurve(5, 5, stpos.begin(),
							    spos.begin(), kdim));
	    qbound[8*kj+1] =
		shared_ptr<SplineCurve>(new SplineCurve(kkcrt2, kkcrt2,
							    stcrt1.begin(),
							    scrt1.begin(), kdim));
	    qbound[8*ki+6] =
		shared_ptr<SplineCurve>(new SplineCurve(5, 5, stpos.begin(),
							    spos.begin(), kdim));
	    qbound[8*ki+7] =
		shared_ptr<SplineCurve>(new SplineCurve(4, 4, stcrt2.begin(),
							    scrt2.begin(), kdim));
      
	    /* Split current boundary of the vertex region.  */
	    qbound[8*ki+4] = 
#if ((_MSC_VER > 0) && (_MSC_VER < 1300))
		shared_ptr<SplineCurve>(static_cast<SplineCurve*>(qc->subCurve(qc->startparam(), tpar,
										   knot_diff_tol)));
#else
		shared_ptr<SplineCurve>(qc->subCurve(qc->startparam(), tpar,
						       knot_diff_tol));
#endif


	    qbound[8*kj+2] =
#if ((_MSC_VER > 0) && (_MSC_VER < 1300))
		shared_ptr<SplineCurve>(static_cast<SplineCurve*>(qc->subCurve(tpar, qc->endparam(),
										   knot_diff_tol)));
#else
		shared_ptr<SplineCurve>(qc->subCurve(tpar, qc->endparam(),
						       knot_diff_tol));
#endif

	    //	    s1710(curves[2*kj+1],tpar,qbound+8*ki+5,qbound+8*kj+3,&kstat);


	    qbound[8*ki+5] =
#if ((_MSC_VER > 0) && (_MSC_VER < 1300))
		shared_ptr<SplineCurve>(static_cast<SplineCurve*>(curves[2*kj+1]->
					  subCurve(curves[2*kj+1]->startparam(),
						   tpar, knot_diff_tol)));
#else
		shared_ptr<SplineCurve>(curves[2*kj+1]->
					  subCurve(curves[2*kj+1]->startparam(),
						   tpar, knot_diff_tol));
#endif


	    qbound[8*kj+3] =
#if ((_MSC_VER > 0) && (_MSC_VER < 1300))
		shared_ptr<SplineCurve>(static_cast<SplineCurve*>(curves[2*kj+1]->
								      subCurve(tpar,
									       curves[2*kj+1]->endparam(),
									       knot_diff_tol)));
#else
		shared_ptr<SplineCurve>(curves[2*kj+1]->
					  subCurve(tpar,
						   curves[2*kj+1]->endparam(),
						   knot_diff_tol));
#endif


	    // VSK, 0203. Start from zero.
	    double endpar = qbound[8*ki+4]->endparam();
  	    qbound[8*ki+4]->setParameterInterval(0.0, endpar);
  	    qbound[8*kj+2]->setParameterInterval(0.0, endpar);
  	    qbound[8*ki+5]->setParameterInterval(0.0, endpar);
  	    qbound[8*kj+3]->setParameterInterval(0.0, endpar);

// 	    // For debugging
// 	    std::ofstream dump3("data/bnd_part_debugdump.g2");
// 	    qbound[8*kj]->writeStandardHeader(dump3);
// 	    qbound[8*kj]->write(dump3);
// 	    qbound[8*ki+6]->writeStandardHeader(dump3);
// 	    qbound[8*ki+6]->write(dump3);
// 	    qbound[8*ki+4]->writeStandardHeader(dump3);
// 	    qbound[8*ki+4]->write(dump3);
// 	    qbound[8*kj+2]->writeStandardHeader(dump3);
// 	    qbound[8*kj+2]->write(dump3);
// 	    // end of debugging

	    /* Turn direction of curves corresponding to standard edge 3.  */
	    //	    s1706(qbound[8*ki+4]);
	    qbound[8*ki+4]->reverseParameterDirection();
	    //	    s1706(qbound[8*ki+5]);
	    qbound[8*ki+5]->reverseParameterDirection();

	    /* Adjust length of derivative curves. */

//  	    for (kk=0; kk<kdim*(qbound[8*ki+5]->numCoefs()); kk++)
//  		*(qbound[8*ki+5]->coefs_begin()+kk) *= 0.5;
      
//  	    for (kk=0; kk<kdim*(qbound[8*kj+3]->numCoefs()); kk++)
//  		*(qbound[8*kj+3]->coefs_begin()+kk) *= 0.5;

// 	    /* Reparametrize bound curves in order to get correct tangent length.*/

// 	    for (kh=2; kh<6; kh++)
// 		{
// 		    kk = (kh > 3) ? 8*ki : 8*kj;
// 		    kk += kh;

// 		    // We may not alter knots as we whish. Extract knots.
// 		    vector<double> knots(qbound[kk]->basis().begin(),
// 					 qbound[kk]->basis().end());
// 		    for (st = knots.begin(), tstart = st[0],
// 			     st2 = st + qbound[kk]->numCoefs() + qbound[kk]->order();
// 			 st < st2; st++)
// 			*st = 2.0*(*st - tstart);
// 		    // We're done altering the knots, time to alter curve.
// 		    qbound[kk] = shared_ptr<SplineCurve>
// 			(new SplineCurve(qbound[kk]->numCoefs(),
// 					   qbound[kk]->order(), knots.begin(),
// 					   qbound[kk]->coefs_begin(),
// 					   qbound[kk]->dimension()));
// 		}
      
	}
  
    for (int ki=0; ki<icurv; ki++)
	{
	    /* Perform rectangular blending.  */

	    //	    s1401(qbound+8*ki,stwist2+4*ki*kdim,surfaces+ki,&kstat);  
	    // Go routine assumes bnd curves and cross curves are separated:
	    vector<shared_ptr<SplineCurve> > bnd_crvs, cross_crvs;
	    for (int kj = 0; kj < 4; ++kj) {
		bnd_crvs.push_back(qbound[8*ki+2*kj]);
		cross_crvs.push_back(qbound[8*ki+2*kj+1]);
	    }

// 	    // For debugging
// 	    std::ofstream dump("data/bnd_debugdump.g2");
// 	    std::ofstream dump2("data/cross_debugdump.g2");
// 	    for (kj = 0; kj < 4; ++kj) {
// 		bnd_crvs[kj]->writeStandardHeader(dump);
// 		bnd_crvs[kj]->write(dump);
// 		cross_crvs[kj]->writeStandardHeader(dump2);
// 		cross_crvs[kj]->write(dump2);
// 	    }
// 	    // end of debugging

	    // We need to unify spine space of corresponding curves.
	    vector<shared_ptr<SplineCurve> > dummy_vector_u(4), dummy_vector_v(4);
	    for (int kj = 0; kj < 2; ++kj) {
		dummy_vector_u[kj*2] = bnd_crvs[kj*2];
		dummy_vector_u[kj*2+1] = cross_crvs[kj*2];
		dummy_vector_v[kj*2] = bnd_crvs[kj*2+1];
		dummy_vector_v[kj*2+1] = cross_crvs[kj*2+1];
	    }
	    unifyCurveSplineSpace(dummy_vector_u, knot_diff_tol);
	    unifyCurveSplineSpace(dummy_vector_v, knot_diff_tol);
	    // As objects may have changed, we must extract the curves.
	    for (int kj = 0; kj < 2; ++kj) {
		bnd_crvs[kj*2] = dummy_vector_u[kj*2];
		cross_crvs[kj*2] = dummy_vector_u[kj*2+1];
		bnd_crvs[kj*2+1] = dummy_vector_v[kj*2];
		cross_crvs[kj*2+1] = dummy_vector_v[kj*2+1];
	    }

	    // As Go routine assumes curves form a counter clockwise loop, we turn:
	    bnd_crvs[2]->reverseParameterDirection();
	    cross_crvs[2]->reverseParameterDirection();
	    bnd_crvs[3]->reverseParameterDirection();
	    cross_crvs[3]->reverseParameterDirection();
	    // As our curves have been preprocessed, we do no not need to alter
	    // the curves; hence this version of CoonsPatch
#if ((_MSC_VER > 0) && (_MSC_VER < 1300))
	    surfaces[ki] = shared_ptr<ParamSurface>
		(static_cast<ParamSurface*>(CoonsPatchGen::createCoonsPatch
					      (bnd_crvs, cross_crvs)));
#else
	    surfaces[ki] = shared_ptr<ParamSurface>
		(CoonsPatchGen::createCoonsPatch(bnd_crvs, cross_crvs));
#endif

	    /* s1390(qbound+8*ki,surfaces+ki,lder,&kstat); */
	}
  
    /* Blending performed.  */

//     // For debugging
//     std::ofstream dump4("data/surfaces_debugdump.g2");
//     for (ki = 0; ki < surfaces.size(); ++ki) {
// 	surfaces[ki]->writeStandardHeader(dump4);
// 	surfaces[ki]->write(dump4);
//     }
//     // end of debugging

  
    return surfaces;
}


namespace
{


//static void
// sh1461_s9hermit(double econd[],int icond,int idim,int *jstat)
void s9hermit(vector<double>& econd, int icond, int idim, double int_length)
// /*
// *********************************************************************
// *                                                                   
// * PURPOSE    : Hermite interpolation of position and icond-3 derivatives
// *              in one endpoint and position and derivative in the other
// *              endpoint, represented as a Bezier curve on the interval [0,1].
// *
// *
// *
// * INPUT      : icond      - Number of interpolation conditions. 
// *                           icond = 4 or icond = 5.
// *              idim       - Dimension of geometry space.
// *
// *
// * INPUT/OUTPUT : econd    - Interpolation conditions as input, Bezier coefficients
// *                           as output. The dimension is icond*idim.
// *                       
// *
// * OUTPUT     : jstat      - status messages  
// *                                         > 0      : warning
// *                                         = 0      : ok
// *                                         < 0      : error
// *
// *
// *********************************************************************
// */
{ 
    int ki;    /* Index.  */

    /* Test input. The number of conditions has to be 4 or 5.  */

    ALWAYS_ERROR_IF((icond != 4) && (icond != 5),
		"Wrong input.");
  
    if (icond == 4)
	{
	    /* Hermit interpolation with Bezier curve of order 4.  */

	    for (ki=0; ki<idim; ki++)
		{
		    econd[idim+ki] = econd[idim+ki]*int_length/3.0 + econd[ki];
		    econd[2*idim+ki] =
			econd[2*idim+ki]*int_length/3.0 + econd[3*idim+ki];
		}
	}
  
    if (icond == 5)
	{
	    /* Hermit interpolation with Bezier curve of order 5.  */

	    for (ki=0; ki<idim; ki++)
		{
		    econd[idim+ki] = econd[idim+ki]*int_length/4.0 + econd[ki];
		    econd[2*idim+ki] = econd[2*idim+ki]*int_length*int_length/12.0 
			+ 2.0*econd[idim+ki] - econd[ki];
		    econd[3*idim+ki] =
			econd[3*idim+ki]*int_length/4.0 + econd[4*idim+ki];	  
		}
	}
   
	return;
}


// static void
//   sh1461_s9coef(double evec1[],double evec2[],double evec3[],
// 		double evec4[],int idim,double *cro00,double *cro01,
// 		double *cro10,double *cro11,int *jstat)
void s9coef(Vector3D& evec1, Vector3D& evec2,
	    Vector3D& evec3, Vector3D& evec4,
	    int idim, double *cro00,double *cro01,
	    double *cro10, double *cro11)
/*
*********************************************************************
*                                                                   
* PURPOSE    : Compute factors of expression to find derivatives of 
*              patches in the midpoint of the vertex region.
*
*
*
* INPUT      : evec1   - Derivative in first parameter direction of first patch.
*              evec2   - Derivative in second parameter direction of first patch.
*              evec3   - Derivative in first parameter direction of current patch.
*              evec4   - Derivative in second parameter direction of current patch.
*                        NB! All these derivative vectors are expected to lie in
*                            the same plane.
*              idim    - Dimension of geometry space.
*
*
* OUTPUT     : cro00   - First factor.
*              cro01   - Second factor.
*              cro10   - Third factor.
*              cro11   - Fourth factor.
*              jstat   - status messages  
*                                         > 0      : warning
*                                         = 0      : ok
*                                         < 0      : error
*
*
* CALLS      : s6length  - Length of vector.
*              s6scpr    - Scalar product between two vectors.
*              s6crss    - Cross product between two vectors.
*
*********************************************************************
*/
{
    int ksin,ksin1,ksin2,ksin3,ksin4;    /* Sign of sinus of angles.      */
    double tang,tang1,tang2,tang3,tang4; /* Angles between input vectors. */
    double tl1,tl2,tl3,tl4;              /* Lengths of input vectors.     */
    double tsin,tsin1,tsin2,tsin3,tsin4; /* Sinus of angles.              */
    Vector3D snorm;                     /* Normal of plane spanned by
					    tangent vectors.              */
    Vector3D svec;                      /* Vector in tangent plane normal
					    to a specific tangent vector. */
    double tlvec;                        /* Length of the vector svec.    */
  
    /* Compute the normal to the tangent plane spanned by the input vectors.  */

    //    s6crss(evec1,evec2,snorm);
    snorm = evec1.cross(evec2);  

    /* Compute the lengths of the vectors evec1 to evec4. */
    tl1 = evec1.length();
    tl2 = evec2.length();
    tl3 = evec3.length();
    tl4 = evec4.length();
  
    /* Compute the vector lying in the tangent plane normal to evec1. */
    svec = snorm.cross(evec1);
    tlvec = svec.length();

    /* Compute the sinus of angles between evec1 and the other input vectors. */
    tsin = svec*evec2 / (tlvec*tl2);
    tsin1 = svec*evec3 / (tlvec*tl3);
    tsin2 = svec*evec4 / (tlvec*tl4);
  
    /* Compute the vector lying in the tangent plane normal to evec2. */
    svec = snorm.cross(evec2);
    tlvec = svec.length();

    /* Compute the sinus of angles between evec2 and the later input vectors. */
    tsin3 = svec*evec3 / (tlvec*tl3);
    tsin4 = svec*evec4 / (tlvec*tl4);

    /* Fetch the sign of the sinuses of the angles.  */
    ksin = (tsin < 0) ? -1 : 1;
    ksin1 = (tsin1 < 0) ? -1 : 1;
    ksin2 = (tsin2 < 0) ? -1 : 1;
    ksin3 = (tsin3 < 0) ? -1 : 1;
    ksin4 = (tsin4 < 0) ? -1 : 1;
  
    /* Compute the angles in a more stable way. */
    tang = s9ang(evec1,evec2,idim);
    tang1 = s9ang(evec1,evec3,idim);
    tang2 = s9ang(evec1,evec4,idim);
    tang3 = s9ang(evec2,evec3,idim);
    tang4 = s9ang(evec2,evec4,idim);
  
    /* Compute the sinuses of the angles.  */
    tsin = ksin*sin(tang);
    tsin1 = ksin1*sin(tang1);
    tsin2 = ksin2*sin(tang2);
    tsin3 = ksin3*sin(tang3);
    tsin4 = ksin4*sin(tang4);
  
    /* Compute the output factors.  */
    *cro00 = (tl3*tsin1)/(tl2*tsin);
    *cro01 = (tl4*tsin2)/(tl2*tsin);
    *cro10 = -(tl3*tsin3)/(tl1*tsin);
    *cro11 = -(tl4*tsin4)/(tl1*tsin);
    
    return;
}


// static void
//   sh1461_s9chcoor(double eblend[],int iordblend,double epos[],
// 		  int iordpos,double ecrt1[],int iordcrt1,
// 		  double evec1[],double evec2[],double evec3[],
// 		  int idim,double ecrt2[],int *jordcrt2,int *jstat)
void s9chcoor(vector<double>& eblend, int iordblend,
				 vector<double>& epos, int iordpos,
				 vector<double>& ecrt1, int iordcrt1,
				 double length,
				 Vector3D& evec1, Vector3D& evec2,
				 Vector3D& evec3, int idim,
				 vector<double>& ecrt2,int *jordcrt2)
/*
*********************************************************************
*                                                                   
* PURPOSE    : Compute that cross tangent curve belonging to one of the 
*              blending patches, that is a mapping of the derivative
*              curves of the previous patch.
*
*
*
* INPUT      : eblend    - Vertices of blending function. The dimension
*                          of the array is iordblend.
*              iordblend - Number of vertices of blending function. 
*                          iordblend = 4.
*              epos      - Vertices of position curve. This curve is
*                          between the corners (0,0) and (0,1) of the
*                          previous patch, and the corners (0,0) and 
*                          (1,0) of the current patch.
*              iordpos   - Number of vertices of the position curve.
*                          iordpos = 5.
*              ecrt1     - Vertices of input cross tangent curve. This
*                          curve is the derivative curve in the 1. parameter
*                          direction along the edge from (0,0) to (0,1)
*                          of the previous patch.
*              iordcrt1  - Number of vertices of input cross tangent curve.
*                          iordcrt1 = 4.
*              length    - Length of the parameter interval of the curves.
*              evec1     - Tangent vector at the corner (0,0) in the 1. 
*                          parameter direction of the previous patch.
*              evec2     - Tangent vector at the corner (0,0) in the 2. 
*                          parameter direction of the previous patch, and
*                          in the 1. parameter direction of the current patch.
*              evec3     - Tangent vector at the corner (0,0) in the 2. 
*                          parameter direction of the current patch.
*              idim      - Dimension of geometry space.
*
*
* OUTPUT     : ecrt2     - Vertices of the produced cross tangent curve. 
*                          This curve is the derivative curve in the 2.
*                          parameter direction along the edge from (0,0) to
*                          (1,0) of the current patch.
*              jordcrt2  - Number of vertices of produced cross tangent curve.
*                          *jordcrt2 = 7.
*              jstat     - status messages  
*                                         > 0      : warning
*                                         = 0      : ok
*                                         < 0      : error
*
*
* METHOD     : Let pc1 be the input cross tangent curve, pc2 be the
*              derivative curve along the position curve and alpha the 
*              blending function. Then the output cross tangent curve is 
*              given by :
*              pc3(s) = ((my+1)*alpha(s)-1)*pc1(s) + lambda*alpha(s)*pc2(s)
*              my and lambda are constants depending on the lengths of 
*              input vectors, evec1, evec2 and evec3, and the angles 
*              between them. See the paper of Hahn.
*
* CALLS     : s6length - Length of vector.
*
*********************************************************************
*/
{
    double tl1,tl2,tl3;              /* Lengths of the input vectors. */
    double tang1,tang2;              /* Angles between input vectors. */
    double tsin1,tsin2,tsin3;        /* Sinus of angles between input
					vectors.                      */
    double tmy,tmy1;                 /* Coefficients dependant on the
					input vectors.                */
    double tlambda;                  /* Coefficient dependant on the
					input vectors.                */
    vector<double> scoef(12);                /* Vertices of the derivative
					curve along the position curve. */
    vector<double> sc1(21), sc2(21); /* Coefficients of products of curves. */
    vector<double> sblend2(4);               /* Coefficients of curve equal to
					a factor times the blending 
					curve minus 1.                 */
  
    //    if (iordblend != 4 || iordcrt1 != 4) goto err002;   
    ALWAYS_ERROR_IF((iordblend != 4) || (iordcrt1 != 4),
		"Wrong input!");

    *jordcrt2 = iordblend + iordcrt1 - 1;
  
    /* Compute constants.  */
    tl1 = evec1.length();
    tl2 = evec2.length();
    tl3 = evec3.length();

    tang1 = s9ang(evec1, evec2, idim);
    tang2 = s9ang(evec2, evec3, idim);
    tsin1 = sin(tang1);
    tsin2 = sin(tang2);
    tsin3 = sin(tang1+tang2);
    tmy = - (tl3*tsin2)/(tl1*tsin1);
    tlambda = (tl3*tsin3)/(tl2*tsin1);
    tmy1 = tmy + 1.0;
  
    /* Compute coefficients of the curve (my+1)alpha(s)-1.  */
    for (int ki=0; ki<4; ki++) sblend2[ki] = tmy1*eblend[ki] - 1.0;
  
    /* Compute coefficients of derivative curve along position curve. */
    for (int ki=0; ki<4*idim; ki++)
	scoef[ki] = 4.0*(epos[idim+ki] - epos[ki])/length;

    /* Compute first part of the second cross tangent curve. */
    s9mult(eblend,scoef,4,idim,sc1);

    /* Compute second part of the second cross tangent curve. */
    s9mult(sblend2,ecrt1,4,idim,sc2);
  
    /* Compute cross tangent curve.  */
    for (int ki=0; ki<7*idim; ki++)
	ecrt2[ki] = (tlambda*sc1[ki] + sc2[ki]);
  
	return;
}



// static void
//   sh1461_s9mult(double eblend[],double ecoef[],int iord,
// 		int idim,double ecoefnew[],int *jstat)
void s9mult(vector<double>& eblend, vector<double>& ecoef,
			       int iord, int idim, vector<double>& ecoefnew)
/*
*********************************************************************
*                                                                   
* PURPOSE    : Compute the product of two Bezier curves of order 4,
*              one of which is a blending function of dimension 1.
*
*
*
* INPUT      : eblend    - Vertices of blending function. The dimension
*                          of the array is iord, i.e. 4.
*              ecoef     - Vertices of the other Bezier curve. The 
*                          dimension of the array is iord*idim, i.e. 4*idim.
*              iord      - Order of the Bezier curves. iord = 4.
*              idim      - Dimension of geometry space.
*
*
* OUTPUT     : ecoefnew  - Vertices of the product curve. The array is
*                          allocated outside this routine, and has
*                          dimension (2*iord-1)*idim, i.e. 7*idim
*              jstat     - status messages  
*                                         > 0      : warning
*                                         = 0      : ok
*                                         < 0      : error
*
* USE        : Used in Hahn's method when computing the cross tangent
*              curve which is a mapping of the derivative curves of the
*              previous patch. Called from s9chcoor.
*
*********************************************************************
*/
{
    int ki;   /* Index.   */
  
    /* Test if the order of the curves is equal to 4.  */

    //    if (iord != 4) goto err001;
    ALWAYS_ERROR_IF(iord != 4,
		"Order of Bezier curve must be 4!");

    for (ki=0; ki<idim; ki++)
	{
	    /* Compute the vertices of the product curve. */

	    ecoefnew[ki] = eblend[0]*ecoef[ki];
	    ecoefnew[idim+ki] = (eblend[0]*ecoef[idim+ki] 
				 + eblend[1]*ecoef[ki])/2.0;
	    ecoefnew[2*idim+ki] = (eblend[0]*ecoef[2*idim+ki]
				   + 3.0*eblend[1]*ecoef[idim+ki]
				   + eblend[2]*ecoef[ki])/5.0;
	    ecoefnew[3*idim+ki] = (eblend[0]*ecoef[3*idim+ki] + eblend[3]*ecoef[ki]
				   + 9.0*(eblend[1]*ecoef[2*idim+ki] + 
						  eblend[2]*ecoef[idim+ki]))/20.0;
	    ecoefnew[4*idim+ki] = (eblend[1]*ecoef[3*idim+ki]
				   + 3.0*eblend[2]*ecoef[2*idim+ki]
				   + eblend[3]*ecoef[idim+ki])/5.0;      
	    ecoefnew[5*idim+ki] = (eblend[2]*ecoef[3*idim+ki] 
				   + eblend[3]*ecoef[2*idim+ki])/2.0;
	    ecoefnew[6*idim+ki] = eblend[3]*ecoef[3*idim+ki];
	}

	return;
}



// static double
//   sh1461_s9ang(double evec1[],double evec2[],int idim)
double s9ang(Vector3D& evec1, Vector3D& evec2, int idim)
/*
*********************************************************************
*
*********************************************************************
*                                                                   
* PURPOSE    : Compute the angle (in radians) between two vectors
*
*
*
* INPUT      : evec1   - First vector 
*              evec2   - Second vector 
*              idim    - Dimension of the space in which the vectors lie.
*
*
*
* OUTPUT     : s9ang   - Angle in radians between vectors
*
*
* METHOD     : Make cosine of the angle by computing the scalar product,
*              then divide by the length of the two vectors.
*
* REFERENCES :
*
*-
* CALLS      : s6scpr   - Scalar product between two vectors.
*              s6length - Length of vector.
*
* WRITTEN BY : Tor Dokken SI, 88-07.
*              Arne Laksaa SI, 89-07.
*              Vibeke Skytt SI, 90-04.
*
*********************************************************************
*/                                     
{
    double tscpr,tang,tlength1,tlength2,tcos;
    tscpr = evec1*evec2;
    tlength1 = evec1.length();
    tlength2 = evec2.length();
  
    if ((tlength1 == 0) || (tlength2 == 0))
	tang = 0;
    else
	{
	    tcos = tscpr/(tlength1*tlength2);
	    tcos = min(1.0, tcos);
	    tcos = max(-1.0, tcos);
	    tang = acos(tcos);
	}
      
    return(tang);
}



// static  void
//   sh1461_s9comder(int ider1,int ider2,double ederprev[],int idim,
// 		  double aro00,double aro01,double aro10,
// 		  double aro11,double eder[],int *jstat)
void s9comder(int ider1, int ider2,
				 vector<double>::iterator ederprev, int idim,
				 double aro00, double aro01, double aro10,
				 double aro11, vector<double>::iterator eder)
/*
*********************************************************************
*                                                                   
* PURPOSE    : Compute given 2. derivative of a blending surface at the 
*              midpoint of the vertex region when the 2. derivatives of the 
*              first blending patch are known.
*
*
*
* INPUT      : ider1    - Order of differentiation in first parameter direction.
*              ider2    - Order of differentiation in second parameter direction.
*              ederprev - 2. derivatives of 1. patch stored in the following 
*                         order, (2,0)-, (1,1)- and (0,2)-derivative.
*              idim     - Dimension of geometry space.
*              aro00    - First factor of coordinate transformation.
*              aro01    - Second factor of coordinate transformation.
*              aro10    - Third factor of coordinate transformation.
*              aro11    - Fourth factor of coordinate transformation.
*
*
* OUTPUT     : eder    - Actual 2. derivative.
*              jstat   - status messages  
*                                         > 0      : warning
*                                         = 0      : ok
*                                         < 0      : error
*
*
*********************************************************************
*/
{
    double num_tol = 1e-16; // @@ Hardcoded tolerance.

    int kj,kk,kh;            /* Counters.  */
    double t00,t01,t10,t11;  /* Coefficients used to compute
				derivatives at midpoint of region. */
    double tfac1,tfac2;      /* Factor used to compute derivatives 
				at midpoint of region.  */
  
    /* Test input.  */

    //    if (ider1 + ider2 != 2) goto err001;
    ALWAYS_ERROR_IF(ider1 + ider2 != 2,
		"Sum of order of partial derivatives order equal 2!");

    /* Compute requested 2. derivative.  */

    t00 = t10 = 1.0;
    for (kj=0; kj<ider1; kj++)
	t00 *= aro00;

    for (kj=0; kj<=ider1; kj++)
	{
	    tfac1 = (ider1 > 0) ? ((kj % ider1) + 1) : 1.0;

	    t01 = t11 = 1.0;
	    for (kk=0; kk<ider2; kk++) t01*=aro01;

	    for (kk=0; kk<=ider2; kk++)
		{
		    tfac2 = (ider2 > 0) ? ((kk % ider2) + 1) : 1.0;

		    for (kh=0; kh<idim; kh++)
			eder[kh] +=
			    tfac1*tfac2*t00*t01*t10*t11*ederprev[(2-kj-kk)*idim+kh];

		    if (fabs(aro01 - num_tol) > num_tol)
			t01 /= aro01;
		    else if (kk == ider2-1)
			t01 = 1.0;
		    t11 *= aro11;
		}
	    if (fabs(aro00 - num_tol) > num_tol)
		t00 /= aro00;
	    else if (kj == ider1-1)
		t00 = 1.0;
	    t10 *= aro10;	  
	}
  
	return;
}


void
gregoryCharrotFunction3(vector<shared_ptr<SplineCurve> >&
					   ecurve,
					   const vector<double>& etwist, int ider,
					   const vector<double>& ebar,
					   vector<double>& eval)
    /*
*********************************************************************
*                                                                   
* PURPOSE    : Given the barycentric coordinates of a point in a 3-sided
*              vertex region, evaluate the value of the ideal blending 
*              surface of the vertex region in this point.
*
*
*
* INPUT      : ecurve - Position and cross-tangent curves around the vertex
*                       region. For each edge of the region position and cross-
*                       tangent curves are given. The curves follow each other
*                       around the region and are oriented counter-clock-wise.
*                       The dimension of the array is 6.
*              etwist - Twist-vectors of the corners of the vertex region. The
*                       first element of the array is the twist in the corner
*                       before the first edge, etc. The dimension of the array
*                       is 3*kdim.
*              ider   - Number of derivatives to compute. Directions of 
*                       differentiation is that of the two first barycentric
*                       coordinates. 0 <= ider <= 2.
*              ebar   - Barycentric coordinates of the point to be evaluated.
*                       The dimension of the array is 3.
*                       
*
* OUTPUT     : eval   - Value and derivatives of ideal blending surface in the 
*                       given point. Dimension of the array is 3*(1+..+(ider+1)).
*              jstat  - status messages  
*                                         > 0      : warning
*                                         = 0      : ok
*                                         < 0      : error
*
*
* METHOD     : A functional description of the ideal surface is given as
*              a blend between three surfaces, each of which fulfill the
*              continuity requirements over two edges.
*
* REFERENCES : Gregory and Charrot : A C1 Triangular Interpolation Patch for
*                                    Computer Aided Geometric Design
*
* USE        : 3D geometry only.
*
*-
* CALLS      : s1221 - Evaluate curve at a given parameter value.
*
* WRITTEN BY : Vibeke Skytt, SI, Dec 89.
*
*********************************************************************
*/
{
    int ki,kj,kk,kh;     /* Counters.                                       */
    int kder;            /* Number of derivatives to evaluate.              */
    int kdim = 3;        /* Dimension of geometry.                          */
    int kwarn = 0;       /* Indicates if a warning is to be sendt.          */
    int knmb;            /* Number of doubles pr derivative.                */
    int kl = 0;          /* Number of derivatives in the output array.      */
    double tpar1;        /* Parameter value of edge between actual corner 
			    and next corner.                                */
    double tpar2;        /* Parameter value of edge between actual corner
			    and previous corner.                            */
    double tlambi;       /* First barycentric coordinate.                   */
    double tlambj;       /* Second barycentric coordinate.                  */
    double tlambk;       /* Third barycentric coordinate.                   */
    vector<double> salpha(3); /* Weight of actual blending surface.              */
    vector<double> sp(9);     /* Value of function interpolating two sides.      */
    vector<double> sstart(3); /* Start-parameters of edge-curves.                */
    vector<double> send(3);    /* End-parameters of edge-curves.                  */
    vector<double> sint(3);    /* Parameter intervals of edge-curves.             */
    vector<double> spos1(27);  /* Position of edge-curve in tpar1.                */
    vector<double> sder1(27);  /* Cross tangent in tpar1.                         */
    vector<double> spos2(27);  /* Position of edge-curve in tpar2.                */
    vector<double> sder2(27);  /* Cross tangent in tpar2.                         */
    vector<double> scorn(9);   /* Position of edge-curves in actual corner of 
			    vertex region.                                  */
    vector<double> scornder1(9);/* Tangent in corner along next edge.              */
    vector<double> scornder2(9);/* Tangent in corner along previous edge.          */

    /* Test input.  */

    if (ider > 2) kwarn = 1;
  
    /* Initialise.  */

    kder = ider;
    knmb = kdim*(ider+1);
    for (ki=0; ki<ider; ki++) kl += ider + 1;
  
    /* Initiate output array to zero.  */

    for (kh=0; kh<kl*kdim; kh++) eval[kh] = 0.0;
        
    /* Get endpoints of parameter intervals of edge curves.  */

    for (ki=0; ki<3; ki++)
	{
	    sstart[ki] = *(ecurve[2*ki]->basis().begin() + ecurve[2*ki]->order()-1);
	    send[ki] = *(ecurve[2*ki]->basis().begin() + ecurve[2*ki]->numCoefs());
	    sint[ki] = send[ki] - sstart[ki];
	}

    /* Evaluate position and cross-tangent curves at points on
       the edges needed when evaluating surface.        */
  
    for (ki=0; ki<3; ki++)
	{
	    kj = (ki+1) % 3;
	    kk = (ki+2) % 3;

	    /* Copy barycentric coordinates to local variables.  */

	    tlambi = ebar[ki];      
	    tlambj = ebar[kj];
	    tlambk = ebar[kk];
      
	    /* Find parameter values of points on edges used to evaluate
	       actual blending surface.  */

	    tpar1 = (1.0 - tlambj)*sstart[ki] + tlambj*send[ki];
	    tpar2 = tlambk*sstart[kk] + (1.0 - tlambk)*send[kk];
 
	    /* Evaluate position and cross-tangent curves at first
	       found parameter value.  */

	    vector<Point> pts(kder+1);
	    ecurve[2*ki]->point(pts, tpar1, kder);
	    for (size_t j = 0; j < pts.size(); ++j)
		copy(pts[j].begin(), pts[j].end(),
		     spos1.begin() + knmb*ki + j*kdim);
      
	    ecurve[2*ki+1]->point(pts, tpar1, kder);
	    for (size_t j = 0; j < pts.size(); ++j)
		copy(pts[j].begin(), pts[j].end(),
		     sder1.begin() + knmb*ki + j*kdim);

	    /* Evaluate position and cross-tangent curves at second
	       found parameter value.  */
      	    ecurve[2*kk]->point(pts, tpar2, kder);
	    for (size_t j = 0; j < pts.size(); ++j)
		copy(pts[j].begin(), pts[j].end(),
		     spos2.begin() + knmb*ki + j*kdim);
      
	    ecurve[2*kk+1]->point(pts, tpar2, kder);
	    for (size_t j = 0; j < pts.size(); ++j)
		copy(pts[j].begin(), pts[j].end(),
		     sder2.begin() + knmb*ki + j*kdim);

	    /* Evaluate position and both cross-tangents at the corner nr ki.  */
	    Point pt(3);
	    ecurve[2*ki]->point(pt, sstart[ki]);
	    copy(pt.begin(), pt.end(), scorn.begin() + ki*kdim);

	    ecurve[2*kk+1]->point(pt, send[kk]);
	    copy(pt.begin(), pt.end(), scornder1.begin() + ki*kdim);

	    ecurve[2*ki+1]->point(pt, sstart[ki]);
	    copy(pt.begin(), pt.end(), scornder2.begin() + ki*kdim);

	    /* Compute the weigth of the actual blending surface.  */

	    salpha[ki] = tlambi*tlambi*(3.0 - 2.0*tlambi +
					6.0*tlambj*tlambk);
      
	    /* Add the contribution of the value of this blending surface to the
	       value of the ideal surface.   */

	    for (kh=0; kh<kdim; kh++)
		{
		    sp[ki*kdim+kh] = spos1[knmb*ki+kh] + tlambk*sder1[knmb*ki+kh] +
			spos2[knmb*ki+kh] + tlambj*sder2[knmb*ki+kh] -
			scorn[ki*kdim+kh] - tlambj*scornder1[ki*kdim+kh] 
			- tlambk*scornder2[ki*kdim+kh] -
			tlambj*tlambk*etwist[ki*kdim+kh];
	  
		    eval[kh] += salpha[ki]*sp[ki*kdim+kh];
		}
	}
  
    if (ider >= 1)
	{
	    /* Compute first derivatives of the Gregory Charrot function. */

	    double tl1,tl2,tl3;            /* Barycentric coordinates.  */
	    vector<double> sd1alpha(3), sd2alpha(3);  // 1. derivative of weight
						     // functions.
	    vector<double> sd1p(9),sd2p(9);          // 1. derivative of
	                                             // blending functions.
      
	    /* Copy barycentric coordinates to local variables.  */

	    tl1 = ebar[0];
	    tl2 = ebar[1];
	    tl3 = ebar[2];
      
	    /* Compute the 1. derivatives of the weight functions.  */

	    sd1alpha[0] = 6.0*tl1*(1.0 - tl1 - tl1*tl2 + 
				   2.0*tl2*tl3);
	    sd2alpha[0] = 6.0*tl1*tl1*(tl3 - tl2);

	    sd1alpha[1] = 6.0*tl2*tl2*(tl3 - tl1);
	    sd2alpha[1] = 6.0*tl2*(1.0 - tl2 - tl1*tl2 + 
				   2.0*tl1*tl3);

	    sd1alpha[2] = 6.0*tl3*(-1.0 + tl3 + tl2*tl3 - 
				   2.0*tl1*tl2);
	    sd2alpha[2] = 6.0*tl3*(-1.0 + tl3 + tl1*tl3 - 
				   2.0*tl1*tl2);
      
	    /* Compute 1. derivatives of the functions which blends two sides
	       of the region.  */

	    for (kh=0; kh<kdim; kh++)
		{
		    sd1p[kh] = (spos2[kdim+kh] + tl2*sder2[kdim+kh])*sint[2]
			- sder1[kh] + scornder2[kh] + tl2*etwist[kh];
		    sd2p[kh] = (spos2[kdim+kh] + tl2*sder2[kdim+kh])*sint[2] 
			+ (spos1[kdim+kh] + tl3*sder1[kdim+kh])*sint[0]
			+ sder2[kh] - sder1[kh] - scornder1[kh]
			+ scornder2[kh] + (tl2 - tl3)*etwist[kh];

		    sd1p[kdim+kh] = -sder2[knmb+kh] + sder1[knmb+kh]
			- (spos2[knmb+kdim+kh] + tl3*sder2[knmb+kdim+kh])*sint[0] 
			- (spos1[knmb+kdim+kh] + tl1*sder1[knmb+kdim+kh])*sint[1] 
			+ scornder1[kdim+kh] - scornder2[kdim+kh]
			+ (tl1 - tl3)*etwist[kdim+kh];
		    sd2p[kdim+kh] = - sder2[knmb+kh]
			- (spos1[knmb+kdim+kh] + tl1*sder1[knmb+kdim+kh])*sint[1]
			+ scornder1[kdim+kh] + tl1*etwist[kdim+kh];

		    sd1p[2*kdim+kh] = sder2[2*knmb+kh]
			+ (spos1[2*knmb+kdim+kh] + tl2*sder1[2*knmb+kdim+kh])*sint[2]
			- scornder1[2*kdim+kh] - tl2*etwist[2*kdim+kh];
		    sd2p[2*kdim+kh] = sder1[2*knmb+kh]
			- (spos2[2*knmb+kdim+kh] + tl1*sder2[2*knmb+kdim+kh])*sint[1]
			- scornder2[2*kdim+kh] - tl1*etwist[2*kdim+kh];
	  
		    /* Compute the first derivative of the Gregory Charrot
		       function. */

		    for (ki=0; ki<3; ki++)
			{
			    eval[kdim+kh] += sd1alpha[ki]*sp[ki*kdim+kh]
				+ salpha[ki]*sd1p[ki*kdim+kh];

			    eval[2*kdim+kh] += sd2alpha[ki]*sp[ki*kdim+kh]
				+ salpha[ki]*sd2p[ki*kdim+kh];
			}
		}
      
	    if (ider >= 2)
		{
		    /* 2. derivatives of blending patches. */
		    vector<double> sd11alpha(3), sd12alpha(3), sd22alpha(3);
		    /* 2. derivatives of weight function. */
		    vector<double> sd11p(9), sd12p(9), sd22p(9);
	  
		    /* Compute second derivatives of the Gregory Charrot function. */

		    /* Compute the 2. derivatives of the weight functions.  */

		    sd11alpha[0] = 6.0 - 12.0*tl1 
			- 24.0*tl1*tl2 + 12.0*tl2*tl3;
		    sd12alpha[0] = tl1*(12.0 - 18.0*tl1
					- 24.0*tl2);
		    sd22alpha[0] = -12.0*tl1*tl1;
	  
		    sd11alpha[1] = -12.0*tl2*tl2; 
		    sd12alpha[1] = tl2*(12.0 - 18.0*tl2
					- 24.0*tl1);
		    sd22alpha[1] = 6.0 - 12.0*tl2 
			- 24.0*tl1*tl2 + 12.0*tl1*tl3;

		    sd11alpha[2] = 6.0 - 12.0*tl3 
			- 24.0*tl2*tl3 + 12.0*tl1*tl2;
		    sd12alpha[2] = 6.0 + 6.0*tl3*tl3
			+ 12.0*(-tl3 + tl1*tl2 - tl1*tl3 - tl2*tl3);
		    sd22alpha[2] = 6.0 - 12.0*tl3 
			- 24.0*tl1*tl3 + 12.0*tl1*tl2;

		    /* Compute 2. derivatives of the functions which blends two sides
		       of the region.  */

		    for (kh=0; kh<kdim; kh++)
			{
			    sd11p[kh] = (spos2[2*kdim+kh]+
					 tl2*sder2[2*kdim+kh])*sint[2]*sint[2];
			    sd12p[kh] = ((spos2[2*kdim+kh] +
					  tl2*sder2[2*kdim+kh])*sint[2]
					 + sder2[kdim+kh]*sint[2] -
					 sder1[kdim+kh])*sint[0] + etwist[kh];
			    sd22p[kh] = ((spos2[2*kdim+kh] +
					  tl2*sder2[2*kdim+kh])*sint[2]
					 + 2.0*sder2[kdim+kh])*sint[2]
				+ ((spos1[2*kdim+kh] + tl3*sder1[2*kdim+kh])*sint[0]
				   - 2.0*sder1[kdim+kh])*sint[0]
				+ 2.0*etwist[kh];
	      
			    sd11p[kdim+kh] = ((spos2[knmb+2*kdim+kh] 
					       + tl3*sder2[knmb+2*kdim+kh])*sint[0]
					      + 2.0*sder2[knmb+kdim+kh])*sint[0]
				+ ((spos1[knmb+2*kdim+kh] 
				    + tl1*sder1[knmb+2*kdim+kh])*sint[1]
				   - 2.0*sder1[knmb+kdim+kh])*sint[1]
				+ 2.0*etwist[kdim+kh];
			    sd12p[kdim+kh] = ((spos1[knmb+2*kdim+kh] 
					       + tl1*sder1[knmb+2*kdim+kh])*sint[1]
					      - sder1[knmb+kdim+kh])*sint[1]
				+ sder2[knmb+kdim+kh]*sint[0]
				+ etwist[kdim+kh];
			    sd22p[kdim+kh] =
				(spos1[knmb+2*kdim+kh]
				 + tl1*sder1[knmb+2*kdim+kh])*sint[1]*sint[1];
	      
			    sd11p[2*kdim+kh] =
				(spos1[2*knmb+2*kdim+kh]
				 + tl2*sder1[2*knmb+2*kdim+kh])*sint[2]*sint[2];
			    sd12p[2*kdim+kh] =
				-sder2[2*knmb+kdim+kh]*sint[1] 
				+ sder1[2*knmb+kdim+kh]*sint[2] - etwist[2*kdim+kh];
			    sd22p[2*kdim+kh] =
				(spos2[2*knmb+2*kdim+kh]
				 + tl1*sder2[2*knmb+2*kdim+kh])*sint[1]*sint[1];

			    /* Compute the 2. derivative of the Gregory Charrot
			       function. */

			    for (ki=0; ki<3; ki++)
				{
				    eval[3*kdim+kh] += sd11alpha[ki]*sp[ki*kdim+kh]
					+ 2.0*sd1alpha[ki]*sd1p[ki*kdim+kh]
					+ salpha[ki]*sd11p[ki*kdim+kh];
		  
				    eval[4*kdim+kh] += sd12alpha[ki]*sp[ki*kdim+kh]
					+ sd1alpha[ki]*sd2p[ki*kdim+kh] 
					+ sd2alpha[ki]*sd1p[ki*kdim+kh]
					+ salpha[ki]*sd12p[ki*kdim+kh];
		  
				    eval[5*kdim+kh] += sd22alpha[ki]*sp[ki*kdim+kh]
					+ 2.0*sd2alpha[ki]*sd2p[ki*kdim+kh]
					+ salpha[ki]*sd22p[ki*kdim+kh];
				}
			}
		}
	}
  
//     /* Ideal surface evaluated.  */

//     *jstat = kwarn;
//     goto out;


//     /* Error in lower level function.  */

//     error :
// 	*jstat = kstat;
//     goto out;

//     out :
	return;
}


			      
// void
//   sh1467(SISLCurve *ecurve[],double etwist[],int ider,double ebar[],
// 	 double eval[],int *jstat)
void
gregoryCharrotFunction5(vector<shared_ptr<SplineCurve> >&
					   ecurve,
					   vector<double>& etwist, int ider,
					   vector<double>& ebar,
					   vector<double>& eval)
    /*
*********************************************************************
*                                                                   
* PURPOSE    : Given the barycentric coordinates of a point in a 5-sided
*              vertex region, evaluate the value of the ideal blending 
*              surface of the vertex region in this point.
*
*
*
* INPUT      : ecurve - Position and cross-tangent curves around the vertex
*                       region. For each edge of the region position and cross-
*                       tangent curves are given. The curves follow each other
*                       around the region and are oriented counter-clock-wise.
*                       The dimension of the array is 10.
*              etwist - Twist-vectors of the corners of the vertex region. The
*                       first element of the array is the twist in the corner
*                       before the first edge, etc. The dimension of the array
*                       is 5*kdim.
*              ider   - Number of derivatives to compute. Directions of 
*                       differentiation is that of the two first barycentric
*                       coordinates. 0 <= ider <= 2.
*              ebar   - Generalized barycentric coordinates of the point to be 
*                       evaluated. The dimension of the array is 5.
*                       
*
* OUTPUT     : eval   - Value and derivatives of ideal blending surface in the 
*                       parameter directions of the two first coordinates in the 
*                       given point. Dimension of the array is 3*(1+..+(ider+1)).
*              jstat  - status messages  
*                                         > 0      : warning
*                                         = 0      : ok
*                                         < 0      : error
*
*
* METHOD     : A functional description of the ideal surface is given as
*              a blend between five surfaces, each of which fulfill the
*              continuity requirements over two edges.
*
* REFERENCES : Gregory and Charrot : A pentagonal surface patch for
*                                    computer aided geometric design
*
* USE        : 3D geometry only.
*
*-
* CALLS      : s1221 - Evaluate curve at a given parameter value.
*              sh1467_s9fac1() - Compute value and derivative of expression 
*			         in the barycentric coordinates.  
*              sh1467_s9fac2() - Compute value, first and second derivative.
*
* WRITTEN BY : Vibeke Skytt, SI, 04.90.
*
*********************************************************************
*/
{
    int ki,kj,kk,kh;     /* Counters.                                       */
    int kder;            /* Number of derivatives to evaluate.              */
    int kdim = 3;        /* Dimension of geometry.                          */
    int kwarn = 0;       /* Indicates if a warning is to be sendt.          */
    int knmb;            /* Number of doubles pr derivative.                */
    int kl = 0;          /* Number of derivatives in the output array.      */
    double tdum;         /* Product of three barycentric coordinates.       */
    double tdenom;       /* Denominator of weight function.                 */
    vector<double> salpha(5); /* Weigth of actual blending surface.              */
    vector<double> sp(15);    /* Value of function interpolating two sides.      */
    vector<double> sstart(5); /* Start-parameters of edge-curves.                */
    vector<double> send(5);   /* End-parameters of edge-curves.                  */
    vector<double> sint(5);   /* Parameter intervals of edge-curves.             */
    vector<double> sx1(5);    /* First local coordinates of the blending patches.  */
    vector<double> sx2(5);    /* Second local coordinates of the blending patches. */
    vector<double> spar(5);   /* Parameter values of edges.                      */
    vector<double> spos(45);  /* Position of edge-curve.                         */
    vector<double> sder(45);  /* Value of oss tangent.                           */
    vector<double> scorn(15); /* Position of edge-curves in actual corner of 
			         vertex region.                                  */
    vector<double> scornder1(15);/* Tangent in corner along next edge.           */
    vector<double> scornder2(15); /* Tangent in corner along previous edge       */

    /* Test input.  */

    if (ider > 2) kwarn = 1;
  
    /* Initialise.  */

    kder = ider;
    knmb = kdim*(ider+1);
    for (ki=0; ki<ider; ki++) kl += ider + 1;
  
    /* Initiate output array to zero.  */

    for (kh=0; kh<kl*kdim; kh++) eval[kh] = 0.0;
        
    for (ki=0; ki<5; ki++)
	{

	    /* Get endpoints of parameter intervals of edge curves.  */

	    sstart[ki] = *(ecurve[2*ki]->basis().begin() + ecurve[2*ki]->order()-1);
	    send[ki] = *(ecurve[2*ki]->basis().begin() + ecurve[2*ki]->numCoefs());
	    sint[ki] = send[ki] - sstart[ki];

	    /* Compute local coordinates corresponding to one blending patch. */
      
	    kj = (ki + 2) % 5;
	    kk = (ki > 0) ? ki-1 : 4;
	    sx1[ki] = ebar[kj]/(ebar[kk] + ebar[kj]);
      
	    kj = (ki > 1) ? ki-2 : 3+ki;
	    kk = (ki + 1) % 5;
	    sx2[ki] = ebar[kj]/(ebar[kk] + ebar[kj]);
      
	    /* Compute parameter value corresponding to the current edge. */

	    spar[kj] = (1.0 - sx1[ki])*sstart[ki] + sx1[ki]*send[ki];

	    /* Evaluate position and cross-tangent curves at points on
	       the edges needed when evaluating surface.        */

	    /* Evaluate position and cross-tangent curves at parameter value.  */
	    vector<Point> pts(kder+1);
	    ecurve[2*ki]->point(pts, spar[kj], kder);
	    for (size_t j = 0; j < pts.size(); ++j)
		copy(pts[j].begin(), pts[j].end(),
		     spos.begin() + knmb*kj + j*kdim);
      
	    ecurve[2*ki+1]->point(pts, spar[kj], kder);
	    for (size_t j = 0; j < pts.size(); ++j)
		copy(pts[j].begin(), pts[j].end(),
		     sder.begin() + knmb*kj + j*kdim);
 
	    /* Evaluate position and both cross-tangents at the corner nr ki.  */

	    kk = (ki > 0) ? ki-1 : 4;
	    Point pt(3);
	    ecurve[2*ki]->point(pt, sstart[ki]);
	    copy(pt.begin(), pt.end(), scorn.begin() + ki*kdim);

	    ecurve[2*kk+1]->point(pt, send[kk]);
	    copy(pt.begin(), pt.end(), scornder1.begin() + ki*kdim);

	    ecurve[2*ki+1]->point(pt, sstart[ki]);
	    copy(pt.begin(), pt.end(), scornder2.begin() + ki*kdim);
	}
  
    /* Compute position of the Gregory Charrot patch at the given parameter
       value.  First compute factors in the expression of the weight function. */

    tdenom = 0.0;
    for (kh=0; kh<5; kh++)
	{
	    kj = (kh + 1) % 5;
	    kk = (kh > 0) ? kh-1 : 4;	  

	    tdum = ebar[kk]*ebar[kh]*ebar[kj];
	    tdenom += tdum*tdum;
	}

    for (ki=0; ki<5; ki++)
	{

	    kj = (ki + 1) % 5;
	    kk = (ki > 0) ? ki-1 : 4;

	    /* Compute the weigth of the actual blending surface.  */

	    tdum = ebar[kk]*ebar[ki]*ebar[kj];
	    salpha[ki] = tdum*tdum/tdenom;

	    /* Compute the blending surface.  */

	    kj = (ki + 2) % 5;
	    kk = (ki > 1) ? ki-2 : 3+ki;

	    for (kh=0; kh<kdim; kh++)
		{
		    sp[ki*kdim+kh] = spos[knmb*kj+kh] + sx1[ki]*sder[knmb*kj+kh] +
			spos[knmb*kk+kh] + sx2[ki]*sder[knmb*kk+kh] -
			scorn[ki*kdim+kh] - sx1[ki]*scornder1[ki*kdim+kh] 
			- sx2[ki]*scornder2[ki*kdim+kh] -
			sx1[ki]*sx2[ki]*etwist[ki*kdim+kh];
	  
		    /* Add the contribution of the value of this blending surface
		       to the value of the ideal surface.   */

		    eval[kh] += salpha[ki]*sp[ki*kdim+kh];
		}
	}
  
    if (ider >= 1)
	{
	    /* Compute first derivatives of the Gregory Charrot function. */

	    double tconst = (sqrt(5.0) - 1.0)/2.0; /* Constant. */
	    double tdum1,tdum2; /* Factor in expression for derivative of
				   weight function.*/
	    double tnom1,tnom2; /* Factor in expression for derivative of
				   weight function.*/
	    vector<double> sd1alpha(5), sd2alpha(5); /* 1. derivative of
						       weight functions. */
	    vector<double> sd1p(15), sd2p(15);       /* 1. derivative of
						       blending functions. */
	    vector<double> sd1bar(5), sd2bar(5);     /* 1. derivative of
						       barycentric coordinates.*/
	    vector<double> sd1x1(5), sd2x1(5);       /* 1. derivative of
						       local coordinates.  */      
	    vector<double> sd1x2(5), sd2x2(5);       /* 1. derivative of
						       local coordinates.  */      
	    vector<double> sd1par(5), sd2par(5);     /* 1. derivative of
						       parameter value at edge. */

	    /* Differentiate barycentric coordinates.  */

	    sd1bar[0] = 1.0;
	    sd1bar[1] = 0.0;
	    sd1bar[2] = -1.0;
	    sd1bar[3] = -tconst;
	    sd1bar[4] = tconst;

	    sd2bar[0] = 0.0;
	    sd2bar[1] = 1.0;
	    sd2bar[2] = tconst;
	    sd2bar[3] = -tconst;
	    sd2bar[4] = -1.0;

	    for (ki=0; ki<5; ki++)
		{
		    /* Differentiate local coordinates.  */      

		    kj = (ki + 2) % 5;
		    kk = (ki > 0) ? ki-1 : 4;
		    tdum = ebar[kk] + ebar[kj];
		    sd1x1[ki] = sd1bar[kj]/tdum
			- ebar[kj]*(sd1bar[kk] + sd1bar[kj])/(tdum*tdum);
		    sd2x1[ki] = sd2bar[kj]/tdum
			- ebar[kj]*(sd2bar[kk] + sd2bar[kj])/(tdum*tdum);

		    kj = (ki > 1) ? ki-2 : 3+ki;
		    kk = (ki + 1) % 5;
		    tdum = ebar[kk] + ebar[kj];
		    sd1x2[ki] = sd1bar[kj]/tdum
			- ebar[kj]*(sd1bar[kk] + sd1bar[kj])/(tdum*tdum);
		    sd2x2[ki] = sd2bar[kj]/tdum
			- ebar[kj]*(sd2bar[kk] + sd2bar[kj])/(tdum*tdum);

		    /* Compute derivative of parameter value corresponding to
		       the current edge.  */

		    sd1par[kj] = sint[ki]*sd1x1[ki];
		    sd2par[kj] = sint[ki]*sd2x1[ki];
		}
      
	    /* Compute factors used in the expression for the 1. derivatives of
	       the weight function.  */

	    tnom1 = tnom2 = 0.0;
	    for (kh=0; kh<5; kh++)
		{
		    kj = (kh + 1) % 5;
		    kk = (kh > 0) ? kh-1 : 4;	  
	      
		    sh1467_s9fac1(ebar,sd1bar,sd2bar,kk,kh,kj,&tdum,&tdum1,&tdum2);

		    tnom1 += 2.0*tdum*tdum1;
		    tnom2 += 2.0*tdum*tdum2;
		}
      

	    for (ki=0; ki<5; ki++)
		{
		    /* Compute the 1. derivatives of the weight functions.  */

		    kj = (ki + 1) % 5;
		    kk = (ki > 0) ? ki-1 : 4;	  
	  
		    sh1467_s9fac1(ebar,sd1bar,sd2bar,kk,ki,kj,&tdum,&tdum1,&tdum2);

		    sd1alpha[ki] = tdum*(2.0*tdum1 - tdum*tnom1/tdenom)/tdenom;
		    sd2alpha[ki] = tdum*(2.0*tdum2 - tdum*tnom2/tdenom)/tdenom;
      
		    /* Compute 1. derivatives of the functions which blends two sides
		       of the region.  */

		    kj = (ki + 2) % 5;
		    kk = (ki > 1) ? ki-2 : 3+ki;	  

		    for (kh=0; kh<kdim; kh++)
			{
			    sd1p[ki*kdim+kh] = spos[kj*knmb+kdim+kh]*sd1par[kj]
				+ sd1x1[ki]*sder[kj*knmb+kh] 
				+ sx1[ki]*sder[kj*knmb+kdim+kh]*sd1par[kj]
				+ spos[kk*knmb+kdim+kh]*sd1par[kk]
				+ sd1x2[ki]*sder[kk*knmb+kh]
				+ sx2[ki]*sder[kk*knmb+kdim+kh]*sd1par[kk]
				- sd1x1[ki]*scornder1[ki*kdim+kh]
				- sd1x2[ki]*scornder2[ki*kdim+kh]
				- (sd1x1[ki]*sx2[ki] 
				   + sx1[ki]*sd1x2[ki])*etwist[ki*kdim+kh];

			    sd2p[ki*kdim+kh] = spos[kj*knmb+kdim+kh]*sd2par[kj]
				+ sd2x1[ki]*sder[kj*knmb+kh] 
				+ sx1[ki]*sder[kj*knmb+kdim+kh]*sd2par[kj]
				+ spos[kk*knmb+kdim+kh]*sd2par[kk]
				+ sd2x2[ki]*sder[kk*knmb+kh]
				+ sx2[ki]*sder[kk*knmb+kdim+kh]*sd2par[kk]
				- sd2x1[ki]*scornder1[ki*kdim+kh]
				- sd2x2[ki]*scornder2[ki*kdim+kh]
				- (sd2x1[ki]*sx2[ki] 
				   + sx1[ki]*sd2x2[ki])*etwist[ki*kdim+kh];

			    /* Add the contribution from this blending surface to
			       the first derivative of the Gregory Charrot function*/

			    eval[kdim+kh] += sd1alpha[ki]*sp[ki*kdim+kh]
				+ salpha[ki]*sd1p[ki*kdim+kh]; 

			    eval[2*kdim+kh] += sd2alpha[ki]*sp[ki*kdim+kh]
				+ salpha[ki]*sd2p[ki*kdim+kh]; 
			}
		}
      
	    if (ider >= 2)
		{
		    /* 2. derivatives of weight function.  */
		    vector<double> sd11alpha(5), sd12alpha(5), sd22alpha(5);
		    /* 2. derivatives of blending surfaces.                  */
		    vector<double> sd11p(15), sd12p(15), sd22p(15);
		    /* 2. derivatives of local coordinates.               */
		    vector<double> sd11x1(5), sd12x1(5), sd22x1(5);
		    /* 2. derivatives of local coordinates.               */
		    vector<double> sd11x2(5), sd12x2(5), sd22x2(5);
		    /* 2. derivatives of parameter value at edge. */
		    vector<double> sd11par(5), sd12par(5), sd22par(5);
								
		    double tdum11,tdum12,tdum22;  /* Factor in expression for 2.
						     derivative of weight function.*/
		    double tnom11,tnom12,tnom22;  /* Factor in expression for 2.
						     derivative of weight function.*/
	  
		    /* Compute second derivatives of the Gregory Charrot function. */

		    for (ki=0; ki<5; ki++)
			{
			    /* Differentiate local coordinates.  */      

			    kj = (ki + 2) % 5;
			    kk = (ki > 0) ? ki-1 : 4;
			    tdum = ebar[kk] + ebar[kj];
			    tdum1 = sd1bar[kk] + sd1bar[kj];
			    tdum2 = sd2bar[kk] + sd2bar[kj];
	  
			    sd11x1[ki] = 2.0*(ebar[kj]*tdum1*tdum1/tdum
					      - sd1bar[kj]*tdum1)/(tdum*tdum);
			    sd12x1[ki] =
				(2.0*ebar[kj]*tdum1*tdum2/tdum
				 - sd1bar[kj]*tdum2 - sd2bar[kj]*tdum1)/(tdum*tdum);
			    sd22x1[ki] = 2.0*(ebar[kj]*tdum2*tdum2/tdum
					      - sd2bar[kj]*tdum2)/(tdum*tdum);

			    kj = (ki > 1) ? ki-2 : 3+ki;
			    kk = (ki + 1) % 5;
			    tdum = ebar[kk] + ebar[kj];
			    tdum1 = sd1bar[kk] + sd1bar[kj];
			    tdum2 = sd2bar[kk] + sd2bar[kj];
	  
			    sd11x2[ki] = 2.0*(ebar[kj]*tdum1*tdum1/tdum
					      - sd1bar[kj]*tdum1)/(tdum*tdum);
			    sd12x2[ki] =
				(2.0*ebar[kj]*tdum1*tdum2/tdum
				 - sd1bar[kj]*tdum2 - sd2bar[kj]*tdum1)/(tdum*tdum);
			    sd22x2[ki] = 2.0*(ebar[kj]*tdum2*tdum2/tdum
					      - sd2bar[kj]*tdum2)/(tdum*tdum);

			    /* Compute derivative of parameter value corresponding to
			       the current edge.  */

			    sd11par[kj] = sint[ki]*sd11x1[ki];
			    sd12par[kj] = sint[ki]*sd12x1[ki];
			    sd22par[kj] = sint[ki]*sd22x1[ki];
			}

		    /* Compute factors used in the expression for the
		       2. derivatives of the weight function.  */

		    tnom11 = tnom12 = tnom22 = 0.0;
		    for (kh=0; kh<5; kh++)
			{
			    kj = (kh + 1) % 5;
			    kk = (kh > 0) ? kh-1 : 4;	  
	      
			    sh1467_s9fac2(ebar,sd1bar,sd2bar,kk,kh,kj,
					  &tdum,&tdum1,&tdum2,
					  &tdum11,&tdum12,&tdum22);
	      
			    tnom11 += 2.0*(tdum1*tdum1 + tdum*tdum11);
			    tnom12 += 2.0*(tdum1*tdum2 + tdum*tdum12);
			    tnom22 += 2.0*(tdum2*tdum2 + tdum*tdum22);

			}

		    for (ki=0; ki<5; ki++)
			{
	      
			    /* Compute the 2. derivatives of the weight functions. */

			    kj = (ki + 1) % 5;
			    kk = (ki > 0) ? ki-1 : 4;	  
	      
			    sh1467_s9fac2(ebar,sd1bar,sd2bar,kk,ki,kj,
					  &tdum,&tdum1,&tdum2,
					  &tdum11,&tdum12,&tdum22);
	      
			    sd11alpha[ki] = 2.0*(tdum1*tdum1 + tdum*tdum11)/tdenom
				- 4.0*tdum*tdum1*tnom1/(tdenom*tdenom)
				- tdum*tdum*tnom11/(tdenom*tdenom) 
				+ 2.0*tdum*tdum*tnom1*tnom1/(tdenom*tdenom*tdenom);

			    sd12alpha[ki] =
				2.0*(tdum1*tdum2 + tdum*tdum12)/tdenom
				- 2.0*tdum*(tdum1*tnom2 +
					    tdum2*tnom1)/(tdenom*tdenom)
				- tdum*tdum*tnom12/(tdenom*tdenom)
				+ 2.0*tdum*tdum*tnom1*tnom2/(tdenom*tdenom*tdenom);

			    sd22alpha[ki] = 2.0*(tdum2*tdum2 + tdum*tdum22)/tdenom
				- 4.0*tdum*tdum2*tnom2/(tdenom*tdenom)
				- tdum*tdum*tnom22/(tdenom*tdenom)
				+ 2.0*tdum*tdum*tnom2*tnom2/(tdenom*tdenom*tdenom);
	      
			    /* Compute 2. derivatives of the functions which blends
			       two sides of the region.  */

			    kj = (ki + 2) % 5;
			    kk = (ki > 1) ? ki - 2 : 3 + ki;
	      
			    for (kh=0; kh<kdim; kh++)
				{
				    sd11p[ki*kdim+kh] =
					spos[kj*knmb+2*kdim+kh]*sd1par[kj]*sd1par[kj]
					+ spos[kj*knmb+kdim+kh]*sd11par[kj]
					+ sd11x1[ki]*sder[kj*knmb+kh]
					+ 2.0*sd1x1[ki]*
					sder[kj*knmb+kdim+kh]*sd1par[kj]
					+ sx1[ki]*sder[kj*knmb+2*kdim+kh]*
					sd1par[kj]*sd1par[kj]
					+ sx1[ki]*sder[kj*knmb+kdim+kh]*sd11par[kj]
					+ 2.0*sd1x2[ki]*sder[kk*knmb+kdim+kh]*
					sd1par[kk]
					+ sx2[ki]*sder[kk*knmb+2*kdim+kh]*
					sd1par[kk]*sd1par[kk]
					+ spos[kk*knmb+2*kdim+kh]*sd1par[kk]*
					sd1par[kk]
					+ spos[kk*knmb+kdim+kh]*sd11par[kk]
					+ sd11x2[ki]*sder[kk*knmb+kh]
					+ sx2[ki]*sder[kk*knmb+kdim+kh]*sd11par[kk]
					- sd11x1[ki]*scornder1[ki*kdim+kh]
					- sd11x2[ki]*scornder2[ki*kdim+kh]
					- (sd11x1[ki]*sx2[ki] 
					   + 2.0*sd1x1[ki]*sd1x2[ki]
					   + sx1[ki]*sd11x2[ki])
					*etwist[ki*kdim+kh];
	      
				    sd12p[ki*kdim+kh] =
					spos[kj*knmb+2*kdim+kh]*sd1par[kj]*sd2par[kj]
					+ spos[kj*knmb+kdim+kh]*sd12par[kj]
					+ sd12x1[ki]*sder[kj*knmb+kh]
					+ sd1x1[ki]*sder[kj*knmb+kdim+kh]*sd2par[kj]
					+ sd2x1[ki]*sder[kj*knmb+kdim+kh]*sd1par[kj]
					+ sx1[ki]*sder[kj*knmb+2*kdim+kh]*
					sd1par[kj]*sd2par[kj]
					+ sx1[ki]*sder[kj*knmb+kdim+kh]*sd12par[kj]
					+ sx2[ki]*sder[kk*knmb+2*kdim+kh]*
					sd1par[kk]*sd2par[kk]
					+ spos[kk*knmb+2*kdim+kh]*sd1par[kk]*
					sd2par[kk]
					+ spos[kk*knmb+kdim+kh]*sd12par[kk]
					+ sd12x2[ki]*sder[kk*knmb+kh]
					+ sd1x2[ki]*sder[kk*knmb+kdim+kh]*sd2par[kk]
					+ sd2x2[ki]*sder[kk*knmb+kdim+kh]*sd1par[kk]
					+ sx2[ki]*sder[kk*knmb+kdim+kh]*sd12par[kk]
					- sd12x1[ki]*scornder1[ki*kdim+kh]
					- sd12x2[ki]*scornder2[ki*kdim+kh]
					- (sd12x1[ki]*sx2[ki] 
					   + sd1x1[ki]*sd2x2[ki]
					   + sd2x1[ki]*sd1x2[ki]
					   + sx1[ki]*sd12x2[ki])
					*etwist[ki*kdim+kh];

				    sd22p[ki*kdim+kh] =
					spos[kj*knmb+2*kdim+kh]*sd2par[kj]*sd2par[kj]
					+ spos[kj*knmb+kdim+kh]*sd22par[kj]
					+ sd22x1[ki]*sder[kj*knmb+kh]
					+ 2.0*sd2x1[ki]*sder[kj*knmb+kdim+kh]*
					sd2par[kj]
					+ sx1[ki]*sder[kj*knmb+2*kdim+kh]*
					sd2par[kj]*sd2par[kj]
					+ sx1[ki]*sder[kj*knmb+kdim+kh]*sd22par[kj]
					+ 2.0*sd2x2[ki]*sder[kk*knmb+kdim+kh]*
					sd2par[kk]
					+ sx2[ki]*sder[kk*knmb+2*kdim+kh]*
					sd2par[kk]*sd2par[kk]
					+ spos[kk*knmb+2*kdim+kh]*sd2par[kk]*
					sd2par[kk]
					+ spos[kk*knmb+kdim+kh]*sd22par[kk]
					+ sd22x2[ki]*sder[kk*knmb+kh]
					+ sx2[ki]*sder[kk*knmb+kdim+kh]*sd22par[kk]
					- sd22x1[ki]*scornder1[ki*kdim+kh]
					- sd22x2[ki]*scornder2[ki*kdim+kh]
					- (sd22x1[ki]*sx2[ki] 
					   + 2.0*sd2x1[ki]*sd2x2[ki]
					   + sx1[ki]*sd22x2[ki])
					*etwist[ki*kdim+kh];
	      

				    /* Add the current contribution to the 2.
				       derivative of the Gregory Charrot function. */

				    eval[3*kdim+kh] += sd11alpha[ki]*sp[ki*kdim+kh]
					+ 2.0*sd1alpha[ki]*sd1p[ki*kdim+kh]
					+ salpha[ki]*sd11p[ki*kdim+kh];
		  
				    eval[4*kdim+kh] += sd12alpha[ki]*sp[ki*kdim+kh]
					+ sd1alpha[ki]*sd2p[ki*kdim+kh] 
					+ sd2alpha[ki]*sd1p[ki*kdim+kh]
					+ salpha[ki]*sd12p[ki*kdim+kh];
		  
				    eval[5*kdim+kh] += sd22alpha[ki]*sp[ki*kdim+kh]
					+ 2.0*sd2alpha[ki]*sd2p[ki*kdim+kh]
					+ salpha[ki]*sd22p[ki*kdim+kh];
				}
			}
		}
	}
  
//     /* Ideal surface evaluated.  */

//     *jstat = kwarn;
//     goto out;


//     /* Error in lower level function.  */

//     error :
// 	*jstat = kstat;
//     goto out;

//     out :
	return;
}



// static void
//   sh1467_s9fac1(double ebar[],double ed1bar[],double ed2bar[],
// 		int i1,int i2,int i3,double *cfac,double *cfac1,
// 		double *cfac2)
void sh1467_s9fac1(vector<double>& ebar, vector<double>& ed1bar,
				      vector<double>& ed2bar, int i1, int i2, int i3,
				      double *cfac, double *cfac1, double *cfac2)
    /*
*********************************************************************
*                                                                   
* PURPOSE    : Compute value and 1. derivative of product of barycentric
*              coordinates.
*
*
*
* INPUT      : ebar   - Array containing value of barycentric coordinates.
*              ed1bar - Array containing derivative in first parameter 
*                       direction of barycentric coordinates.
*              ed2bar - Array containing derivative in second parameter 
*                       direction of barycentric coordinates.
*              i1     - Index of first barycentric coordinate in product.
*              i2     - Index of second barycentric coordinate in product.
*              i3     - Index of third barycentric coordinate in product.
*                       
*
* OUTPUT     : cfac   - Product of barycentric coordinates.
*              cfac1  - Derivative in first parameter direction of product.
*              cfac2  - Derivative in second parameter direction of product.
*
*
* WRITTEN BY : Vibeke Skytt, SI, 04.90.
*
*********************************************************************
*/
{
    *cfac = ebar[i1]*ebar[i2]*ebar[i3];
    *cfac1 = ed1bar[i1]*ebar[i2]*ebar[i3] + ebar[i1]*ed1bar[i2]*ebar[i3]
	+ ebar[i1]*ebar[i2]*ed1bar[i3];
    *cfac2 = ed2bar[i1]*ebar[i2]*ebar[i3] + ebar[i1]*ed2bar[i2]*ebar[i3]
	+ ebar[i1]*ebar[i2]*ed2bar[i3];

    return;
}

  

// static void
//   sh1467_s9fac2(double ebar[],double ed1bar[],double ed2bar[],
// 		int i1,int i2,int i3,double *cfac,double *cfac1,
// 		double *cfac2,double *cfac11,double *cfac12,
// 		double *cfac22)
void sh1467_s9fac2(vector<double>& ebar, vector<double>& ed1bar,
				      vector<double>& ed2bar, int i1, int i2, int i3,
				      double *cfac, double *cfac1, double *cfac2,
				      double *cfac11, double *cfac12, double *cfac22)
    /*
*********************************************************************
*                                                                   
* PURPOSE    : Compute value and 1. and 2. derivatives of product of 
*              barycentric coordinates.
*
*
*
* INPUT      : ebar   - Array containing value of barycentric coordinates.
*              ed1bar - Array containing derivative in first parameter 
*                       direction of barycentric coordinates.
*              ed2bar - Array containing derivative in second parameter 
*                       direction of barycentric coordinates.
*              i1     - Index of first barycentric coordinate in product.
*              i2     - Index of second barycentric coordinate in product.
*              i3     - Index of third barycentric coordinate in product.
*                       
*
* OUTPUT     : cfac   - Product of barycentric coordinates.
*              cfac1  - Derivative in first parameter direction of product.
*              cfac2  - Derivative in second parameter direction of product.
*              cfac11 - 2. derivative in first parameter direction of product.
*              cfac11 - Mixed derivative of product.
*              cfac22 - 2. derivative in second parameter direction of product.
*
*
* WRITTEN BY : Vibeke Skytt, SI, 04.90.
*
*********************************************************************
*/
{
  
    *cfac = ebar[i1]*ebar[i2]*ebar[i3];
    *cfac1 = ed1bar[i1]*ebar[i2]*ebar[i3] + ebar[i1]*ed1bar[i2]*ebar[i3]
	+ ebar[i1]*ebar[i2]*ed1bar[i3];
    *cfac2 = ed2bar[i1]*ebar[i2]*ebar[i3] + ebar[i1]*ed2bar[i2]*ebar[i3]
	+ ebar[i1]*ebar[i2]*ed2bar[i3];
    *cfac11 = ed1bar[i1]*(ed1bar[i2]*ebar[i3] + ebar[i2]*ed1bar[i3])
	+ ed1bar[i2]*(ed1bar[i1]*ebar[i3] + ebar[i1]*ed1bar[i3])
	+ ed1bar[i3]*(ed1bar[i1]*ebar[i2] + ebar[i1]*ed1bar[i2]);
    *cfac12 = ed1bar[i1]*(ed2bar[i2]*ebar[i3] + ebar[i2]*ed2bar[i3])
	+ ed1bar[i2]*(ed2bar[i1]*ebar[i3] + ebar[i1]*ed2bar[i3])
	+ ed1bar[i3]*(ed2bar[i1]*ebar[i2] + ebar[i1]*ed2bar[i2]);
    *cfac22 = ed2bar[i1]*(ed2bar[i2]*ebar[i3] + ebar[i2]*ed2bar[i3])
	+ ed2bar[i2]*(ed2bar[i1]*ebar[i3] + ebar[i1]*ed2bar[i3])
	+ ed2bar[i3]*(ed2bar[i1]*ebar[i2] + ebar[i1]*ed2bar[i2]);
}


void
midpoint3(vector<shared_ptr<SplineCurve> >& curves, int icurv,
			     vector<double>& etwist, vector<double>& etang,
			     vector<double>& eder)
    /*
*********************************************************************
*                                                                   
* PURPOSE    : Given a three sided vertex region, evaluate the first
*              blending surface in the corner lying in the middle of the
*              vertex region. Compute the tangent vectors in the middle 
*              vertex along the inner boundaries of the region.
*
*
*
* INPUT      : fshape  - Application driven routine that gives the user an
*                        ability to change the middle point of the region
*                        (the vertex at which the blending surfaces meet),
*                        and the tangent vectors in the middle point along
*                        the curves which divedes the region. 
*              vboundc - Position and cross-tangent curves around the vertex
*                        region. For each edge of the region position and cross-
*                        tangent curves are given. The curves follow each other
*                        around the region and are oriented counter-clock-wise.
*                        The dimension of the array is 6, i.e. 2*icurv.
*              icurv   - Number of sides. icurv = 3.
*              etwist  - Twist-vectors of the corners of the vertex region. The
*                        first element of the array is the twist in the corner
*                        before the first edge, etc. The dimension of the array
*                        is 3*kdim.
*                       
*
* OUTPUT     : etang   - Tangent vectors at the midpoint of the vertex region.
*                        The dimension is icurv*idim.
*              eder    - Value, first and second derivative of the first blending
*                        surface in the corner at the midpoint. The sequence is the
*                        following : Value, 1. derivative in 1. parameter direction,
*                        1. derivative in the 2. parameter direction, 2. derivative
*                        in the 1. parameter direction, mixed derivative and 2.
*                        derivative in the 2. parameter direction. Dimension 6*idim.
*              jstat   - status messages  
*                                         > 0      : warning
*                                         = 0      : ok
*                                         < 0      : error
*
*
* METHOD     : Evaluate the Gregory Charrot function in the midpoint of the
*              vertex region. Compute the wanted derivatives.
*
* REFERENCES : 
*
* USE        : 3D geometry only.
*
*-
* CALLS      : sh1466  - Evaluate the Gregory Charrot function.  
*
* WRITTEN BY : Vibeke Skytt, SI, 03.90.
*
*********************************************************************
*/
{
    ALWAYS_ERROR_IF(icurv != 3,
		"Routine assumes polygonial surface is 3-sided!");

    int kdim = 3;          /* Dimension of geometry space.  */
    ALWAYS_ERROR_IF((int(etwist.size()) != icurv*kdim) || (int(etang.size()) != icurv*kdim) ||
		(int(eder.size()) < 6*kdim),
		"Wrong size of input vector!");

    int kder = 2;          /* Number of derivatives to evaluate.  */
    int ki;                /* Counter.  */
    double tonethird = 1.0/3.0;  /* 1/3     */
    double tonesixth = 1.0/6.0;  /* 1/6     */
    vector<double> sbar(3);        /* Barycentric coordinates of point to evaluate. */
    vector<double> sder(18);       /* Value and derivatives of blending. */

    /* Set up the barycentric coordinates of the midpoint of the region. */

    sbar[0] = sbar[1] = sbar[2] = tonethird;
  
    /* Evaluate the Gregory Charrot function at the midpoint. */

    //    sh1466(curves,etwist,kder,sbar,sder);
    gregoryCharrotFunction3(curves, etwist, kder, sbar, sder);
    //    if (kstat < 0) goto error;
   
    /* Compute tangent vectors.  */

    for (ki=0; ki<kdim; ki++)
	{
	    etang[ki] = -sder[kdim+ki]*tonethird + sder[2*kdim+ki]*tonesixth;
	    //	    etang[ki] *= 0.4;
	    etang[kdim+ki] = sder[kdim+ki]*tonesixth - sder[2*kdim+ki]*tonethird;
	    //	    etang[kdim+ki] *= 0.4;
	    etang[2*kdim+ki] = sder[kdim+ki]*tonesixth + sder[2*kdim+ki]*tonesixth;
	    //	    etang[2*kdim+ki] *= 0.4;
	}

// Ported version does not give the user any choice!
//     /* Application driven routine to alter the midpoint and tangents in the
//        midpoint.  */

//     fshape(sder,etang,kdim,icurv,&kstat);
//     if (kstat < 0) goto error;
  

    /* Copy value and 1. derivatives of first patch.  */

//     memcopy(eder,sder,kdim,DOUBLE);
	    copy(sder.begin(), sder.begin() + kdim, eder.begin());
//     memcopy(eder+kdim,etang+2*kdim,kdim,DOUBLE);
	    copy(etang.begin() + 2*kdim, etang.begin() + 3*kdim,
		 eder.begin() + kdim);
//     memcopy(eder+2*kdim,etang,kdim,DOUBLE);
	    copy(etang.begin(), etang.begin() + kdim, eder.begin() + 2*kdim);
  
    /* Compute 2. derivatives.  */

    for (ki=0; ki<kdim; ki++)
	{
	    eder[3*kdim+ki] = sder[3*kdim+ki]*tonesixth*tonesixth 
		+ 2.0*sder[4*kdim+ki]*tonesixth*tonesixth 
		+ sder[5*kdim+ki]*tonesixth*tonesixth;
	    //	    eder[3*kdim+ki] *= 0.1;
	    eder[4*kdim+ki] = -sder[3*kdim+ki]*tonesixth*tonethird 
		+ sder[4*kdim+ki]*tonesixth*(tonesixth - tonethird)
		+ sder[5*kdim+ki]*tonesixth*tonesixth;
	    //	    eder[4*kdim+ki] *= 0.1;
	    eder[5*kdim+ki] = sder[3*kdim+ki]*tonethird*tonethird 
		- 2.0*sder[4*kdim+ki]*tonethird*tonesixth 
		+ sder[5*kdim+ki]*tonesixth*tonesixth;
	    //	    eder[5*kdim+ki] *= 0.1;
	}
  
//     *jstat = 0;
//     goto out;
  
//     /* Error in a lower level function.  */

//  error:
//     *jstat = kstat;
//     goto out;
  
//     out :
	return;
}



// void
//   sh1463(fshapeProc fshape,
// 	 SISLCurve *vboundc[],int icurv,double etwist[],
// 	 double etang[],double eder[],int *jstat)
void
midpoint4(vector<shared_ptr<SplineCurve> >& curves, int icurv,
			     vector<double>& etwist, vector<double>& etang,
			     vector<double>& eder,
			     double neighbour_tol, double kink_tol)
/*
*********************************************************************
*                                                                   
* PURPOSE    : Given a four sided vertex region, evaluate the first
*              blending surface in the corner lying in the middle of the
*              vertex region. Compute the tangent vectors in the middle 
*              vertex along the inner boundaries of the region.
*
*
*
* INPUT      : fshape  - Application driven routine that gives the user an
*                        ability to change the middle point of the region
*                        (the vertex at which the blending surfaces meet),
*                        and the tangent vectors in the middle point along
*                        the curves which divedes the region. 
*              vboundc - Position and cross-tangent curves around the vertex
*                        region. For each edge of the region position and cross-
*                        tangent curves are given. The curves follow each other
*                        around the region and are oriented counter-clock-wise.
*                        The dimension of the array is 8.
*              icurv   - Number of sides. icurv = 4.
*              etwist  - Twist-vectors of the corners of the vertex region. The
*                        first element of the array is the twist in the corner
*                        before the first edge, etc. The dimension of the array
*                        is 4*kdim.
*                       
*
* OUTPUT     : etang   - Tangent vectors at the midpoint of the vertex region.
*                        The dimension is icurv*idim.
*              eder    - Value, first and second derivative of the first blending
*                        surface in the corner at the midpoint. The sequence is the
*                        following : Value, 1. derivative in 1. parameter direction,
*                        1. derivative in the 2. parameter direction, 2. derivative
*                        in the 1. parameter direction, mixed derivative and 2.
*                        derivative in the 2. parameter direction. Dimension 6*idim.
*              jstat   - status messages  
*                                         > 0      : warning
*                                         = 0      : ok
*                                         < 0      : error
*
*
* METHOD     : Perform rectangular blending of the whole region. Evaluate
*              the resulting surface in the midpoint.
*
* REFERENCES : 
*
* USE        : 3D geometry only.
*
*-
* CALLS      : s1401 - Rectangular blending. 
*             s1421        - Surface evaluation.   
*             s1706        - Turn orientation of curve. 
*             freeCurve    - Free space occupied by a curve.
*             newCurve     - Create new curve.     
*
* WRITTEN BY : Vibeke Skytt, SI, 05.90.
*
*********************************************************************
*/
{
    ALWAYS_ERROR_IF(icurv != 4,
		"Routine assumes polygonial surface is 4-sided!");

    int kdim = 3;          /* Dimension of geometry space.  */
    ALWAYS_ERROR_IF((int(etwist.size()) != icurv*kdim) || (int(etang.size()) != icurv*kdim) ||
		(int(eder.size()) < 6*kdim),
		"Wrong size of input vector!");

    int kder = 2;      /* Number of derivatives of surface to evaluate. */
    int ki;            /* Counters.         */
    vector<double> spar(2);    /* Parameter value of midpoint of rectanular patch. */
    vector<double> sder(18);   /* Value and derivatives of patch in the midpoint.  */
    vector<double> snorm(3);   /* Normal of patch in the midpoint.      */
    SplineSurface* qsurf = new SplineSurface(); // = SISL_NULL;
                            /* Rectangular blending patch.  */
    vector<shared_ptr<ParamCurve> > qc_bnd(4); /* Copy of edge curves.         */
    vector<shared_ptr<ParamCurve> > qc_cross(4); /* Copy of edge curves.       */
    shared_ptr<SplineCurve> qpt_bnd;          /* Pointer to curve.            */
    shared_ptr<SplineCurve> qpt_cross;          /* Pointer to curve.            */

    /* Make a copy of the edge curves.  */

    for (ki=0; ki<4; ki++)
	{
	    qpt_bnd = curves[2*ki];
	    qpt_cross = curves[2*ki+1];
      
	    /* Test dimension of curves. */

	    //	    if (qpt->idim != kdim) goto err104;
	    ALWAYS_ERROR_IF((qpt_bnd->dimension() != kdim) ||
			(qpt_cross->dimension() != kdim),
			"All input curves must be 3-dimensional!");

	    //qc[ki] =
	    // newCurve(qpt->in,qpt->ik,qpt->et,qpt->ecoef,qpt->ikind,kdim,1);
	    qc_bnd[ki] = shared_ptr<ParamCurve>
		(new SplineCurve(qpt_bnd->numCoefs(), qpt_bnd->order(),
				   qpt_bnd->basis().begin(),
				   qpt_bnd->coefs_begin(), kdim));
	    qc_cross[ki] = shared_ptr<ParamCurve>
		(new SplineCurve(qpt_cross->numCoefs(), qpt_cross->order(),
				   qpt_cross->basis().begin(),
				   qpt_cross->coefs_begin(), kdim));

	    //	    if (qc[ki] == SISL_NULL) goto err101;
	}
    
    /* Turn the orientation of the curves at the 3. and 4. edge. */
    /* Compute Coon's patch.  */
    // As CoonsPatchGen assumes curves form a loop (when traversing from start to end),
    // we do not have to turn the curves.
    qsurf = CoonsPatchGen::createCoonsPatch(qc_bnd, qc_cross,
					      neighbour_tol, kink_tol);

//     // For debugging
//     std::ofstream dump("data/surf_debugdump.g2");
//     qsurf->writeStandardHeader(dump);
//     qsurf->write(dump);
//     // end of debugging

//     if (kstat < 0) goto error;
   
    /* Evaluate the surface in the middle.  */

    spar[0] = 0.5*(qsurf->startparam_u() + qsurf->endparam_u());
    spar[1] = 0.5*(qsurf->startparam_v() + qsurf->endparam_v());

    vector<Point> pts((kder + 1)*(kder + 2)/2);
    qsurf->point(pts, spar[0], spar[1], kder);


    // We must transfer the values to sder:
    for (size_t i = 0; i < pts.size(); ++i)
	for (int j = 0; j < pts[i].size(); ++j)
	    sder[i*kdim+j] = pts[i][j];
  
    /* Compute tangent vectors at the midpoint of the vertex region as 
       half the tangent vectors at the rectangular patch.  */

    for (ki=0; ki<kdim; ki++)
	{
// 	    etang[ki] = (double)0.5*sder[2*kdim+ki];
// 	    etang[kdim+ki] = (double)0.5*sder[kdim+ki];
// 	    etang[2*kdim+ki] = -(double)0.5*sder[2*kdim+ki];
// 	    etang[3*kdim+ki] = -(double)0.5*sder[kdim+ki];

	    etang[ki] = sder[kdim+ki];
	    etang[kdim+ki] = sder[2*kdim+ki];
	    etang[2*kdim+ki] = -sder[kdim+ki];
	    etang[3*kdim+ki] = -sder[2*kdim+ki];
	}

  
//     /* Application driven routine to alter the midpoint and tangents in the
//        midpoint.  */

//     fshape(sder,etang,kdim,icurv,&kstat);
//     if (kstat < 0) goto error;
  
    /* Copy value and 1. derivatives of first patch.  */

//     memcopy(eder,sder,kdim,DOUBLE);
    copy(sder.begin(), sder.begin() + kdim, eder.begin());
//     memcopy(eder+kdim,etang+3*kdim,kdim,DOUBLE);
    copy(etang.begin() +3*kdim, etang.begin() + 4*kdim, eder.begin() + kdim);
//     memcopy(eder+2*kdim,etang,kdim,DOUBLE);
    copy(etang.begin(), etang.begin() + kdim, eder.begin() + 2*kdim);
  
    /* Compute 2. derivatives.  */

    for (ki=0; ki<kdim; ki++)
	{
	    eder[3*kdim+ki] = sder[5*kdim+ki];
	    eder[4*kdim+ki] = sder[4*kdim+ki];
	    eder[5*kdim+ki] = sder[3*kdim+ki];
	}
  
//     *jstat = 0;
//     goto out;
  
//     /* Error in scratch allocation. */

//     err101 :
// 	*jstat = -101;
//     goto out;
 
//     /* Error in input. Dimension not equal to 3.  */

//     err104 :
// 	*jstat = -104;
//     goto out;
   
//     /* Error in a lower level function.  */

//  error:
//     *jstat = kstat;
//     goto out;
  
//     out :

// 	/* Free space occupied by local curves and surface. */

// 	for (ki=0; ki<8; ki++)
// 	    if (qc[ki] != SISL_NULL) freeCurve(qc[ki]);
//     if (qsurf != SISL_NULL) freeSurf(qsurf);
  
    return;
}



// void
//   sh1464(fshapeProc fshape,SISLCurve *vboundc[],int icurv,
// 	 double etwist[],double etang[],double eder[],int *jstat)
void
midpoint5(vector<shared_ptr<SplineCurve> >& curves, int icurv,
			     vector<double>& etwist, vector<double>& etang,
			     vector<double>& eder)
    /*
*********************************************************************
*                                                                   
* PURPOSE    : Given a five sided vertex region, evaluate the first
*              blending surface in the corner lying in the middle of the
*              vertex region. Compute the tangent vectors in the middle 
*              vertex along the inner boundaries of the region.
*
*
*
* INPUT      : fshape  - Application driven routine that gives the user an
*                        ability to change the middle point of the region
*                        (the vertex at which the blending surfaces meet),
*                        and the tangent vectors in the middle point along
*                        the curves which divedes the region. 
*              vboundc - Position and cross-tangent curves around the vertex
*                        region. For each edge of the region position and cross-
*                        tangent curves are given. The curves follow each other
*                        around the region and are oriented counter-clock-wise.
*                        The dimension of the array is 10.
*              icurv   - Number of sides. icurv = 5.
*              etwist  - Twist-vectors of the corners of the vertex region. The
*                        first element of the array is the twist in the corner
*                        before the first edge, etc. The dimension of the array
*                        is icurve*kdim.
*                       
*
* OUTPUT     : etang   - Tangent vectors at the midpoint of the vertex region.
*                        The dimension is icurv*idim.
*              eder    - Value, first and second derivative of the first blending
*                        surface in the corner at the midpoint. The sequence is the
*                        following : Value, 1. derivative in 1. parameter direction,
*                        1. derivative in the 2. parameter direction, 2. derivative
*                        in the 1. parameter direction, mixed derivative and 2.
*                        derivative in the 2. parameter direction. Dimension 6*idim.
*              jstat   - status messages  
*                                         > 0      : warning
*                                         = 0      : ok
*                                         < 0      : error
*
*
* METHOD     : Evaluate the Gregory Charrot function in the midpoint of the
*              vertex region. Compute the wanted derivatives.
*
* REFERENCES : 
*
* USE        : 3D geometry only.
*
*-
* CALLS      : sh1467  - Evaluate Gregory Charrot function over 
*			 5-sided region.                       
*
* WRITTEN BY : Vibeke Skytt, SI, 03.90.
*
*********************************************************************
*/
{
    const double PI = 3.14159265358979323;

    ALWAYS_ERROR_IF(icurv != 5,
		"Routine assumes polygonial surface is 5-sided!");

    int kdim = 3;          /* Dimension of geometry space.  */
    ALWAYS_ERROR_IF((int(etwist.size()) != icurv*kdim) || (int(etang.size()) != icurv*kdim) ||
		(int(eder.size()) < 6*kdim),
		"Wrong size of input vector!");

    int kder = 2;          /* Number of derviatives to evaluate.  */
    int ki;                /* Counter.  */
    double tlambda = 1.0/sqrt(5.0);           /* Constant.  */
    double tl1 = 2.0*tlambda*tan(PI/5.0);     /* Constant.  */
    double tl2 = sin(0.3*PI);                         /* Constant.  */
    double tconst1 = tl1/2.0 - tlambda;               /* Constant.  */
    double tconst2 = tl2 - tlambda;                           /* Constant.  */
    vector<double> sbar(5);  /* Barycentric coordinates of the blending function. */
    vector<double> sder(18);  /* Value and derivatives of blending function.       */

    /* Set up the barycentric coordinates of the midpoint of the region. */

    sbar[0] = sbar[1] = sbar[2] = sbar[3] = sbar[4] = tlambda;
  
    /* Evaluate the Gregory Charrot function at the midpoint. */
    //    sh1467(curves,etwist,kder,sbar,sder,&kstat);
    gregoryCharrotFunction5(curves,etwist, kder, sbar, sder);
    //    if (kstat < 0) goto error;
   
    /* Compute tangent vectors.  */

    for (ki=0; ki<kdim; ki++)
	{
	    etang[ki] = sder[kdim+ki]*tconst1 + sder[2*kdim+ki]*tconst2;
	    //	    etang[ki] *= 0.1;
	    etang[kdim+ki] = -sder[kdim+ki]*tlambda + sder[2*kdim+ki]*tconst1;
	    //	    etang[kdim+ki] *= 0.1;
	    etang[2*kdim+ki] = sder[kdim+ki]*tconst1 - sder[2*kdim+ki]*tlambda;
	    //	    etang[2*kdim+ki] *= 0.1;
	    etang[3*kdim+ki] = sder[kdim+ki]*tconst2 + sder[2*kdim+ki]*tconst1;
	    //	    etang[3*kdim+ki] *= 0.1;
	    etang[4*kdim+ki] = sder[kdim+ki]*tconst2 + sder[2*kdim+ki]*tconst2;
	    //	    etang[4*kdim+ki] *= 0.1;
	}
  
//     /* Application driven routine to alter the midpoint and tangents in the
//        midpoint.  */

//     fshape(sder,etang,kdim,icurv,&kstat);
//     if (kstat < 0) goto error;
  
    /* Copy value and 1. derivatives of first patch.  */

//     memcopy(eder,sder,kdim,DOUBLE);
    copy(sder.begin(), sder.begin() + kdim, eder.begin());
//     memcopy(eder+kdim,etang+4*kdim,kdim,DOUBLE);
    copy(etang.begin() + 4*kdim, etang.begin() + 5*kdim, eder.begin()+ kdim);
//     memcopy(eder+2*kdim,etang,kdim,DOUBLE);
    copy(etang.begin(), etang.begin() + kdim, eder.begin() + 2*kdim);

    /* Compute 2. derivatives.  */

    for (ki=0; ki<kdim; ki++)
	{
	    eder[3*kdim+ki] = sder[3*kdim+ki]*tconst2*tconst2 
		+ 2.0*sder[4*kdim+ki]*tconst2*tconst2 
		+ sder[5*kdim+ki]*tconst2*tconst2;
	    //	    eder[3*kdim+ki] *= 0.1;
	    eder[4*kdim+ki] = sder[3*kdim+ki]*tconst1*tconst2 
		+ sder[4*kdim+ki]*tconst2*(tconst1+tconst2)
		+ sder[5*kdim+ki]*tconst2*tconst2;
	    //	    eder[4*kdim+ki] *= 0.1;
	    eder[5*kdim+ki] = sder[3*kdim+ki]*tconst1*tconst1 
		- 2.0*sder[4*kdim+ki]*tconst1*tconst2 
		+ sder[5*kdim+ki]*tconst2*tconst2;
	    //	    eder[5*kdim+ki] *= 0.1;
	}
  
//     *jstat = 0;
//     goto out;
  
//     /* Error in a lower level function.  */

//  error:
//     *jstat = kstat;
//     goto out;
  
//     out :
	return;
}



// void
//   sh1465(fshapeProc fshape,SISLCurve *vboundc[],int icurv,
// 	 double etwist[],double etang[],double eder[],int *jstat)
void
midpoint6(vector<shared_ptr<SplineCurve> >& curves, int icurv,
			     vector<double>& etwist, vector<double>& etang,
			     vector<double>& eder)
    /*
*********************************************************************
*                                                                   
* PURPOSE    : Given a vertex region whith an equal number of sides, 
*              evaluate the first blending surface in the corner lying in 
*              the middle of the vertex region. Compute the tangent vectors 
*              in the middle vertex along the inner boundaries of the region.
*
*
*
* INPUT      : fshape  - Application driven routine that gives the user an
*                        ability to change the middle point of the region
*                        (the vertex at which the blending surfaces meet),
*                        and the tangent vectors in the middle point along
*                        the curves which divedes the region. 
*              vboundc - Position and cross-tangent curves around the vertex
*                        region. For each edge of the region position and cross-
*                        tangent curves are given. The curves follow each other
*                        around the region and are oriented counter-clock-wise.
*                        The dimension of the array is 2*icurv.
*              icurv   - Number of sides.
*              etwist  - Twist-vectors of the corners of the vertex region. The
*                        first element of the array is the twist in the corner
*                        before the first edge, etc. The dimension of the array
*                        is icurv*kdim.
*                       
*
* OUTPUT     : etang   - Tangent vectors at the midpoint of the vertex region.
*                        The dimension is icurv*idim.
*              eder    - Value, first and second derivative of the first blending
*                        surface in the corner at the midpoint. The sequence is the
*                        following : Value, 1. derivative in 1. parameter direction,
*                        1. derivative in the 2. parameter direction, 2. derivative
*                        in the 1. parameter direction, mixed derivative and 2.
*                        derivative in the 2. parameter direction. Dimension 6*idim.
*              jstat   - status messages  
*                                         > 0      : warning
*                                         = 0      : ok
*                                         < 0      : error
*
*
* METHOD     : 
*
* REFERENCES : 
*
* USE        : 3D geometry only.
*
*-
* CALLS      : s6norm   - Normalize vector.  
*              s6scpr   - Scalar product between two vectors.  
*              s6length - Lenght of vector.  
*              s6dist   - Distance between two vectors. 
*              s6crss   - Cross product between two vectors. 
*              s6diff   - Difference between two vectors.  
*              s6curvrad - Estimate curvature radius in point. 
*              s1325    - Calculate tangent length given opening
*                         angle of curve segment and curvature radius.
*              s1221    - Evaluate curve.  
*              sh1465_s9der2 - Compute 2. derivatives of the first 
*			       blending patch.    
*              
*
* WRITTEN BY : Vibeke Skytt, SI, 04.90.
*
*********************************************************************
*/
{
    // @@sbr Hardcoded values! Not sure what values are appropriate.
    const double ANGULAR_TOLERANCE = 0.1; // In radians.
    const double equal_tol = 1e-15; // Tolerance for equalty of values.

    ALWAYS_ERROR_IF(icurv != 6,
		"Routine assumes polygonial surface is 6-sided!");

    int kdim = 3;          /* Dimension of geometry space.  */
    ALWAYS_ERROR_IF((int(etwist.size()) != icurv*kdim) || (int(etang.size()) != icurv*kdim) ||
		(int(eder.size()) < 6*kdim),
		"Wrong size of input vector!");

    int kder = 0;       /* Number of derivatives of curve to evaluate. */
    int kder1 = 1;      /* Number of derivatives of curve to evaluate. */
    int ki,kj;          /* Counters.  */
    int kcurv2;         /* Number of edges diveded in two.  */
    double tpar;        /* Parameter value at which to evaluate curve. */
    double t1;          /* Scalar product between tangent in midpoint
			   and normal of vertex region in midpoint.    */
    double trad1,trad2; /* Estimate of curvature radius of curves in
			   endpoints.   */
    double ta1;         /* Opening angle of curve segment.  */
    double tb1,tb2;     /* Length of original tangents in endpoints
			   of curve segments.               */
    double tang1,tang2; /* Length of tangents based on curvature radius. */
    double tscal1;      /* Scaler product between unit tangent vectors.  */
    double tdist;       /* Distance between endpoints of curve segment.  */
    Point smid(3);     /* Midpoint of vertex region.       */
    Point snorm1(3);   /* Cross product of two tangents in the midpoint. */
    Point snorm(3);    /* Normal of vertex region in the midpoint.      */
    vector<Point> svec; /* Tangent along curve in the midpoint of the two
			   first position curves.           */
    for (ki = 0; ki < 2; ++ki)
	svec.push_back(Point(kdim));
    //double *sder = NULL;//SISL_NULL;/* Value of boundary curves at the midpoints of
    //				   the curves.                   */
    vector<Point> sder;
    for (ki = 0; ki < 2*icurv; ++ki)
	sder.push_back(Point(kdim));
    //    double *stang = NULL; //SISL_NULL;  /* Tangent vectors in the
                                        // midpoint of the region. */
    vector<Point> stang;
    for (ki = 0; ki < icurv; ++ki)
	stang.push_back(Point(kdim));
    shared_ptr<SplineCurve> qc;      /* Local pointer to edge curve.  */
  
    kcurv2 = icurv/2;

    /* Allocate scratch for values on edge curves. */

//     if ((sder = newarray(2*icurv*kdim,DOUBLE)) == SISL_NULL) goto err101;
//     if ((stang = newarray(icurv*kdim,DOUBLE)) == SISL_NULL) goto err101;

    for (ki=0; ki<icurv; ki++)
	{  
	    qc = curves[2*ki];

	    /* Test dimension of curve. */

	    //	    if (qc->idim != kdim) goto err102;
	    ALWAYS_ERROR_IF(qc->dimension() != kdim,
			"Input curves must be 3-dimensional.");

	    /* Evaluate boundary curves of current edge in the midpoint.
	       First find midpoint.  */

	    //	    tpar = (*(qc->et + qc->ik - 1) + *(qc->et + qc->in))/2.0;
      	    tpar = (*(qc->basis().begin() + qc->order() - 1) +
		    *(qc->basis().begin() + qc->numCoefs()))/2.0;

	    /* Evaluate position curve. */
	    vector<Point> pts(kder1+1);
	    qc->point(pts, tpar, kder1);
	    // We must transfer the values to sder+2*ki*kdim:
	    for (size_t i = 0; i < pts.size(); ++i)
		sder[2*ki+i] = pts[i];
      
	    if (ki < 2)
		svec[ki] = sder[2*ki+1];
	    if (ki == 1)
		kder1 = 0;
      
	    /* Evaluate cross derivative curve. */
	    pts.resize(kder+1);
	    curves[2*ki+1]->point(pts, tpar, kder);
	    // We must transfer the values to sder+2*ki*kdim:
	    for (size_t i = 0; i < pts.size(); ++i)
		sder[2*ki+1+i] = pts[i];
	    //	    if (kstat < 0) goto error;
	}
  
    for (ki=0; ki<kcurv2; ki++)
	{
	    /* Set new tangent length based on curvature radius. First estimate
	       curvature radius.  */
//       s6curvrad(sder+2*ki*kdim,sder+2*(ki+kcurv2)*kdim,sder+(2*ki+1)*kdim,
// 		kdim,&trad1,&kstat);
	    curvrad(sder[2*ki], sder[2*(ki+kcurv2)], sder[2*ki+1], kdim, &trad1);
//       s6curvrad(sder+2*ki*kdim,sder+2*(ki+kcurv2)*kdim,
//              sder+(2*(ki+kcurv2)+1)*kdim,kdim,&trad2,&kstat);
	    curvrad(sder[2*ki], sder[2*(ki+kcurv2)], sder[(2*(ki+kcurv2)+1)],
		    kdim, &trad2);
  
	    /* Normalize tangents. */

// 	    tb1 = s6norm(sder+(2*ki+1)*kdim,kdim,sder+(2*ki+1)*kdim,&kstat);
	    sder[2*ki+1].normalize();
	    tb1 = sder[2*ki+1].length(); // @@sbr 191202 Hmm... Length prior to normalization?
// 	    tb2 = s6norm(sder+(2*(ki+kcurv2)+1)*kdim,kdim,
// 			 sder+(2*(ki+kcurv2)+1)*kdim,&kstat);
	    sder[2*(ki+kcurv2)+1].normalize();
	    tb2 = sder[2*(ki+kcurv2)+1].length();

	    /* Compute distance between endpoints of curve.  */

	    //	    tdist = s6dist(sder+2*ki*kdim,sder+2*(ki+kcurv2)*kdim,kdim);
      	    tdist = sder[2*ki].dist(sder[2*(ki+kcurv2)]);

	    /* Find opening angle of curve segment.  */

	    // tscal1  = s6scpr(sder+(2*ki+1)*kdim,sder+(2*(ki+kcurv2)+1)*kdim,kdim);
	    tscal1 = sder[2*ki+1]*sder[2*(ki+kcurv2)+1];


	    if (tscal1 >= 0.0)
		tscal1  = min(1.0,tscal1);
	    else
		tscal1  = max(-1.0,tscal1);

	    ta1 = acos(tscal1);

	    if (fabs(ta1) < ANGULAR_TOLERANCE) ta1 = 0.0;

	    //	    if (DNEQUAL(ta1,0.0))
	    if (fabs(ta1) > equal_tol)
		{
		    /*  Make tangents based on radius of curvature */

		    tang1 = local_s1325(trad1,ta1);
		    tang2 = local_s1325(trad2,ta1);
		}

	    /* Test if the found tangent length can be used. Otherwise
	       adjust the length.   */

	    //	    if (DEQUAL(ta1,0.0) || trad1 < 0)
	    if ((fabs(ta1) < equal_tol) || trad1 < 0)
		tang1 = tdist/3.0;
	    //	    if (DEQUAL(ta1,0.0) || trad2 < 0)
	    if ((fabs(ta1) < equal_tol) || trad2 < 0)
		tang2 = tdist/3.0;
	    if (tang1 > 0.5*tdist || tang2 > 0.5*tdist) 
		{
		    tang1 = tb1;
		    tang2 = tb2;
		}
      
	    /* Set tangent length of tangents in endpoints of the curves to find. */

	    for (kj=0; kj<kdim; kj++)
		{
		    sder[(2*ki+1)][kj] *= tang1;
		    sder[(2*(ki+kcurv2)+1)][kj] *= tang2;
		} 
	}  
  

    /* Estimate midpoint of region and tangent in midpoint.  */

    for (kj=0; kj<kdim; kj++) 
	{
	    smid[kj] = 0.0;
	    snorm[kj] = 0.0;
      
	    for (ki=0; ki<kcurv2; ki++)
		{
		    /* Estimate midpoint.  */

		    smid[kj] +=
			(sder[2*ki][kj] + sder[2*(ki+kcurv2)][kj])/2.0 +
			(sder[2*ki+1][kj] +
			 sder[2*(ki+kcurv2)+1][kj])/8.0;
	  
		    /* Compute tangent.  */

		    stang[(ki+kcurv2)][kj] = 1.5*(sder[2*(ki+kcurv2)][kj]
						      - sder[2*ki][kj])
			+ 0.25*(sder[2*(ki+kcurv2)+1][kj]
				- sder[(2*ki+1)][kj]);
		    stang[ki][kj] = - stang[(ki+kcurv2)][kj];
		}
	    smid[kj] /= kcurv2;
	}

    /* Find medium plane given by the tangents.  */

    for (ki=0; ki<kcurv2; ki++)
	{
	    //	    s6crss(stang+ki*kdim,stang+(ki+1)*kdim,snorm1);
	    snorm1 = stang[ki]%stang[ki+1];
  
	    for (kj=0; kj<kdim; kj++) 
		snorm[kj] += snorm1[kj]/kcurv2;
	}
  
    //    (void)s6norm(snorm,kdim,snorm,&kstat);
    snorm.normalize();

    /* Project the tangents into this plane.  */

    for (ki=0; ki<icurv; ki++)
	{
	    //	    t1 = -s6scpr(stang+ki*kdim,snorm,kdim);
	    t1 = -stang[ki]*snorm;

	    for (kj=0; kj<kdim; kj++)
		stang[ki][kj] += t1*snorm[kj];
	}
  
//     /* Application driven routine to alter the midpoint and tangents in the
//        midpoint.  */

//     fshape(smid,stang,kdim,icurv,&kstat);
//     if (kstat < 0) goto error;
  
    /* Compute second order derivatives of first surface in the midpoint.  */

    //    sh1465_s9der2(sder,smid,stang,snorm,svec,icurv,kdim,eder+3*kdim,&kstat);
    vector<double> eder_part(eder.begin() + 3*kdim, eder.begin() + 6*kdim);
    midpoint6_s9der2(sder, smid, stang, snorm, svec, icurv, kdim, eder_part);
    copy(eder_part.begin(), eder_part.begin() + 3*kdim, eder.begin() + 3*kdim);

    //    if (kstat < 0) goto error;
  
    /* Set lengths of tangents and second derivatives according to patches
       whith side length equal to half the curvelengths considered in this
       routine.  */

    for (ki=0; ki<icurv; ki++)
	for (kj = 0; kj < kdim; ++kj)
	    stang[ki][kj] *= 0.5;
    for (ki=3*kdim; ki<6*kdim; ki++)
	eder[ki] *= 0.25;
  
    /* Copy position and tangent information into the output array giving
       derivatives of the first blending surface in the midpoint. */

//     memcopy(eder,smid,kdim,DOUBLE);
    copy(smid.begin(), smid.begin() + kdim, eder.begin());
//     memcopy(eder+kdim,stang,2*kdim,DOUBLE);
    copy(stang[0].begin(), stang[0].begin() + kdim, eder.begin() + kdim);
    copy(stang[1].begin(), stang[1].begin() + kdim, eder.begin() + 2*kdim);
  
    /* Copy tangents into output array containing tangents.  */

//     memcopy(etang,stang+kdim,(icurv-1)*kdim,DOUBLE);
    for (ki = 0; ki < (icurv-1); ++ki)
	copy(stang[1+ki].begin(), stang[1+ki].begin() + kdim,
	     etang.begin() + ki*kdim);
//     memcopy(etang+(icurv-1)*kdim,stang,kdim,DOUBLE);
    copy(stang[0].begin(), stang[0].begin() + kdim, etang.begin() + (icurv-1)*kdim);
					
//     *jstat = 0;
//     goto out;
  
//     /* Error in scratch allocation.  */

//     err101 :
// 	*jstat = -101;
//     goto out;
  
//     /* Error in input. Dimension not equal to 3.  */

//     err102 :
// 	*jstat = -102;
//     goto out;
  
//     /* Error in a lower level function.  */

//  error:
//     *jstat = kstat;
//     goto out;
  
//     out :

// 	/* Free space occupied by local arrays.  */

// 	if (sder) freearray(sder);
//     if (stang) freearray(stang);
  
    return;
}



// static void
//   sh1465_s9der2(double ebound[],double epoint[],double etang[],
// 		double enorm[],double evec[],int icurv,
// 		int idim,double eder2[],int *jstat)
void
midpoint6_s9der2(vector<Point>& ebound, Point& epoint,
				    vector<Point>& etang, Point& enorm,
				    vector<Point>& evec,
				    int icurv, int idim, vector<double>& eder2)
    /*
*********************************************************************
*                                                                   
* PURPOSE    : Given a vertex region with an equal number of sides, 
*              estimate the 2. derivatives of the first blending surface 
*              in the midpoint of the region.
*
*
*
* INPUT      : ebound  - Position of boundary curve and cross tangent at the 
*                        midpoint of each edge. Dimension is 2*icurv*idim.
*              epoint  - The midpoint of the vertex region. Dimension is idim.
*              etang   - Tangents in epoint, pointing towards the midpoints of
*                        the edges. Dimension is icurv*idim.
*              enorm   - Normal to surface in midpoint of region. Dimension is idim.
*              evec    - The tangent at the midpoint of the boundary curves at the
*                        two first edges. Dimension is 2*idim.
*              icurv   - Number of sides. icurv is an equal number.
*              idim    - Dimension of the geometry space.
*                       
*
* OUTPUT     : eder2   - Second derivative of the first blending surface in 
*                        the corner at the midpoint. The sequence is the
*                        following : 2. derivative in the 1. parameter direction, 
*                        mixed derivative and 2. derivative in the 2. parameter 
*                        direction. Dimension is 3*idim.
*              jstat   - status messages  
*                                         > 0      : warning
*                                         = 0      : ok
*                                         < 0      : error
*
*
* METHOD     : 
*
* REFERENCES : 
*
* USE        : 3D geometry only.
*
*-
* CALLS      : s1221     - Curve evaluator.  
*              s1334     - Curve interpolation. 
*              s6lufacp  - LU-factorizing of matrix. 
*              s6lusolp  - Solve to LU-factorized equation system. 
*              s6curvature - Compute curvature vector of curve. 
*              s6scpr    - Scalar product between two vectors. 
*              s6length  - Length of vector.  
*              
*
* WRITTEN BY : Vibeke Skytt, SI, 04.90.
*
*********************************************************************
*/
{
    // @@sbr Hardcoded value! Not sure what value is appropriate.
    const double equal_tol = 1e-15; // Tolerance for equalty of values.
    const double PI = 3.14159265358979323;

    int ki,kj;         /* Counters.         */
    int kcurv2;        /* Number of edges divided into 2. */
    int kord = 6;      /* Order of interpolated curve.    */
    vector<int> ll(3);         /* Pivoting array in solving equation system. */
    double tpar;       /* Parameter value at which to evaluate. */
    double tncurv;     /* Normal curvature at midpoint of curve. */
    double te,tf,tg;   /* Coefficients of first fundamental form of surface. */
    double tl,tm,tn;   /* Coefficients of second fundamental form of surface. */
    double talfa;      /* Angle between the two parameter directions of a
			  rectangular patch in the parameter area of the 
			  n-sided vertex region.       */
    double tcos;       /* Cosinus of the angle talfa.  */
    double tdudt;      /* Factor in parameter change.  */
    double tdvdt;      /* Factor in parameter change.  */
    double tform1;     /* First fundamental form.      */
    //    double spoint[18]; /* Interpolation conditions.    */
    vector<double> spoint(9);
    vector<double> spoint_der(9);
    vector<double> stype(6);   /* Type of interpolation conditions. */
    // double *spar 0; // = SISL_NULL; Parameter value of interpolation conditions.
    vector<double> spar;
    vector<Point> sder;   /* Value and derivatives of curve in the midpoint. */
    for (ki = 0; ki < 6; ++ki)
	sder.push_back(Point(3));
    //    double scurv[3];   /* Curvature vector in midpoint of curve. */
    Point scurv;
    vector<double> smat(9); // Matrix in equation system to compute mixed
			    // derivative.
    SplineCurve *qc; // = SISL_NULL;  /* Curve across the vertex region
	                   // between midpoint of edge position curves through
	                   // the midpoint of the region.      */
  
    kcurv2 = icurv/2;
  
    /* Test input.  */

    ALWAYS_ERROR_IF(idim != 3,
		"Input dimension must be 3.");
  
    /* Make curves that limits first blending surface. */

    for (ki=0; ki<2; ki++)
	{
	    /* Set up interpolation conditions of curve. */

	    for (kj=0; kj<idim; kj++)
		{
		    spoint[kj] = ebound[2*(kcurv2+ki)][kj];
		    spoint_der[kj] = ebound[(2*(kcurv2+ki)+1)][kj];
		    spoint[idim+kj] = epoint[kj];
		    spoint_der[idim+kj] = etang[ki][kj];
		    spoint[2*idim+kj] = ebound[2*ki][kj];
		    spoint_der[2*idim+kj] = -ebound[(2*ki+1)][kj];
		}
      
	    /* Type of interpolation conditions.  */

	    stype[0] = 1.0;
	    stype[1] = 4.0;
	    stype[2] = 1.0;
	    stype[3] = 4.0;
	    stype[4] = 1.0;
	    stype[5] = 4.0;

	    /* Interpolate curve.  */
	    SplineInterpolator interpolator;
	    vector<double> params;
	    params.push_back(0.0);
	    params.push_back(0.5);
	    params.push_back(1.0);
	    vector<int> tangent_index;
	    tangent_index.push_back(1);
	    tangent_index.push_back(3);
	    tangent_index.push_back(5);
	    vector<double> coefs;
	    interpolator.interpolate(params, spoint, tangent_index,
				     spoint_der, kord, coefs);

	    int nmb = (int)coefs.size() / idim;
	    spar.insert(spar.end(),
			interpolator.basis().begin(), interpolator.basis().end());
	    qc = new SplineCurve(nmb, kord, spar.begin(), coefs.begin(), idim);

	    /* Evaluate curve in midpoint.  */
				
	    tpar = spar[1];
	    qc->point(sder, tpar, 2);
	}
  
    /* Copy 2. derivatives of the two curves to the 2. derivatives of the
       surface in the midpoint in the 1. and 2. parameter direction. */

//     memcopy(eder2,sder+2*idim,idim,DOUBLE);
    copy(sder[2].begin(), sder[2].begin() + idim, eder2.begin());
//     memcopy(eder2+2*idim,sder+5*idim,idim,DOUBLE);
    copy(sder[5].begin(), sder[5].begin() + idim, eder2.begin() + 2*idim);
  
    /* Compute curvature vector in the midpoint of the first curve.  */

    //    s6curvature(sder,idim,scurv,&kstat);
    vector<Point> temp_vec(sder.begin(), sder.begin() + idim);
    curvature(temp_vec, idim, scurv);

    /* Compute normal curvature of the first curve at the midpoint of the region. */
  
    //    tncurv = s6scpr(scurv,enorm,idim);
    tncurv = scurv*enorm;

    /* Compute parameter direction of the curve compared to that of the
       first blending surface. */

    //    talfa = TWOPI/icurv;
    talfa = 2*PI/icurv;
    tcos = cos(talfa);
    //    tdudt = (DEQUAL(tcos+1.0,1.0)) ? 0.0 : 1.0/tcos;
    tdudt = (fabs(tcos) < equal_tol) ? 0.0 : 1.0/tcos;
    tdvdt = 1.0;
  
    /* Compute coefficients of the first fundamental form of the surface. */

//     te = s6scpr(sder+idim,sder+idim,idim);
    te = sder[1]*sder[1];
//     tf = s6scpr(sder+idim,sder+4*idim,idim);
    tf = sder[1]*sder[4];
//     tg = s6scpr(sder+4*idim,sder+4*idim,idim);
    tg = sder[4]*sder[4];
  
    /* Compute the first fundamental form.  */

    tform1 = te*tdudt*tdudt + 2.0*tf*tdudt*tdvdt + tg*tdvdt*tdvdt;
  
    /* Compute 1. and 3. coefficient of the second fundamental form of the surface. */

//     tl = s6scpr(sder+2*idim,enorm,idim);
    tl = sder[2]*enorm;
//     tn = s6scpr(sder+5*idim,enorm,idim);
    tn = sder[5]*enorm;
  
    /* Compute 2. coefficient of the second fundamental form which is equal to
       the length of the component of the twist along the surface normal.  */

    eder2[idim] = tm = (tncurv*tform1 - tl*tdudt*tdudt - tn*tdvdt*tdvdt)/2.0;
  
    /* Set the length of the component of the twist along the derivative in
       the first parameter direction equal to zero.   */

    eder2[idim+1] = 0.0;
  
    /* Set the length of the component of the twist along the derivative 
       in the second parameter direction equal to zero.   */

    eder2[idim+2] = 0.0;
  
    /* Compute twist vector at the midpoint.  */

//     memcopy(smat,enorm,idim,DOUBLE);
    copy(enorm.begin(), enorm.end(), smat.begin());
//     memcopy(smat+idim,sder+idim,idim,DOUBLE);
    copy(sder[1].begin(), sder[1].end(), smat.begin() + idim);
//     memcopy(smat+2*idim,sder+4*idim,idim,DOUBLE);
    copy(sder[4].begin(), sder[4].end(), smat.begin() + 2*idim);


    vector<vector<double> > smatrix(3);
    for (int ii = 0; ii < 3; smatrix[ii++].resize(3)) {}
    for (ki = 0; ki < idim; ++ki) {
	for (kj = 0; kj < idim; ++kj) {
	    smatrix[ki][kj] = smat[ki*idim + kj];
	}
    }
    vector<double> sol(3);
    for (ki = 0; ki < idim; ++ki) {
	sol[ki] = eder2[idim + ki];
    }
    LUsolveSystem(smatrix, (int)smatrix.size(), &sol[0]);
    for (ki = 0; ki < idim; ++ki) {
	eder2[idim + ki] = sol[ki];
    }

// -- the below block of code is 'newmat' dependent.  Replaced by above eq. ---
//     Matrix smatrix(3, 3);
//     for (ki = 0; ki < idim; ++ki)
// 	for (kj = 0; kj < idim; ++kj)
// 	    smatrix.element(ki, kj) = smat[ki*idim+kj];

//     //    s6lufacp(smat,ll,3,&kstat);
//     //    if (kstat < 0) goto error;
//     CroutMatrix smatrixLUfact = smatrix;
//     ALWAYS_ERROR_IF(smatrixLUfact.IsSingular(),
// 		"System couldn't be solved; should not happen!");

// 	    //    s6lusolp(smat,eder2+idim,ll,3,&kstat);
//     ColumnVector unknowns(3);
//     for (ki = 0; ki < idim; ++ki)
// 	unknowns.element(ki) = eder2[idim + ki];
//     ColumnVector solution = smatrixLUfact.i() * unknowns;
//     for (ki = 0; ki < idim; ++ki)
// 	eder2[idim + ki] = solution.element(ki);


    //    if (kstat < 0) goto error;
  
//     *jstat = 0;
//     goto out;
  
//     /* Error in input. Dimension not equal to 3.  */

//     err104 :
// 	*jstat = -104;
//     goto out;
  
//     /* Error in lower level routine.  */

//     error :
// 	*jstat = kstat;
//     goto out;
  
//     out :
	return;
}     



// void
//       s6curvrad(double epnt1[],double epnt2[],double etang[],int idim,
// 		double *crad,int *jstat)
void curvrad(const Point& epnt1, const Point& epnt2,
				const Point& etang, int idim, double *crad)
/*
*********************************************************************
*                                                                   
* PURPOSE    : Given both endpoints of a curve segment and the tangent in
*              one of the endpoints, estimate the curvature radius of the 
*              curve segment in the endpoint where the tangent is given.
*
*
* INPUT      : epnt1   - First endpoint.
*              epnt2   - Second endpoint.
*              etang   - Given tangent.
*              idim    - Dimension of geometry space.
*                       
*
* OUTPUT     : crad    - Estimated curvature radius.
*              jstat   - status messages  
*                                         = 1      : Tangent vector zero.
*                                         = 0      : ok
*                                         < 0      : error
*
*
* METHOD     : Pass a circle through the points having the given tangent
*              in the given endpoint. Compute the radius of this circle.
*
* REFERENCES : 
*
* USE        : 
*
*-
* CALLS      : s6dist - Distance between two points.  
*              s6diff - Difference vector between two vectors. 
*
* WRITTEN BY : Vibeke Skytt, SI, 05.90.
*
*********************************************************************
*/
{
    double REL_COMP_RES = 0.000000000000001;

    double tdist;         /* Distance between the endpoints.   */
    double tdot;          /* Scalar product between tangent and vector
                           between endpoints.                */
    double tlmid;         /* Length of tangent vector.         */
    double tcos;          /* Cosinus to the angle between the tangent
                           and the vector between the endpoints. */
    double tang;          /* Angle at the centre of the circle
			     between the vectors from origo to
			     the given endpoints.              */
    double tdum;          /* Denominator in expression to find 
			     the  radius of the circle.        */
    double trad;          /* The curvature radius, i.e. the radius
			     of the circle.                    */
    //  double sdiff[3];      /* Difference vector between the endpoints. */
    Point sdiff(3);

    /* Test input.  */

    ALWAYS_ERROR_IF(idim != 3,
		"Expecting 3-dimensional points.");
  
    /* Estimate curvature radius based on endpoints and tangent.  */

    //   tdist = s6dist(epnt1,epnt2,idim);
    tdist = epnt1.dist(epnt2);
    //   s6diff(epnt2,epnt1,idim,sdiff);
    sdiff = epnt2 - epnt1;

    //   tdot = s6scpr(etang,sdiff,idim);
    tdot = etang*sdiff;
    //   tlmid = s6length(etang,idim,&kstat);
    tlmid = etang.length();

    tcos = (tlmid*tdist != 0.0) ? fabs(tdot/(tlmid*tdist)) : fabs(tdot);
    tcos = min(1.0,tcos);
  
    tang = 2*acos(tcos);
    tdum = sqrt(2-2*cos(tang));
    trad = (tdum > REL_COMP_RES) ? tdist/tdum : -1.0;
  
    /* Set curvature radius. */
  
    *crad = trad;
    //   *jstat = 0;
    //   goto out;

    //   /* Error in input.  */

    //   err104 :
    //     *jstat = -104;
    //   goto out;

    //   out :
    return;
}    



// double 
// s1325(double aradiu,double angle)
double local_s1325(double aradiu,double angle)
/*
*********************************************************************
*                                                                   
* PURPOSE    : To create the tangent length for interpolating a
*              circular arc with an almost equi-oscillating Hermit qubic
*
* INPUT      : aradiu  - The radius of the circular arc
*              angle   - The opening angle of the circular arc
*
* OUTPUT     : s1325   - The proposed tangent length
*
* METHOD     : A second degree equation giving the tanget length is
*              solved
*
* WRITTEN BY : Tor Dokken, SI, Oslo, Norway. 30. June 1988
*                                  
*********************************************************************
*/
{
  double tcos,tsin;          /* Dummy variables                     */
  double ta,tb,tc,tl;        /* Dummy variables                     */
  double tconst = 1.85530139760811990992528773586425;
                             /* Constant used in the calculation    */
  
  tcos = cos(angle);
  tsin = sin(angle);
  
  /*  Calculate length of tangents
   *   tconst = (3-2sqrt(2))**1/3 + (3+2sqrt(2))**1/3 - 0.5 */
  
  ta     = (double)0.6*tconst - (double)0.9*tcos;
  tb     = ((double)0.4*tconst+(double)1.8)*tsin;
  tc     = ((double)0.4*tconst+(double)1.0)
           * tcos - (double)0.4*tconst - (double)1.0;
  tl     = aradiu*(-tb+sqrt(tb*tb-4*ta*tc))/((double)2.0*ta);

  return(tl);
}


// void
//       s6curvature(double eder[],int idim,double ecurv[],int *jstat)
void
curvature(const vector<Point>& eder, int idim, Point& ecurv)
/*
*********************************************************************
*                                                                   
* PURPOSE    : Given position, first and second derivative of a curve at
*              a point, compute curvature vector.
*
*
*
* INPUT      : eder    - Array containing position, 1. and 2. derivative of
*                        curve. Dimension is 3*idim.
*              idim    - Dimension of geometry space.
*                       
*
* OUTPUT     : ecurv   - Curvature vector.
*              jstat   - status messages  
*                                         = 1      : Tangent vector zero.
*                                         = 0      : ok
*                                         < 0      : error
*
*
* METHOD     : Express the curve using cord length parametrisation. Then the
*              curvature is equal to the 2. derivative of the curve.
*
* REFERENCES : 
*
* USE        : 
*
*-
* CALLS      : s6length - Length of vector.   
*              s6scpr   - Scalar product between two vectors.  
*
* WRITTEN BY : Vibeke Skytt, SI, 04.90.
*
*********************************************************************
*/
{
    ALWAYS_ERROR_IF(ecurv.dimension() != idim,
		    "Dimension of input point does not correspond with space.");


    int ki;           /* Counter.   */
    double tleng;     /* Length of first derivative.  */
    double tleng2;    /* Square of length.  */
    double tdot;      /* Scalar product between 1. and 2. derivative. */
  
    /* Compute length of 1. derivative. */

    //  tleng = s6length(eder+idim,idim,&kstat);
    tleng = eder[1].length();

    //   if (kstat == 0)
    //     {
    //       /* The first derivative is zero. */

    //       for (ki=0; ki<idim; ki++) ecurv[ki] = (double)0.0;
    //       goto warn1;
    //     }
    //   else
    {
	tleng2 = tleng*tleng;
	//      tdot = s6scpr(eder+idim,eder+2*idim,idim);
	tdot = eder[1]*eder[2];

	for (ki=0; ki<idim; ki++)
	    //    ecurv[ki] = (eder[2*idim+ki] - eder[idim*ki]*tdot/tleng2)/tleng2;
	    ecurv[ki] = (eder[2][ki] - eder[ki][0]*tdot/tleng2)/tleng2;
    }
}

} // end anonymous namespace
