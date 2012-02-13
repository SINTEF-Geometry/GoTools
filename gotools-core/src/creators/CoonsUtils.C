//===========================================================================
//                                                                           
// File: CoonsUtils.C                                                      
//                                                                           
// Created: Fri Aug  3 17:34:36 2001                                         
//                                                                           
// Author: Sverre Briseid <Sverre.Briseid@math.sintef.no>
//                                                                           
// Revision: $Id: CoonsUtils.C,v 1.5 2009-02-18 11:48:18 vsk Exp $
//                                                                           
// Description: Utility functions for modifying/adding/testing cross tangent
//              curves. Ported from sisl.                                           
//===========================================================================
//
// @@ Remove unused header files.
#include "GoTools/creators/CoonsPatchGen.h"
#include "GoTools/creators/CurveCreators.h"
#include "GoTools/geometry/SplineCurve.h"
#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/geometry/SplineInterpolator.h"
#include "GoTools/creators/HermiteAppC.h"
#include "GoTools/creators/CrossTanOffDist.h"
#include "GoTools/creators/CrossTangentOffset.h"
#include "GoTools/geometry/GeometryTools.h"
#include <algorithm>
#include <functional>
#include "GoTools/geometry/Utils.h"
//#include "newmat.h"
#include "GoTools/utils/LUDecomp.h"
#include <vector>
#include <numeric>
#include <fstream> // for debugging
#include <cstdio> // for debugging
#include <iterator>


using namespace Go;
using std::back_inserter;
using std::vector;
using std::swap;

namespace {

    /// Estimate the derivatives in one endpoint of the blending 
    /// functions meeting in a corner of a vertex region.
    /// \param ea matrix containing coefficients in equation
    ///           system representing conditions on the derivatives.
    /// \param eb rigth side of equation system.
    /// \param ix 
    /// \param ieq 
    /// \param irang 
    /// \param starti 
    /// \param aendi 
    /// \param astartj 
    /// \param aendj 
    /// \param ealfai 
    /// \param ebetai 
    /// \param nmbcoef1 
    /// \param nmbcoef2 
    /// \param ealfaj 
    /// \param ebetaj 
    /// \param nmbcoef2 
    void blendder(double ea[],double eb[],int ix,int ieq,
		  int irang,double astarti,double aendi,
		  double astartj,double aendj,double ealfai[],
		  double ebetai[],int nmbcoef1,double ealfaj[],
		  double ebetaj[],int nmbcoef2);


    /// Test function: Check continuity between two adjacent surfaces.
    void checkContinuity(const SplineCurve* pos,
			 const SplineCurve* der1, 
			 const SplineCurve* der2,
			 const SplineCurve* cross,
			 double& minAng, double& medAng, 
			 double& maxAng);

};// end anonymous namespace


//===========================================================================
void
Go::CoonsPatchGen::blendcoef(double evecu[],double evecv[],double etang[],
			       int idim,int isign,double *coef1, double *coef2)
//===========================================================================
/*
*********************************************************************
*                                                                   
* PURPOSE    : Given three vectors, evecu, evecv and etang, find the 
*              coefficients, coef1 and coef2, such that the vector
*              coef1*evecu + coef2*evecv, is as close as possible
*              to the vector etang. If the three vectors lie in a
*              plane, coef1*evecu + coef2*evecv = etang.
*              
*
* INPUT      : evecu      - First vector.
*              evecv      - Second vector.
*              etang      - Vector to approximate.
*              idim       - Dimension of geometry space.
*              isign      - Sign with wich etang is to be multiplied.
*
*
* OUTPUT     : coef1      - First factor.
*              coef2      - Second factor.
*              jstat      - status messages  
*                                         > 0      : warning
*                                         = 0      : ok
*                                         < 0      : error
*
*
* METHOD     : Minimize the square of the expression 
*                   dist(coef1*evecu+coef2*evecv,isign*etang)
*              over coef1 and coef2.
*              The expression is differentiated and set equal to
*              zero. Then this equation system of 2 equations
*              with two unknowns is solved.
*
*********************************************************************
*/
{
  double epsco = 1e-14;
  double REL_PAR_RES = 1e-6;
  
  // int kstat = 0;           /* Status variable.                    */
  int ki;                  /* Counter.                            */
  double tdotuu;           /* Scalar product of evecu and evecu.  */
  double tdotuv;           /* Scalar product of evecu and evecv.  */
  double tdotutang;        /* Scalar product of evecu and etang.  */
  double tdotvv;           /* Scalar product of evecv and evecv.  */
  double tdotvtang;        /* Scalar product of evecv and etang.  */
  double tdiv;             /* Determinant of equation system.     */
  
  /* Set output to zero. */
  //  *jstat = 0;
  *coef1 = (double)0.0;
  *coef2 = (double)0.0;
  
  /* Compute coefficients of equation system.  */


  //  tdotuu = s6scpr(evecu,evecu,idim);
  tdotuu = inner(evecu, evecu + idim, evecu);
  //tdotuu = std::inner_product(evecu, evecu + idim, evecu, 0.0);
  //  tdotuv = s6scpr(evecu,evecv,idim);
  tdotuv = inner(evecu, evecu + idim, evecv);
  //tdotuv = std::inner_product(evecu, evecu + idim, evecv, 0.0);
  //  tdotutang = (double)isign*s6scpr(evecu,etang,idim);
  tdotutang = (double)isign*inner(evecu, evecu + idim, etang);
  //tdotutang = (double)isign*(std::inner_product(evecu, evecu + idim, 
  //					etang, 0.0));
  //  tdotvv = s6scpr(evecv,evecv,idim);
  tdotvv = inner(evecv, evecv + idim, evecv);
  //tdotvv = std::inner_product(evecv, evecv + idim, evecv, 0.0);
  //  tdotvtang = (double)isign*s6scpr(evecv,etang,idim);
  tdotvtang = (double)isign*inner(evecv, evecv + idim, etang);
  //tdotvtang = (double)isign*(std::inner_product(evecv, evecv + idim, 
  //					etang, 0.0));

  tdiv = tdotuv*tdotuv - tdotuu*tdotvv;
  if (fabs(tdiv) < epsco)
    {
	if (fabs(tdotuu)<epsco && fabs(tdotvv)<epsco);
	else
	    if (fabs(tdotuu)<epsco)
	      {
		//	  *coef2 = s6length(etang,idim,&kstat)/sqrt(tdotvv);
		*coef2 = sqrt(sum_squared(etang,etang + idim)/tdotvv);
		if (inner(evecv, evecv + idim, etang) < 0.0)
		  *coef2 = -1.0*(*coef2);
	      }
	    else
	      {
		//	*coef1 = s6length(etang,idim,&kstat)/sqrt(tdotuu);
		*coef1 = sqrt(sum_squared(etang,etang + idim)/tdotuu);
		if (inner(evecu, evecu + idim, etang) < 0.0)
		  *coef1 = -1.0*(*coef1);
	      }
	return;
    }
  
  /* Compute output factors.  */

  *coef1 = (tdotvtang*tdotuv - tdotutang*tdotvv)/tdiv;
  *coef2 = (tdotutang*tdotuv - tdotvtang*tdotuu)/tdiv;

  /* Test result.  */

  for (ki=0; ki<idim; ki++)
    if (fabs((*coef1)*evecu[ki]+(*coef2)*evecv[ki]-isign*etang[ki]) > 
	REL_PAR_RES) 
      break;
  
  if (ki < idim)
    { 
//       std::cout << "Missing equality in cross tangent ";
//       std::cout << (*coef1)*evecu[ki]+(*coef2)*evecv[ki]-isign*etang[ki];
//       std::cout << std::endl;
    }
//        MESSAGE("Equality not achieved.");
  //    *jstat = 2;   // Equality not achieved

  return;
}


//===========================================================================
void
Go::CoonsPatchGen::hermit(double econd[], int icond, bool hasder1,
			  double astart, double aend, int idim)
//===========================================================================
/*
*********************************************************************
*                                                                   
* PURPOSE    : Hermite interpolation of position derivative in the
*              two endpoints represented as a Bezier curve on the interval 
*              [astart,aend]. icond is expected to be less than or equal to 4.
*              
*
*
*
* INPUT      : icond      - Number of interpolation conditions. 
*              hasder1    - Whether derivative is given in the start of the
*                           interval.
*              astart     - Start of parameter interval.
*              aend       - End of parameter interval
*              idim       - Dimension of geometry space.
*
*
* INPUT/OUTPUT : econd    - Interpolation conditions as input, Bezier coefficients
*                           as output. The dimension is icond*idim.
*                       
*
* OUTPUT     : jstat      - status messages  
*                                         > 0      : warning
*                                         = 0      : ok
*                                         < 0      : error
*
*
*********************************************************************
*/
{
  // Test number of conditions

    //  *jstat = 0;
  if (icond < 2 || icond > 4) 
      THROW("Wrong number of conditions.");
      //    {
      // *jstat = -1;
      // return;
      //  }
  
  int ki;
  if (icond == 4)
    {
      // Hermit interpolation with Bezier curve of order 4.  
      for (ki=0; ki<idim; ki++)
	{
	  econd[idim+ki] = (aend-astart)*econd[idim+ki]/3.0 + econd[ki];
	  econd[2*idim+ki] = -(aend-astart)*econd[2*idim+ki]/3.0 + 
	    econd[3*idim+ki];
	}
    }
  else if (icond == 3)
    {
      // Interpolation with Bezier curve of order 3.  
      for (ki=0; ki<idim; ki++)
	econd[idim+ki] = hasder1 ? 0.5*(aend-astart)*econd[idim+ki]+econd[ki] :
                              econd[2*idim+ki]-0.5*(aend-astart)*econd[idim+ki];
    }
  else if (icond == 2) {}  // Nothing to do.
  
  return;
}


// Size of curves = iedge * 3 (i.e. bd_curve, cross_curve, tangent curve).
// and in what ordering they are stored. Same applies for other functions.
// Boundary curves form a loop, and all curves share orientation and parametrization.
// cross-curves point inwards.
//===========================================================================
void
Go::CoonsPatchGen::getTangBlends(vector<shared_ptr<SplineCurve> >& curves,
				   // Boundary curves and tangent curves 
				   int iedge,
				   // Blending functions
				   vector<shared_ptr<SplineCurve> >& blend)
//===========================================================================
//
// PURPOSE: Find blending functions used to blend two derivative
//          along some boundary curves corresponding to a surface,
//          into a cross derivative curve pr edge.
// //
// //
// //
// CALLS:
// //
// //
// WRITTEN BY: Vibeke Skytt, SINTEF, 9801.
// //
// //
// REFERENSE:  sh1262 in SISL.
// //
// //
// REVISED BY:
//
//
////////////////////////////////////////////////////////////////////////
{
  double epsang = 1e-8; // @@@ Check value. Check all '1e-*'.
  double epsco = 1e-14;
  //  double REL_PAR_RES = 1e-8;


  // int kstat = 0;     /* Status variable.     */
  // int kstat1 = 0;    /* Status variable used to check if the requirements
  //			to input is satisfied.         */
  double pihalf = 3.141592653589793/2.0;
  int ki,kj,kk;      /* Counters.            */
  int ki4,ki12;      /* ki4 = 4*ki, and ki12 = 12*ki.  */
  int kj4,kj12;      /* kj4 = 4*kj, and kj12 = 12*kj.  */
  int kdim = 3;      /* Dimension of geometry space.   */
  int kder = 1;      /* Number of derivatives of curve to evaluate. */
  int kncurve = 3*iedge;  /* Number of input curves.   */
  int kncond = 4*iedge;   /* Number of coeffecients of blending
                             functions corresponding to derivative
                             curves in one parameter direction.    */
  int ksign;              /* Indicates if the tangent along a
                             position curve and the cross tangent
                             of the adjacent edge in a corner is
                             to have the same sign.                */
  int krang;              /* Rang of equation system used to find derivatives
                             in endpoints  of blending functions. */
  double ang_tol = epsang; // Tolerance used in computing the
				       // space spanned by tang. vectors. 
  double tfac;            /* Used in linearity test of blending function. */
  double tref=(double)0.0; /* Reference value in equality test.  */
  double tpar1,tpar2;     /* Endpoints of parameter interval of curve. */
  double thelp1,thelp2,thelp3,thelp4,thelp5;   /* Help variables in 
					          computation of rang.    */
  double tnorm1,tnorm2;   /* Lenght of normals.        */
  double *spar1 = NULL;   /* Startpoints of parameter intervals of 
                             input curves.       */
  double *spar2 = NULL;   /* Endpoints of parameter intervals of 
                             input curves.       */
  double *sder = NULL;    /* Result of curve evaluation. The values are
                             stored as follows : Value and first derivative
                             of all curves corresponding to an edge in
                             first endpoint, the same values in the second
                             endpoint. All curves of one edge is treated
                             before the curve of the next edge is treated.
                             In total 48 points = 144 doubles. */
  double sa[12];          /* Matrix of equation system used to find
                             derivatives of blending functions.     */
  double sb[3];           /* Right side of equation system.         */
  double snorm1[3];       /* Normal of tangent plane at a corner of
                             the vertex region.                     */
  double snorm2[3];       /* Normal of tangent plane at a corner of
                             the vertex region.                     */
  //  SISLCurve *qc;              // Pointer to curve.                      
  shared_ptr<SplineCurve> qc;

  // Allocate space for array containing results of curve evaluation.  
  vector<double> scratch(12*kdim*iedge + 2*kncurve + 8*iedge);
  sder = &scratch[0];
  spar1 = sder + 12*kdim*iedge;
  spar2 = spar1 + kncurve;
  double* ecoef = spar2 + kncurve;   // Coefficients of blending curve

  // Number of coefficients for each blending curve
  vector<int> nmbcoef(2*iedge);
  
  // Compute number of coeffients in the blending curves for each boundary.
  for (ki=0; ki<iedge; ki++)
    {
      kj = (ki > 0) ? ki-1 : iedge-1;
      kk = (ki < iedge-1) ? ki+1 : 0;

      //      if (!(curves[3*ki+1].hasCurve()))
      if (curves[3*ki+1].get() == NULL)
      //    if (curves[3*ki+2].get() == NULL)
	nmbcoef[ki] = 0;   // No blending curve for this boundary
      else if ((curves[3*kj+1].get() != NULL) &&
	       (curves[3*kk+1].get() != NULL))
	nmbcoef[ki] = 4;   // Tangent conditions on both adjacent boundaries
      else if ((curves[3*kj+1].get() != NULL) ||
	       (curves[3*kk+1].get() != NULL))
	nmbcoef[ki] = 3;   // Tangent conditions on one adjacent boundary
      else
	nmbcoef[ki] = 2;   // No tangent conditions on adjacent boundaries

//        if (nmbcoef[ki]>0)
//  	nmbcoef[ki]=4;
    }
      
  // Evaluate all boundary curves in the endpoints.  
  for (ki=0; ki<kncurve; ki++)
    {
      kj = ki/3;
      kk = ki % 3;

      qc = curves[ki];
      if (qc.get() == NULL) 
	continue;   // Missing tangent curve.
     
      // Fetch parameter values at endpoints.  
      spar1[ki] = tpar1 = qc->startparam();//*(qc->et + qc->ik - 1);
      spar2[ki] = tpar2 = qc->endparam();//*(qc->et + qc->in);

      // Evaluate the current curve in the startpoint.  
      std::vector<Point> pts(kder + 1);
      qc->point(pts, tpar1, kder);
      // It's now time to transfer the info from 'pts' to 'sder+...'.
      for (size_t i = 0; i < pts.size(); ++i)
	  for (int j = 0; j < qc->dimension(); ++j)
	      sder[(2*(6*kj+kk)+i)*kdim + j] = pts[i][j];

      // Evaluate the curve in the endpoint.  
      qc->point(pts, tpar2, kder);
      // It's now time to transfer the info from 'pts' to 'sder+...'.
      for (size_t i = 0; i < pts.size(); ++i)
	  for (int j = 0; j < qc->dimension(); ++j)
	      sder[(2*(6*kj+kk+3)+i)*kdim + j] = pts[i][j];
    }
  
  // Compute first and last coefficients of blending functions, i.e find
  // value of blending functions in the endpoints.  
  for (ki=0; ki<iedge; ki++)
    {
      if (nmbcoef[ki] == 0) continue;  // No blending functions to compute

      // For each edge, make sure that the endpoints of the cross tangent 
      // curve is equal (possibly expect from a sign) to the tangents of the
      // position curve along the adjacent edges.  

      kj = (ki > 0) ? ki-1 : iedge-1;
      ki4 = 4*ki, ki12 = 12*ki;
      kj4 = 4*kj, kj12 = 12*kj;
      
      // Treat startpoint of edge.  
      ksign = -1;
      blendcoef(sder+(ki12+2)*kdim,sder+(ki12+4)*kdim,sder+(kj12+7)*kdim,
		kdim,ksign,ecoef+ki4,ecoef+kncond+ki4);
      
      // Treat endpoint of edge. 
      ksign = 1;
      blendcoef(sder+(ki12+8)*kdim,sder+(ki12+10)*kdim,
		sder+(((ki+1)%iedge)*12+1)*kdim,
		kdim,ksign,ecoef+ki4+nmbcoef[ki]-1,
		ecoef+kncond+ki4+nmbcoef[ki]-1);
    }
  
  for (tref=(double)0.0, ki=0; ki<iedge; ki++)
    {
      // For each edge, set up an equation system to find the derivatives
      // of the blending functions excisting at the corner lying at the 
      // startpoint of the edge (also the blending functions corresponding
      // to the previous edge). The equation system represent the conditions
      // to be satisfied to have a consistent twist vector in this corner. 

      kj = (ki > 0) ? ki-1 : iedge-1;

      ki4 = 4*ki, ki12 = 12*ki;
      kj4 = 4*kj, kj12 = 12*kj;
      
      // Check if two edges with tangent curves meet in this corner
      if ((curves[3*ki+1].get() == NULL) || (curves[3*kj+1].get() == NULL))
	{
//  	  ecoef[ki4+1] = ecoef[kncond+ki4+1] = ecoef[kj4+2] = 
//  	    ecoef[kncond+kj4+2] = 0.0;
	  continue;    // No twist conditions to satisfy
	  }

      for (kk=0; kk < kdim; kk++)
	{
          sa[kk*4] = sder[(ki12+2)*kdim+kk];
	  sa[kk*4+1] = sder[(ki12+4)*kdim+kk];
	  sa[kk*4+2] = sder[(kj12+8)*kdim+kk];
	  sa[kk*4+3] = sder[(kj12+10)*kdim+kk];

	  sb[kk] = - ecoef[ki4]*sder[(ki12+3)*kdim+kk]
	    - ecoef[kncond+ki4]*sder[(ki12+5)*kdim+kk]
	    - ecoef[kj4+nmbcoef[kj]-1]*sder[(kj12+9)*kdim+kk]
		- ecoef[kncond+kj4+nmbcoef[kj]-1]*sder[(kj12+11)*kdim+kk];
	  
	  tref = std::max(tref,std::max(sa[kk*4],sa[kk*4+1]));
	  tref = std::max(tref,std::max(sa[kk*4+2],sa[kk*4+3]));
	}

      // Find rang of equation system, i.e. check wether the 4
      // tangent vectors are able to span a plane.   

      //      s6crss(sder+(ki12+2)*kdim,sder+(ki12+4)*kdim,snorm1);
      Vector3D x, y, z;
      x.setValue(sder+(ki12+2)*kdim);
      y.setValue(sder+(ki12+4)*kdim);
      z = x.cross(y);
      int i;
      for (i = 0; i < 3; ++i) snorm1[i] = z[i];
      //      snorm1 = z.begin();
      tnorm1 = sqrt(sum_squared(snorm1,snorm1+kdim));

      //      s6crss(sder+(kj12+8)*kdim,sder+(kj12+10)*kdim,snorm2);
      x.setValue(sder+(ki12+8)*kdim);
      y.setValue(sder+(ki12+10)*kdim);
      z = x.cross(y);
      for (i = 0; i < 3; ++i) snorm2[i] = z[i];
      //      snorm2 = z.begin();
      tnorm2 = sqrt(sum_squared(snorm2,snorm2+kdim));

      //      thelp1 = s6ang(snorm1,sder+(kj12+8)*kdim,kdim);
      x.setValue(snorm1);
      y.setValue(sder+(ki12+8)*kdim);
      thelp1 = x.angle(y);
      //      thelp2 = s6ang(snorm1,sder+(kj12+10)*kdim,kdim);
      x.setValue(snorm1);
      y.setValue(sder+(ki12+10)*kdim);
      thelp2 = x.angle(y);

      //      thelp3 = s6ang(snorm2,sder+(ki12+2)*kdim,kdim);
      x.setValue(snorm2);
      y.setValue(sder+(ki12+2)*kdim);
      thelp3 = x.angle(y);
      //      thelp4 = s6ang(snorm2,sder+(ki12+4)*kdim,kdim);
      x.setValue(snorm2);
      y.setValue(sder+(ki12+4)*kdim);
      thelp4 = x.angle(y);

      //      thelp5 = s6ang(snorm1,sb,kdim);
      x.setValue(snorm1);
      y.setValue(sb);
      thelp5 = x.angle(y);

      z.setValue(snorm2);
      if (fabs(x.angle(y)-pihalf) < ang_tol &&
	  fabs(thelp5) > ang_tol)
	krang = 1;
      else if ((fabs(tnorm1) < epsco && 
	   fabs(thelp1-pihalf) < ang_tol &&
	   fabs(thelp2-pihalf) < ang_tol) ||
	  (fabs(tnorm2) < epsco &&
	   fabs(thelp3-pihalf) < ang_tol &&
	   fabs(thelp4-pihalf) < ang_tol) ||
	  (fabs(tnorm1) < epsco &&
	   (fabs(thelp1-pihalf) < ang_tol ||
 	    fabs(thelp2-pihalf) < ang_tol)) ||
	  (fabs(tnorm2) < epsco &&
	   (fabs(thelp3-pihalf) < ang_tol ||
 	    fabs(thelp4-pihalf) < ang_tol)) ||
	  (fabs(tnorm1) < epsco && 
	   fabs(tnorm2) < epsco))
	krang = 2;
      else
	krang = 3;
      
      //      thelp1 = s6ang(snorm1,sder+(ki12+3)*kdim,kdim);
      x.setValue(snorm1);
      y.setValue(sder+(ki12+3)*kdim);
      thelp1 = x.angle(y);
      //      thelp2 = s6ang(snorm1,sder+(ki12+5)*kdim,kdim);
      x.setValue(snorm1);
      y.setValue(sder+(ki12+5)*kdim);
      thelp2 = x.angle(y);
      //      thelp3 = s6ang(snorm1,sder+(kj12+9)*kdim,kdim);
      x.setValue(snorm1);
      y.setValue(sder+(ki12+9)*kdim);
      thelp3 = x.angle(y);
      //      thelp4 = s6ang(snorm1,sder+(kj12+11)*kdim,kdim);
      x.setValue(snorm1);
      y.setValue(sder+(ki12+11)*kdim);
      thelp4 = x.angle(y);

//       for (ki=0; ki<kdim; ki++)
// 	snorm1[ki] += snorm2[ki];
      int kh;
      for (kh=0; kh<kdim; kh++)
	snorm1[kh] += snorm2[kh];

      if (sum_squared(snorm1,snorm1+kdim) != 0)
	  normalize(snorm1,snorm1+kdim);
      //      s6norm(snorm1, kdim, snorm1, &kstat);
      double t1 = std::inner_product(sb, sb+kdim, snorm1, 0.0);
      //      double t1 = s6scpr(sb, snorm1, kdim);

      double sc[3];
//       for (ki=0; ki<kdim; ki++)
// 	sc[ki] = sb[ki] - t1*snorm1[ki];
      for (kh=0; kh<kdim; kh++)
	sc[kh] = sb[kh] - t1*snorm1[kh];

      krang = 2;
      
      // Estimate derivatives in the current endpoints of the blending 
      // functions existing in this corner.  

	blendder(sa,sc,4,kdim,krang,spar1[3*ki],spar2[3*ki],spar1[3*kj],
		   spar2[3*kj],ecoef+ki4,ecoef+kncond+ki4,nmbcoef[ki],
		   ecoef+kj4,ecoef+kncond+kj4,nmbcoef[kj]);
	//      if (kstat < 0) 
	// return; // GR_SISLERR;
	//  kstat1 = max(kstat,kstat1);
    }
   
  for (ki=0; ki<iedge; ki++)
    {
      int ncf = nmbcoef[ki];
      if (ncf == 0) continue;

      kj = (ki > 0) ? ki-1 : iedge-1;

      // Interpolate blending curves.  

      hermit(ecoef+4*ki, nmbcoef[ki], (curves[3*kj+1].get() != NULL),
	     spar1[3*ki], spar2[3*ki], 1);
      //   if (kstat < 0) 
      //  return; // GR_SISLERR;

      hermit(ecoef+kncond+4*ki, nmbcoef[ki], (curves[3*kj+1].get() != NULL),
	     spar1[3*ki], spar2[3*ki], 1);
      //  if (kstat < 0)
      //  return; // GR_SISLERR;
      
      // Test blending curves. If the coeffecients deviate much from
      // those of a linear curve, they are linearized. In addition, the
      // coefficients of the blending curve corresponding to the input
      // cross derivative curve is not allowed to be negative.      
      
      tfac = std::max(std::max(fabs(ecoef[4*ki+nmbcoef[ki]-1]),fabs(ecoef[4*ki])),
		 std::max(fabs(ecoef[kncond+4*ki+nmbcoef[ki]-1]),
		     fabs(ecoef[kncond+4*ki])));
      for (kj=1; kj<ncf-1; kj++)
      {
	 if (ecoef[4*ki+kj] < std::min(ecoef[4*ki+ncf-1],ecoef[4*ki])-tfac ||
	     ecoef[4*ki+kj] > std::max(ecoef[4*ki+ncf-1],ecoef[4*ki])+tfac)
	 {
	    ecoef[4*ki+kj] = ((double)(ncf-1-kj)*ecoef[4*ki] +
	       (double)kj*ecoef[4*ki+ncf-1])/(double)(ncf-1);
	    //  kstat1 = max(kstat1,1);
	 }
	 
	 if (ecoef[kncond+4*ki+kj] < 
	     std::min(ecoef[kncond+4*ki+ncf-1],ecoef[kncond+4*ki])-tfac ||
	     ecoef[kncond+4*ki+kj] > 
	     std::max(ecoef[kncond+4*ki+ncf-1],ecoef[kncond+4*ki])+tfac)
	 {
	    ecoef[kncond+4*ki+kj] = ((double)(ncf-1-kj)*ecoef[kncond+4*ki] +
	       (double)kj*ecoef[kncond+4*ki+3])/(double)(ncf-1);
	    //   kstat1 = max(kstat1,1);
	 }
      }
    }
  
  
  // Coefficients of blending functions found. Make knot vectors.
  double et[8];   // Maximum number of knots.
  for (ki=0; ki<iedge; ki++)
    {
      if (nmbcoef[ki] == 0) 
	{
	  blend.push_back(shared_ptr<SplineCurve>());
	  blend.push_back(shared_ptr<SplineCurve>());
	  continue;
	}

      for (kj=0; kj<nmbcoef[ki]; kj++)
	{
	  et[kj] = spar1[3*ki];
	  et[nmbcoef[ki]+kj] = spar2[3*ki];
	}

      // Create sisl curve.
      //      SISLCurve *qc = newCurve(nmbcoef[ki], nmbcoef[ki], et, 
      //		       ecoef+4*ki, 1, 1, 1);
      shared_ptr<SplineCurve> qc = shared_ptr<SplineCurve>
	  (new SplineCurve(nmbcoef[ki], nmbcoef[ki], et,
			     ecoef+4*ki, 1));

      // Make GRSislCurve
      blend.push_back(qc);

      //      qc = newCurve(nmbcoef[ki], nmbcoef[ki], et, 
      //		    ecoef+kncond+4*ki, 1, 1, 1);
      qc = shared_ptr<SplineCurve>
	  (new SplineCurve(nmbcoef[ki], nmbcoef[ki], et,
			     ecoef+kncond+4*ki, 1));

      // Make GRSislCurve
      blend.push_back(qc);
    }

  return; // (kstat1 == 0) ? GR_OK : GR_CONTINUE;
}


//===========================================================================
void
Go::CoonsPatchGen::fixCrossEndPts(const vector<shared_ptr<SplineCurve> >& bd_curves,
				  const vector<shared_ptr<SplineCurve> >& cross_curves)
//===========================================================================
/*
*********************************************************************
*                                                                   
* PURPOSE    : The produced cross tangent curves are in some corners not
*              consistent with the tangents of the adjacent curves.
*              Change the cross tangent curves in order to achieve this
*              consistency.
*
*
*
* INPUT      : boundary   - Boundary curves.
*              crosstan   - Cross boundary curves.
*                       
*
* OUTPUT     : crosstan   - Cross boundary curves after modification.
*              jstat      - status messages  
*                                         > 0      : warning
*                                         = 0      : ok
*                                         < 0      : error
*
*
* METHOD     : 
*
*
* REFERENCES : 
*              
*
* USE        :
*
*-
* CALLS      : s1221, s6lenght
*              
*
* WRITTEN BY : Vibeke Skytt, SI, 11.92.
*
*********************************************************************
*/
{
    double epsge = 1e-8; //

    int ki,kj,kk,kl;              /* Counters.                    */
    int kn;                       /* Number of vertices of cross tangent curve.*/
    int kn2;                      /* Half the number of vertices. */
    int k1;                       /* Counter.                     */
    int kdim = bd_curves[0]->dimension(); /* Dimension of geometry space. */
    double tpar1,tpar2;           /* Parameter values of endpoints of curves. */
    double tmult;                 /* Multiplicator used to smoothen curve.    */
    double sdiff[3];              /* Difference between wanted and actual
				     cross tangent vector.                    */
    double sder[9];               /* Result of curve evalutation.             */
    shared_ptr<SplineCurve> qcross;  /* Cross tangent curve.                */
    shared_ptr<const SplineCurve> qc; /* Adjacent position curve.                 */
   
#ifdef CREATORS_DEBUG
    std::cout << "Dim CoonsUtils " << kdim << std::endl;
#endif

    // Traverse the edges.
    for (ki=0; ki < (int)bd_curves.size(); ki++) {
	if (cross_curves[ki].get() == NULL) continue;  // No cross tangent curve

	qcross = cross_curves[ki];
      
	// Check if the cross tangent curve is consistent with
	// the adjacent tangent in the startpoint of the curve. 
	// First fetch adjacent curve.   
      
	kj = (ki > 0) ? ki-1 : (int)bd_curves.size()-1;
	qc = bd_curves[kj];
      
	/* Fetch parameter values at the current corner.  */

	//      tpar1 = *(qc->et + qc->in);
	tpar1 = qc->endparam();
	//      tpar2 = *(qcross->et + qcross->ik - 1);
	tpar2 = qcross->startparam();

	// Evaluate curves in corner.  
	std::vector<Point> pts(3);
	qc->point(pts, tpar1, 1);
	// It's now time to transfer the info from 'pts' to 'sder+...'.
	for (size_t i = 0; i < pts.size(); ++i)
	    for (int j = 0; j < kdim; ++j)
		sder[i*kdim + j] = pts[i][j];
      
	qcross->point(pts[2], tpar2);
	// It's now time to transfer the info from 'pts' to 'sder+...'.
	for (int i = 2; i < 3; ++i)
	    for (int j = 0; j < kdim; ++j)
		sder[i*kdim + j] = pts[i][j];
      
	// Check equality in each dimension, opposite sign of vectors.
      	for (kk=0; kk<kdim; kk++) // The tangent of qc points outwards, cross inwards.
	    if (fabs(sder[kdim+kk]+sder[2*kdim+kk]) >= epsge)
		break;
      
	if (kk < kdim) {
	    // Inconsistence. Change the first coefficient of the
	    // cross tangent curve to be equal to the tangent (the cross 
	    // tangent curve is k-regular). Then smooth this change
	    // out along the cross tangent curve.  
	 
	    // The sign of the cross tangent curve and the sign of the
	    // tangent is opposite.   

	    for (kl=0; kl<kdim; kl++)
		sder[kdim+kl] *= -(double)1.0;

	    // Find difference between coefficient of cross tangent
	    // curve and wanted coefficient.  

	    //	 s6diff(sder+kdim,qcross->coefs_begin(),kdim,sdiff);
	    for (int i = 0; i < kdim; ++i)
		sdiff[i] = sder[kdim + i] - qcross->coefs_begin()[i];
	 
	    // Change cross tangent curve, most at the corner, then
	    // smoothen out the effect of the change. 
	    
	    for (kn=qcross->numCoefs(), kn2=std::min(kn/2,2), k1=0, tmult=(double)1.0;
		 k1<kn2; k1++, tmult /= 2.0)
		for (kl=0; kl<kdim; kl++)
		    qcross->coefs_begin()[k1*kdim+kl] += tmult*sdiff[kl];
	}
	    
	// Check if the cross tangent curve is consistent with
	// the adjacent tangent in the endpoint of the curve. 
	// First fetch adjacent curve.   
      
	kj = (ki < ((int)bd_curves.size()-1)) ? ki+1 : 0;
	qc = bd_curves[kj];
      
	// Fetch parameter values at the current corner.  
      
	tpar1 = *(qc->basis().begin() + qc->order() - 1);
	tpar2 = *(qcross->basis().begin() + qcross->numCoefs());
      
	// Evaluate curves in corner.  
	qc->point(pts, tpar1, 1);
	// It's now time to transfer the info from 'pts' to 'sder+...'.
	for (size_t i = 0; i < pts.size(); ++i)
	    for (int j = 0; j < kdim; ++j)
		sder[i*kdim + j] = pts[i][j];

	qcross->point(pts[2], tpar2);
	// It's now time to transfer the info from 'pts' to 'sder+...'.
	for (int i = 2; i < 3; ++i)
	    for (int j = 0; j < kdim; ++j)
		sder[i*kdim + j] = pts[i][j];

	// Check equality in each dimension.
	for (kk=0; kk<kdim; kk++)
	    if (fabs(sder[kdim+kk]-sder[2*kdim+kk]) >= epsge)
		break;
      
	if (kk < kdim) {
	    // Inconsistence. Change the last coefficient of the
	    // cross tangent curve to be equal to the tangent (the cross 
	    // tangent curve is k-regular). Then smooth this change
	    // out along the cross tangent curve.  
	 
	    // Find difference between coefficient of cross tangent
	    // curve and wanted coefficient.  
	 
	    ///	 s6diff(sder+kdim,
	    //qcross->coefs_begin()+(qcross->numCoefs()-1)*kdim,
	    //kdim,sdiff);
	    for (int i = 0; i < kdim; ++i)
		sdiff[i] = sder[kdim + i] -
		    qcross->coefs_begin()[(qcross->numCoefs() - 1) * kdim + i];

	 
	    // Change cross tangent curve, most at the corner, then
	    // smoothen out the effect of the change. 
	    
	    for (kn=qcross->numCoefs(), kn2=std::max(kn/2,kn-3),
		     k1=kn-1, tmult=(double)1.0; k1>kn2; k1--, tmult /= 2.0)
		for (kl=0; kl<kdim; kl++)
		    qcross->coefs_begin()[k1*kdim+kl] += tmult*sdiff[kl];
	}
      
	// Make sure that the 4 coefficients (boundary curve and cross
	// tangent curve) in each corner lies in a plane.
	// First corner
	if (false) {
	    qc = bd_curves[ki];
	    qcross = cross_curves[ki];
	    double vec[3];
	    double norm[3];
	    double t1;
	    //	   s6diff(qc->coefs_begin()+kdim, qc->coefs_begin(), kdim, vec);
	    for (int i = 0; i < kdim; ++i)
		sdiff[i] = qc->coefs_begin()[kdim + i] - qc->coefs_begin()[i];
	    //     s6crss(vec, qcross->coefs_begin(), norm);
	    Vector3D x, y, z;
	    x.setValue(vec);
	    y.setValue(&(*qcross->coefs_begin()));
	    z = x.cross(y);
	    for (int i = 0; i < 3; ++i) norm[i] = z[i];

	    normalize(norm, norm+kdim);
	    //     s6norm(norm, kdim, norm, &kstat);
	    t1 = std::inner_product(qcross->coefs_begin() + kdim,
				    qcross->coefs_begin() + 2*kdim, norm, 0.0);
	    //     t1 = s6scpr(qcross->coefs_begin()+kdim, norm, kdim);
	    for (kl=0; kl<kdim; kl++)
		qcross->coefs_begin()[kdim+kl] -= t1*norm[kl];

	    // Second corner
	    kn = qc->numCoefs();
	    kn2 = qcross->numCoefs();
	    //	   s6diff(qc->coefs_begin()+(kn-2)*kdim,
	    //  qc->coefs_begin()+(kn-1)*kdim, kdim, vec);
	    for (int i = 0; i < kdim; ++i)
		sdiff[i] = qc->coefs_begin()[(kn - 2) * kdim + i] -
		    qc->coefs_begin()[(kn - 1) * kdim + i];
	    //     s6crss(vec, qcross->coefs_begin()+(kn2-1)*kdim, norm);
	    x.setValue(vec);
	    y.setValue(&(*(qcross->coefs_begin()+(kn2-1)*kdim)));
	    z = x.cross(y);
	    for (int i = 0; i < 3; ++i)
		norm[i] = z[i];

	    normalize(norm, norm+kdim);
	    //     s6norm(norm, kdim, norm, &kstat);
	    //	   t1 = s6scpr(qcross->coefs_begin()+(kn2-2)*kdim, norm, kdim);
	    t1 = std::inner_product(qcross->coefs_begin()+(kn-2)*kdim,
				    qcross->coefs_begin()+(kn-2)*kdim+kdim, norm, 0.0);
	    for (kl=0; kl<kdim; kl++)
		qcross->coefs_begin()[(kn2-2)*kdim+kl] -= t1*norm[kl];
	}
    }

    return;
}     
  

//===========================================================================
void
Go::CoonsPatchGen::getCrossTangs(const vector<shared_ptr<SplineCurve> >& curves,
				   vector<shared_ptr<SplineCurve> >&
				   mod_cross_curves,
				   double tol1, double tol2)
//===========================================================================
//getCrossTangs(GRSislCurveArray& curve, // Boundary curves and
//                                                 // tangent curves 
//		  int               iedge,  // # of edges,
//		  GRSislCurveArray& crosstan)  // Cross tangents
////////////////////////////////////////////////////////////////////////
// //
// PURPOSE: Prepare cross tangents for surface creation.
// //
// //
// //
// CALLS:
// //
// //
// WRITTEN BY: Vibeke Skytt, SINTEF, 9801.
// //
// //
// REVISED BY:
// //
// //
// REFERENCE:  sh1263 in SISL.
// //
// METHOD     : If necessary reparametrize the edge-curve in order to avoid 
//              tangent lengths of the position curves which very long.
//              Blend the derivative curves to produce the cross tangent
//              curves in such a way that :
//              1. In a corner of the region, the cross tangent belonging to
//                 one of the adjacent edges is equal to the tangent belonging
//                 to the other adjacent edge (except for a sign).
//              2. The derivatives of the two cross tangent curves meeting in a
//                 corner are equal, i.e. consistent twist.
// //
////////////////////////////////////////////////////////////////////////
{
    int iedge = (int)curves.size()/2;

    //  GRStatus grstat = GR_OK;   // Status variable.  
    // int kstat = 0;        // Status variable.  
    // Just to be safe...
    mod_cross_curves.clear();
    // Copy input curves
    vector<shared_ptr<SplineCurve> > cpCurve;
    cpCurve.reserve(curves.size() * 3 / 2);
  
    for (size_t ki=0; ki<curves.size()/2; ki++) {
	cpCurve.push_back(curves[2*ki]);
	cpCurve.push_back(curves[2*ki + 1]);
	SplineCurve* tangent_crv = curves[2*ki]->derivCurve(1);
	cpCurve.push_back(shared_ptr<SplineCurve>(tangent_crv));
    }

    // Modify cross tangent curves that are zero in the endpoints
    // to avoid problems later.
    // @@sbr This should be handled in a different manner!!!
    Point cross, der;
    std::vector<double>::iterator sc;
    for (int ki=0; ki<iedge; ki++) {
	if (cpCurve[3*ki+1].get() == 0)
	    continue;

	int dim = cpCurve[3*ki+1]->dimension();
	cpCurve[3*ki+1]->point(cross, cpCurve[3*ki+1]->startparam());
	if (cross.length() < tol1) {
	    int kj = (ki == 0) ? iedge-1 : ki-1;
	    cpCurve[3*kj+2]->point(der, cpCurve[3*kj+2]->endparam());
	    der.normalize();
	    cross = -tol1*der;
	    sc = cpCurve[3*ki+1]->coefs_begin();
	    for (int kh=0; kh<dim; kh++)
		sc[kh] = cross[kh];
	}
	cpCurve[3*ki+1]->point(cross, cpCurve[3*ki+1]->endparam());
	if (cross.length() < tol1) {
	    int kj = (ki + 1)%iedge;
	    cpCurve[3*kj+2]->point(der, cpCurve[3*kj+2]->startparam());
	    der.normalize();
	    cross = tol1*der;
	    sc = cpCurve[3*ki+1]->coefs_end();
	    for (int kh=0; kh<dim; kh++)
		sc[kh-dim-1] = cross[kh];
	}
    }

    // Find the blending functions.
    //  GRSislCurveArray blend;
    vector<shared_ptr<SplineCurve> > blend;
    getTangBlends(cpCurve, iedge, blend);
  
    // Compute the cross tangent curves.
    for (int ki=0; ki<iedge; ki++) {
	if (blend[2*ki].get() == NULL) {
	    // Empty cross tangent curve
	    mod_cross_curves.push_back
		(shared_ptr<SplineCurve>());
	    continue;
	}

	// Blend the spline curves.
	shared_ptr<SplineCurve> qc;
	// Make evaluator curve to be approximated
	int kj = (ki > 1) ? ki-2 : iedge-2+ki;
	int kh = (ki+2)%iedge;
	double fac = (kj == kh) ? 1.0 : 0.5;
	int method = (iedge == 4) ? 1 : 2;
	if (method < 10) {
	    shared_ptr<EvalCurve> offsetcrv;
	    if (method == 1) {
		offsetcrv = shared_ptr<CrossTanOffDist>
		    (new CrossTanOffDist(cpCurve[3*ki],
					   cpCurve[3*ki+1], cpCurve[3*ki+2],
					   blend[2*ki], blend[2*ki+1],
					   cpCurve[3*kj], cpCurve[3*kh], fac));
	    } else if (method == 2) {
		offsetcrv = shared_ptr<CrossTangentOffset>
		    (new CrossTangentOffset(cpCurve[3*ki],
					      cpCurve[3*ki+1], cpCurve[3*ki+2],
					      blend[2*ki], blend[2*ki+1]));
	    }

	    //       // DEBUG
	    //       std::ofstream curveout2("data/output/offsetdump2.g2");
	    //       qcoffset->writeStandardHeader(curveout2);
	    //       qcoffset->write(curveout2);
	    //       // END DEBUG
 
	    // Approximate
	    vector<double> initpars;
	    int order = cpCurve[3*ki]->order();
	    int nb_coef = cpCurve[3*ki]->numCoefs();
	    vector<double>::const_iterator knots = cpCurve[3*ki]->basis().begin();
	    initpars.push_back(knots[order-1]);
	    for (kj=order; kj<=nb_coef; kj++)
		if (knots[kj] > initpars[initpars.size()-1])
		    initpars.push_back(knots[kj]);
	    HermiteAppC approximator(offsetcrv.get(), 
				     &initpars[0], (int)initpars.size(),
				       /* tol_.neighbour */ 100.0*tol1,
				       /* tol_.kink */ /*0.9**/ tol2);
	    approximator.refineApproximation();
	    shared_ptr<SplineCurve> qcoffset = approximator.getCurve();

	    //       // DEBUG
	    //       qcoffset->writeStandardHeader(curveout);
	    //       qcoffset->write(curveout);
	    //       // END DEBUG

	    //       // debug
	    //       shared_ptr<SplineCurve> cv;
	    //       Point pt1, pt2;
	    //       double par = (9.0*cpCurve[3*ki]->startparam() + 
	    // 		    cpCurve[3*ki]->endparam())/10.0; 
	    //       cpCurve[3*ki]->point(pt1, par);
	    //       qcoffset->point(pt2, par);
	    //       cv = shared_ptr<SplineCurve>(new SplineCurve(pt1, pt2));
	    //       cv->writeStandardHeader(curveout);
	    //       cv->write(curveout);
	    //       par = (cpCurve[3*ki]->startparam() + 
	    // 	     9.0*cpCurve[3*ki]->endparam())/10.0;
	    //       cpCurve[3*ki]->point(pt1, par);
	    //       qcoffset->point(pt2, par);
	    //       cv = shared_ptr<SplineCurve>(new SplineCurve(pt1, pt2));
	    //       cv->writeStandardHeader(curveout);
	    //       cv->write(curveout);
	    //       // end debug

	    //       sh1261(cpCurve[3*ki+1], cpCurve[3*ki+2], blend[2*ki]->coefs_begin(), 
	    // 	     blend[2*ki]->numCoefs(), blend[2*ki+1]->coefs_begin(),
	    // 	     blend[2*ki+1]->numCoefs(), qc, &kstat);
	    //  if (kstat < 0) 
	    //	  return; // GR_SISLERR;
      
	    // Store the cross tangent curve.
	    //      mod_cross_curves.push_back(GRSislCurve(qc));
	    try {
		qc = shared_ptr<SplineCurve>(GeometryTools::curveSum(*qcoffset, 1,
						      *cpCurve[3*ki], -1));
	    } catch (...) {
		THROW("Failed adding curves.");
	    }
	} else {
	    //      shared_ptr<SplineCurve> qc2;
	    qc = shared_ptr<SplineCurve>
		(CurveCreators::blend(*blend[2 * ki], *cpCurve[3 * ki + 1],
					*blend[2 * ki + 1], *cpCurve[3 * ki + 2]));
	    shared_ptr<SplineCurve> cv;
	    try {
		cv = shared_ptr<SplineCurve>(GeometryTools::curveSum(*qc.get(), 1,
						      *cpCurve[3*ki], 1));
	    } catch (...) {
		THROW("Failed adding curves.");
	    }
	    // 	cv->writeStandardHeader(curveout);
	    // 	cv->write(curveout);
	}
	mod_cross_curves.push_back(qc);

	// Check
	double minang, medang, maxang;
	checkContinuity(cpCurve[3*ki].get(), cpCurve[3*ki+1].get(), 
			cpCurve[3*ki+2].get(), qc.get(), 
			minang, medang, maxang);
	//       std::cout << "Boundary no " << ki << std::endl;
	//       std::cout << "minang : " << minang << ", medang : " << medang;
	//       std::cout << ", maxang : " << maxang << std::endl;
      
#ifdef CREATORS_DEBUG
	      int k1;
	      Point p1, p2;
	      k1 = (ki+1)%iedge;
	      cpCurve[3*k1+2]->point(p1, cpCurve[3*k1+2]->startparam());
	      qc->point(p2, qc->endparam());
	      std::cout << "Derivative in end: " << std::endl;
	      p1.write(std::cout);
	      p2.write(std::cout);

	      k1 = (ki == 0) ? iedge-1 : ki-1;
	      cpCurve[3*k1+2]->point(p1, cpCurve[3*k1+2]->endparam());
	      qc->point(p2, qc->startparam());
	      std::cout << "Derivative in start: " << std::endl;
	      p1.write(std::cout);
	      p2.write(std::cout);
#endif
    }
  
    // grstat = GR_OK;
    // if (grstat == GR_CONTINUE)
    // {
    /* The tangent conditions in some corner(s) of the vertex region
       are not satisfied. Thus, the blending will not be C0 in some
       areas. Make sure that the blend is C0, sacrifising the G1-
       condition.      */

    // @@ sbr To be included when fixed. Make sure that arguments are correct (cpCurve.size()==4)
    //  fixCrossEndPts(cpCurve, mod_cross_curves, iedge);

    return; // grstat;
}


//===========================================================================
void
Go::CoonsPatchGen::addMissingCrossCurves(const vector<shared_ptr<SplineCurve> >&
					   boundary_crvs,
					   vector<shared_ptr<SplineCurve> >&
					   cross_crvs)
//===========================================================================
    //  const GRSislCurveArray& boundary_crvs,
    // GRSislCurveArray& cross_crvs)
////////////////////////////////////////////////////////////////////////
//
// PURPOSE:     Function that generates missing cross boundary curves.
//              In case of missing cross_crv, we're expecting a NULL pointer.
//
//
// CALLS:
//
//
// WRITTEN BY: Vibeke Skytt, Sintef Oslo, 9510.
//
//
// REVISED BY:
//
//
////////////////////////////////////////////////////////////////////////
{
  //  GRStatus grstat = GR_OK;
  int dim = 3;
  int nmbcrvs = (int)cross_crvs.size();
  ALWAYS_ERROR_IF(int(boundary_crvs.size()) != nmbcrvs,
		  "Number of boundary curves must equal number of cross curves.");

//   if (nmbcrvs != 4) // Missing curves are denoted by NULL pointers.
//       //   return GR_ERRINPUT;
//       GO_ERROR("addMissingCrossCurves demands 4 input curves of each type.",
// 	       InputError());

  // Create missing cross boundary curves
  vector<Point> pts(2); // We must save both point and tangent value.
  Point point;
  Point tangent;
  vector<double> points(6);
  vector<double> derivs;
  vector<int> tangent_index;
  vector<double> spar(2);
  vector<double> knots;
  //  int ltype[4];
  int nmbder = 0;
  int nmbpts = 0;
  int num_coefs; // As we're interpolating two points, while paying attention
  int order;     // to end tangents, we use a cubic B-spline.
  double par;
  for (int ki=0; ki<nmbcrvs; ki++)
    if (cross_crvs[ki].get() == NULL)
	//    if (!cross_crvs[ki].hasCurve())
      {
        order = 4;
	num_coefs = 4;
	derivs.clear();
	tangent_index.clear();

	nmbpts = 0;
        nmbder = 0;
        for (int kr=0; kr<2; kr++)
          {
            // Treat start- and end point.
	    // Remember that curves form a loop.

//             if ((kr == 0 && ki <= 1) || (kr == 1 && ki >= 2))
//               kj = (ki==0) ? nmbcrvs-1 : ki-1;
//             else
//               kj = (ki==nmbcrvs-1) ? 0 : ki+1;
	      int kj;
	      if (kr == 0)
		  kj = (ki==0) ? nmbcrvs-1 : ki-1;
	      else
		  kj = (ki==nmbcrvs-1) ? 0 : ki+1;

	    // We must satisfy tangent condition.
            // Point condition
//             par = (ki==0 || ki==3) ? boundary_crvs[kj]->startparam()
//               : boundary_crvs[kj]->endparam();
	      par = (kr == 0 ? boundary_crvs[kj]->endparam() :
		     boundary_crvs[kj]->startparam());
	      boundary_crvs[kj]->point(pts, par, 1);
	      //	      point = pts[0]; // Not really used...
	      tangent = pts[1];
	    //            grstat = boundary_crvs[kj].evalPointAndTangent(par, point,
	    //                                             tangent);
	    //            if (grstat != GR_OK)
	    //return grstat;

	      // As curves form a loop, some tangent vectors must be turned.
	      int ksign = (kr == 0 ? -1 : 1);
	      spar[kr] = (kr==0) ? boundary_crvs[ki]->startparam()
		  : boundary_crvs[ki]->endparam();
	      for (int kh = 0; kh < 4; ++kh)
		  knots.push_back(spar[kr]);

	      for (int kh=0; kh<dim; kh++)
		  points[nmbpts*dim+kh] = ksign*tangent[kh];
	      nmbpts++;

	      if (cross_crvs[kj].get() != NULL)
		  {
		      // We must satisfy twist condition.
		      // Derivative condition in start point
		      cross_crvs[kj]->point(pts, par, 1);
		      //     point = pts[0]; // Not that it is used...
		      tangent = pts[1];
		      //  grstat = cross_crvs[kj].evalPointAndTangent(par, point,
		      //                                        tangent);
		      //  if (grstat != GR_OK)
		      //                  return grstat;

//    		      ksign = ((kr==0 && (ki==0 || ki==3)) ||
//    			       (kr==1 && (ki==1 || ki==2))) ? 1 : -1;
    		      ksign = ((kr==0 && (ki==0 || ki==2)) ||
    			       (kr==1 && (ki==1 || ki==3))) ? -1 : 1;
		      for (int kh=0; kh<dim; kh++)
			  //	derivs[nmbder*dim+kh] = ksign*tangent[kh];
			  derivs.push_back(ksign*tangent[kh]);
		      tangent_index.push_back(nmbpts - 1);
		      nmbder++;
		  }
	      else 
		  {
		      --order; // We decrease the order used in the interpolation.
		      --num_coefs;
		  }
          }

	SplineInterpolator interpolator;
	// Using the information in points, we make a suiting basis.
	interpolator.makeBasis(spar, tangent_index, order);
	vector<double> coefs;
	interpolator.interpolate(spar, points, tangent_index, derivs, coefs);

	//	int num_coefs = coefs.size()/boundary_crvs[ki]->dimension();
	shared_ptr<SplineCurve> qc = shared_ptr<SplineCurve>
	    (new SplineCurve(num_coefs, order,
			       interpolator.basis().begin(),
			       coefs.begin(), boundary_crvs[ki]->dimension()));

	// We must raise the curve to the order of the other cross_curves.
	size_t krsz;
	for (krsz = 0; krsz < cross_crvs.size(); ++krsz)
	    if (cross_crvs[krsz].get() != NULL) break;
	if (krsz < cross_crvs.size())
	    qc->raiseOrder(cross_crvs[krsz]->order() - order);

//         // Interpolate curve
//         double endpar;
//         SISLCurve *qc = NULL;
//         double *spar2 = NULL;
//         int nbpar;
//         s1357(sder, nmbder, dim, ltype, spar, 0, 0, 1,
//               boundary_crvs[ki].ik(), spar[0], &endpar, &qc,
//               &spar2, &nbpar, &kstat);

//         if (spar2) free(spar2);
//         if (kstat < 0)
//           return GR_SISLERR;

        cross_crvs[ki] = qc;
	//        cross_crvs[ki].setCurve(qc, false);
      }

  //  return GR_OK;
}



//===========================================================================
void
Go::CoonsPatchGen::makeLoftParams(vector<shared_ptr<SplineCurve> >::const_iterator
				    first_curve,
				    int nmb_crvs, double param_length,
				    vector<double>& params)
//===========================================================================
{
  params.clear();

  // Put the curves into common basis
  double tolerance = 1e-05;
  // @@ sbr Suppose we should implement unify... with iterator arguments.
  vector<shared_ptr<SplineCurve> > the_curves(first_curve, first_curve + nmb_crvs);
  GeometryTools::unifyCurveSplineSpace(the_curves, tolerance);
  // We're further assuming all curves have the same order.

  // Compute parameterization.
  // For each adjacent pair of curves compute the minimum
  // distance between the coefficients
  int dim = first_curve[0]->dimension();
  int num_coefs = first_curve[0]->numCoefs();
  double largenmb = 100000000.0;
  double mindist, dist;
  params.push_back(0.0);  // Parameter value of first curve.
  double currpar;

  int h, i, j;
  for (h = 1; h < nmb_crvs; ++h) {
      mindist = largenmb;
      for (i = 0; i < num_coefs; ++i) {
	  dist = 0;
	  for (j = 0; j < dim; ++j)
	      dist += (first_curve[h]->coefs_begin()[i*dim+j] -
		  first_curve[h-1]->coefs_begin()[i*dim+j]) *
		  (first_curve[h]->coefs_begin()[i*dim+j] -
		   first_curve[h-1]->coefs_begin()[i*dim+j]);
	  mindist = std::min(sqrt(dist), mindist);
      }
      currpar = params[h-1] + mindist;
      params.push_back(currpar);
  }

  // We make sure the parameters go from 0.0 to param_length.
  double frac = param_length / params[nmb_crvs-1];
  for (h = 1; h < nmb_crvs; ++h)
      params[h] *= frac;
}


// void
//       sh1260(double aconst,SISLCurve *vcurve[],int icurve,int *jstat)
void
Go::CoonsPatchGen::reparamBoundaryCurve(vector<shared_ptr<SplineCurve> >& curves,
					  double aconst)
/*
*********************************************************************
*                                                                   
* PURPOSE    : Check length of tangent vectors at the endpoints of a
*              curve compared to the size of the curve. If the vectors
*              are too long, reparametrize the curve and a number of
*              corresponding curves.
*
*
*
* INPUT      : aconst     - Constant used to check when the tangent
*                           vectors at the endpoints of the curve are
*                           too long.
*              icurve     - Number of curves. icurve >= 1.
*              
*
* INPUT/OUTPUT : vcurve   - Array containing a curve set. The curves
*                           are expected to have the same parametrization.
*                           If the tangents in the endpoints of the
*                           first curve is too long, all curves are
*                           reparametrized. Dimension of the array is icurve.
*                       
*
* OUTPUT     : jstat      - status messages  
*                                         > 0      : warning
*                                         = 0      : ok
*                                         < 0      : error
*
*
* METHOD     : 
*
*
* REFERENCES : 
*              
*
* USE        : 3D geometry only.
*
*-
* CALLS      : s1221    - Evaluate curve.  
*              s6diff   - Difference vector between two vectors.  
*              s6scpr   - Scalar product between two vectors.  
*              s6length - Length of vector.   
*              
*
* WRITTEN BY : Vibeke Skytt, SI, 06.90.
*
*********************************************************************
*/
{
    // Testing input.
    int icurve = (int)curves.size();
    ALWAYS_ERROR_IF(icurve == 0,
		"Input vector was empty!");
    int kdim = curves[0]->dimension();/* Dimension of first curve in the curve set.*/
    double tpar1 = curves[0]->startparam(); /*Start of param interval of 1.curve.*/
    double tpar2 = curves[0]->endparam(); /* End of param interval of 1. curve. */
    int ki;             /* Counter.           */
    for (ki=1; ki<icurve; ki++) {
	ALWAYS_ERROR_IF(curves[ki]->dimension() != kdim,
		    "Dimension of curve spaces differ.");
	ALWAYS_ERROR_IF((curves[ki]->startparam() != tpar1) ||
		    (curves[ki]->endparam() != tpar2),
			"Input curves not defined over same parameter interval.");

    }

    int kder = 1;       /* Number of derivatives of curve to evaluate. */
    //  double *sder1 = SISL_NULL;  /* Value of 1.curve in start of param interval.*/
    vector<Point> sder1(kder+1);
    vector<Point> sder2(kder+1);
    Point sdiff;
    double tdiff;          /* Length of sdiff.  */
    double t1,t2;          /* Length of the components of the tangent vectors in the 
                            endpoints along sdiff, compared with length of sdiff.*/
    double tscal;          /* The factor with which to scale the curves.  */
    double tnewend;        /* New endpoint of parameter interval.         */
    vector<double>::iterator s1; /* Pointer used to traverse knot vector of curve.*/
    vector<double>::iterator s2; /* Pointer used to stop traversing knot vector.  */
    //    shared_ptr<SplineCurve> qcpt;        /* Pointer to curve in curve set. */

    
  
    /* Evaluate the first curve in the endpoints.  */
    curves[0]->point(sder1, tpar1, 1);
    curves[0]->point(sder2, tpar2, 1);

    /* Compute difference vector between endpoints of position curve. */
    sdiff = sder2[0] - sder1[0];

    /* Compute length of difference vector.  */
    tdiff = sdiff.length();
  
    /* Compute length of the component of the tangent in the first endpoint
       along the difference vector compared to the distance between the
       endpoints. The result lies between 0 and 1.  */
    t1 = sder1[1]*sdiff / (tdiff*tdiff);
  
    /* Compute length of the component of the tangent in the second endpoint
       along the difference vector compared to the distance between the
       endpoints. The result lies between 0 and 1.  */
    t2 = sder2[1]*sdiff / (tdiff*tdiff);
  
    /* Check if any of the tangents are too long.  */
    if (std::max(t1,t2) > aconst) {
	/* One of the tangents is too long. Reparametrize to reduce tangent
	   length.   */

	tscal = std::max(t1,t2)/aconst;
      
	/* Find new endpoint of parameter interval. The startpoint is kept. */

	tnewend = tscal*(tpar2-tpar1) + tpar1;

	for (ki=0; ki<icurve; ki++) {
	    vector<double> new_knots(curves[ki]->basis().begin(),
				 curves[ki]->basis().end());

	    //	    qcpt = curves[ki];
	  
	    /* Traverse position and u- and v-derivative curves.  */

// 	    for (s1=qcpt->basis().begin(),
// 		     s2=qcpt->basis().begin() + qcpt->numCoefs() +
// 		     qcpt->order();
	    for (s1 = new_knots.begin(), s2 = new_knots.end(); s1 < s2; s1++)
		*s1 = tpar1 + (*s1 - tpar1)*(tnewend - tpar1)/(tpar2 - tpar1);

	    *curves[ki] = SplineCurve(curves[ki]->numCoefs(), curves[ki]->order(),
					new_knots.begin(), curves[ki]->coefs_begin(),
					curves[ki]->dimension());
	}

    }

    return;
}     


namespace {


//===========================================================================
void
blendder(double ea[],double eb[],int ix,int ieq,
	 int irang,double astarti,double aendi,
	 double astartj,double aendj,double ealfai[],
	 double ebetai[],int nmbcoef1,double ealfaj[],
	 double ebetaj[],int nmbcoef2)
//===========================================================================
/*
*********************************************************************
*                                                                   
* PURPOSE    : Estimate the derivatives in one endpoint of the blending 
*              functions meeting in a corner of a vertex region.
*
*
* INPUT      : ea         - Matrix containg coefficients in equation
*                           system representing conditions on the 
*                           derivatives.
*              eb         - Rigth side of equation system.
*              ix         - Number of coefficients. ix = 4.
*              ieq        - Number of equations. ieq = 3.
*              irang      - Rang of equation system.
*              nmbcoef1   - Number of coefficients of blending functions
*                           along first adjacent edge.
*              nmbcoef2   - Number of coefficients of blending functions
*                           along second adjacent edge.
*
*
* INPUT/OUTPUT : ealfai   - Position and derivatives of first blending 
*                           function along first adjacent edge.
*                           Estimate derivative of first endpoint.
*                ebetai   - Position and derivatives or second blending 
*                           function along first adjacent edge.
*                           Estimate derivative of first endpoint.
*                ealfaj   - Position and derivatives or first blending 
*                           function along second adjacent edge.
*                           Estimate derivative of second endpoint.
*                ebetaj   - Position and derivatives or second blending 
*                           function along second adjacent edge.
*                           Estimate derivative of second endpoint.
*                
* 
* OUTPUT     : jstat      - status messages  
*                                         > 0      : warning
*                                         = 0      : ok
*                                         < 0      : error
*
*
* METHOD     :
*
*
* REFERENCES : 
*              
*
* USE        : 
*
*-
* CALLS      : s6lufacp - LU-factorizing of matrix in equation system. 
*              s6lusolp - Solve equtaion system. Matrix is LU-factorized. 
*              
*
* WRITTEN BY : Vibeke Skytt, SI, 06.90.
*
*********************************************************************
*/
{
  double epsco = 1e-14;
  double REL_PAR_RES = 1e-6;

  // int kstat = 0;          /* Status variable.  */
  int ki,kj,kk,kh;        /* Counters.         */
  int l[2][2];            /* Pointer to lines of equation system to be used
                             when solving with respect to two unknowns.   */
  int kcount;             /* Number of iterations used to find derivatives. */
  double teps = 0.001;    /* Local tolerance.                             */
  double tdiff;           /* Difference between derivatives of linear blending
                             function and current blending function.      */
  double tdum;            /* Used in testing of result.                   */
  double tdet;            /* Current determinant of equation system.      */
  double sdet[2];         /* Determinant of equation system.              */
  double sx[4];           /* Derivative of linear blending funtion.       */
  double sxx[4];          /* Derivative of current cubic blending funtion. */
  double sb[3];           /* Rigth side of local equation system.         */
  
  /* Test input.  */

  if (ix != 4 || ieq != 3) 
    {
	THROW("Wrong dimension of matrix system.");
	//      return;
    }
    
  /* Initialize output. */
  //  *jstat = 0;
  
  /* Set up derivatives of linear blend.  */

  sx[0] = (ealfai[nmbcoef1-1] - ealfai[0])/(aendi - astarti);
  sx[1] = (ebetai[nmbcoef1-1] - ebetai[0])/(aendi - astarti);
  sx[2] = (ealfaj[nmbcoef2-1] - ealfaj[0])/(aendj - astartj);
  sx[3] = (ebetaj[nmbcoef2-1] - ebetaj[0])/(aendj - astartj); 
  
  if (irang == 1)
    {
      /* It is not possible to solve the equation system. Find the
	 best possible solution. */
	vector<vector<double> > smat2(4);
	for (int ii = 0; ii < 4; smat2[ii++].resize(4));
	vector<double> sb2(4);
      for (ki=0; ki<ix; ki++) {
	  sb2[ki] = 0.0;
	  for (kk=0; kk<ieq; kk++)
	      sb2[ki] += ea[kk*ix+ki]*eb[kk];

	  for (kj=0; kj<ix; kj++)
	    {
	      smat2[ki][kj] = 0.0;
	      for (kk=0; kk<ieq; kk++)
		smat2[ki][kj] += ea[kk*ix+ki]*ea[kk*ix+kj];
	    }
      }
      LUsolveSystem(smat2, 4, &sb2[0]);
      for (ki = 0; ki < ix; ++ki) {
	  sxx[ki] = sb2[ki];
      }
    }
  if (irang == 2)
    {
      /* Rang of equation system equal to 2. */
      /* Find equations to be used, i.e. find the lines in the
	 equation system which gives the most stable 2x2 equation
	 system when finding the 2 last derivatives when the 2
	 first are set.  */

      sdet[0] = sdet[1] = (double)0.0;
      
      for (ki=0; ki<ieq; ki++)
	for (kj=ki+1; kj<ieq; kj++)
	  for (kk=0,kh=0; kh<2; kk+=2,kh++)
	    {
	      tdet = ea[4*ki+kk]*ea[4*kj+kk+1] - ea[4*ki+kk+1]*ea[4*kj+kk];
	      if (fabs(tdet) > fabs(sdet[kh]))
		{
		  sdet[kh] = tdet;
		  l[kh][0] = ki, l[kh][1] = kj;
		}
	    }
      for (kh=0; kh<2; kh++)
	if (fabs(sdet[kh])<epsco) 
	  sdet[kh] = (double)1.0;
      
      kcount = 0;
      
      /* Set the two first derivatives equal to those of a corresponding
	 linear blending function.   */

      ki = 0, kk = 2;
      sxx[ki] = sx[ki];
      sxx[ki+1] = sx[ki+1];
      
      while (1)
	{
	  /* Iterate until the conditions on the derivatives is satisfied
	     and the difference from the linear case is eqal for all 
	     derivatives or the number of iterations is equal to 4.   */

	  /* Find right side of local equation system for finding the 
	     remaining two derivatives.   */

	  for (kj=0; kj<ieq; kj++)
	    sb[kj] = eb[kj] - ea[4*kj+ki]*sxx[ki] - ea[4*kj+ki+1]*sxx[ki+1];
	  
	  /* Compute the remaining two derivatives.  */

	  kj = kk/2;
	  sxx[kk] = (-ea[4*l[kj][0]+kk+1]*sb[l[kj][1]] 
		       + ea[4*l[kj][1]+kk+1]*sb[l[kj][0]])/sdet[kj];
	  sxx[kk+1] = (ea[4*l[kj][0]+kk]*sb[l[kj][1]] 
		       - ea[4*l[kj][1]+kk]*sb[l[kj][0]])/sdet[kj];

	  
	  /* Check if the difference between current derivatives and 
	     derivatives of a linear blend is small.  */

	  tdiff = fabs(sxx[0] - sx[0]);
	  for (kh=1; kh<4; kh++)
	    if (fabs(tdiff-fabs(sxx[kh]-sx[kh])) > teps) break;
	  
	  if (kh == 4 || kcount == 4) break;   /* Stop iteration.  */
	  
	  /* Initiate next iteration step.  */

	  sxx[kk] = (sx[kk] + sxx[kk])/(double)2.0;
	  sxx[kk+1] = (sx[kk+1] + sxx[kk+1])/(double)2.0;
	  
	  ki = (ki + 2) % 4;
	  kk = (kk + 2) % 4;

	  kcount++;
	}
    }
  else if (irang == 3)
    {
	vector<vector<double> > smat3(3);
	for (int ii = 0; ii < 3; smat3[ii++].resize(3));
	vector<double> sb3(3);

	/* Rang of equation system equal to 3. */
	ki = 0;
	kcount = 0;
      
	while (1) {
	  /* Iterate until the conditions on the derivatives is satisfied
	     and the difference from the linear case is eqal for all 
	     derivatives or the number of iterations is equal to 4.   */

	  /* Set the current derivative equal to those of a corresponding
	     linear blending function.   */

	  sxx[ki] = sx[ki];
	  
	  /* Set up local equation system.  */

	  for (kj=0; kj<ieq; kj++)
	    {
	      sb3[kj] = eb[kj] - ea[kj*4+ki]*sxx[ki];
	      for (kk=0,kh=0; kk<ix; kk++)
		{
		  if (kk == ki) continue;
		  smat3[kj][kh] = ea[kj*ix+kk];
		  kh++;
		}
	    }
	  
	  /* Compute the remaining three derivatives.  */
	  LUsolveSystem(smat3, (int)smat3.size(), &sb3[0]);
	  for (kk=0,kh=0; kk<ix; kk++)
	    {
	      if (kk == ki) continue;
	      sxx[kk] = sb3[kh++];
	    }
	  
	  /* Check if the difference between current derivatives and 
	     derivatives of a linear blend is small.  */

	  tdiff = fabs(sxx[0] - sx[0]);
	  for (kh=1; kh<4; kh++)
	    if (fabs(tdiff-fabs(sxx[kh]-sx[kh])) > teps) break;
	  
	  if (kh == 4 || kcount == 4) break;   /* Stop iteration.  */
	  
	  /* Initiate next iteration step.  */

	  kj = (ki + 1) % 4;
	  sxx[kj] = (sx[kj] + sxx[kj])/(double)2.0;
	  
	  ki = kj;
	  kcount++;
	}
      
    }
  
  /* Copy derivatives into output array.  */ 

  ealfai[1] = sxx[0];
  ebetai[1] = sxx[1];
  ealfaj[nmbcoef2-2] = sxx[2];
  ebetaj[nmbcoef2-2] = sxx[3];
  
  /* Test result.  */

  for (ki=0; ki<ieq; ki++)
    {
      tdum = (double)0.0;
      for (kj=0; kj<ix; kj++)
	tdum += ea[ki*ix+kj]*sxx[kj];
      if (fabs(tdum - eb[ki]) > REL_PAR_RES) break;
    }
//    if (ki < ieq) 
//      /* Twist requirement not satisfied completely.  */
//        //    *jstat = 1;
//        MESSAGE("Twist requirement not satisfied completely.");

  return;
}

// the below implementation depended on newmat. Use above implementation instead

// //===========================================================================
// void
// blendder(double ea[],double eb[],int ix,int ieq,
// 	 int irang,double astarti,double aendi,
// 	 double astartj,double aendj,double ealfai[],
// 	 double ebetai[],int nmbcoef1,double ealfaj[],
// 	 double ebetaj[],int nmbcoef2)
// //===========================================================================
// /*
// *********************************************************************
// *                                                                   
// * PURPOSE    : Estimate the derivatives in one endpoint of the blending 
// *              functions meeting in a corner of a vertex region.
// *
// *
// * INPUT      : ea         - Matrix containg coefficients in equation
// *                           system representing conditions on the 
// *                           derivatives.
// *              eb         - Rigth side of equation system.
// *              ix         - Number of coefficients. ix = 4.
// *              ieq        - Number of equations. ieq = 3.
// *              irang      - Rang of equation system.
// *              nmbcoef1   - Number of coefficients of blending functions
// *                           along first adjacent edge.
// *              nmbcoef2   - Number of coefficients of blending functions
// *                           along second adjacent edge.
// *
// *
// * INPUT/OUTPUT : ealfai   - Position and derivatives of first blending 
// *                           function along first adjacent edge.
// *                           Estimate derivative of first endpoint.
// *                ebetai   - Position and derivatives or second blending 
// *                           function along first adjacent edge.
// *                           Estimate derivative of first endpoint.
// *                ealfaj   - Position and derivatives or first blending 
// *                           function along second adjacent edge.
// *                           Estimate derivative of second endpoint.
// *                ebetaj   - Position and derivatives or second blending 
// *                           function along second adjacent edge.
// *                           Estimate derivative of second endpoint.
// *                
// * 
// * OUTPUT     : jstat      - status messages  
// *                                         > 0      : warning
// *                                         = 0      : ok
// *                                         < 0      : error
// *
// *
// * METHOD     :
// *
// *
// * REFERENCES : 
// *              
// *
// * USE        : 
// *
// *-
// * CALLS      : s6lufacp - LU-factorizing of matrix in equation system. 
// *              s6lusolp - Solve equtaion system. Matrix is LU-factorized. 
// *              
// *
// * WRITTEN BY : Vibeke Skytt, SI, 06.90.
// *
// *********************************************************************
// */
// {
//   double epsco = 1e-14;
//   double REL_PAR_RES = 1e-6;

//   // int kstat = 0;          /* Status variable.  */
//   int ki,kj,kk,kh;        /* Counters.         */
//   int l[2][2];            /* Pointer to lines of equation system to be used
//                              when solving with respect to two unknowns.   */
//   int kcount;             /* Number of iterations used to find derivatives. */
//   double teps = 0.001;    /* Local tolerance.                             */
//   double tdiff;           /* Difference between derivatives of linear blending
//                              function and current blending function.      */
//   double tdum;            /* Used in testing of result.                   */
//   double tdet;            /* Current determinant of equation system.      */
//   double sdet[2];         /* Determinant of equation system.              */
//   double sx[4];           /* Derivative of linear blending funtion.       */
//   double sxx[4];          /* Derivative of current cubic blending funtion. */
//   double sb[3];           /* Rigth side of local equation system.         */
  
//   /* Test input.  */

//   if (ix != 4 || ieq != 3) 
//     {
// 	THROW("Wrong dimension of matrix system.");
// 	//      return;
//     }
    
//   /* Initialize output. */
//   //  *jstat = 0;
  
//   /* Set up derivatives of linear blend.  */

//   sx[0] = (ealfai[nmbcoef1-1] - ealfai[0])/(aendi - astarti);
//   sx[1] = (ebetai[nmbcoef1-1] - ebetai[0])/(aendi - astarti);
//   sx[2] = (ealfaj[nmbcoef2-1] - ealfaj[0])/(aendj - astartj);
//   sx[3] = (ebetaj[nmbcoef2-1] - ebetaj[0])/(aendj - astartj); 
  
//   if (irang == 1)
//     {
//       /* It is not possible to solve the equation system. Find the
// 	 best possible solution. */
//       Matrix smat2(4,4);
//       ColumnVector sb2(4);
//       for (ki=0; ki<ix; ki++)
// 	{
// 	  sb2.element(ki) = 0.0;
// 	  for (kk=0; kk<ieq; kk++)
// 	    sb2.element(ki) += ea[kk*ix+ki]*eb[kk];

// 	  for (kj=0; kj<ix; kj++)
// 	    {
// 	      smat2.element(ki, kj) = 0.0;
// 	      for (kk=0; kk<ieq; kk++)
// 		smat2.element(ki, kj) += ea[kk*ix+ki]*ea[kk*ix+kj];
// 	    }
// 	}
// //       s6lufacp(smat2, lpiv, ix, &kstat);
// //       if (kstat < 0) 
// // 	{
// // 	  *jstat = kstat;
// // 	  return;
// // 	}

// //       s6lusolp(smat2, sb2, lpiv, ix, &kstat);
// //       if (kstat < 0) 
// // 	{
// // 	  *jstat = kstat;
// // 	  return;
// // 	}
//       CroutMatrix smat2LUfact = smat2;
//       ALWAYS_ERROR_IF(smat2LUfact.IsSingular(),
// 		      "System couldn't be solved; should not happen!");

//       ColumnVector x = smat2LUfact.i() * sb2;

//       for (ki=0; ki<ix; ki++)
// 	sxx[ki] = x.element(ki);
//     }
//   if (irang == 2)
//     {
//       /* Rang of equation system equal to 2. */
//       /* Find equations to be used, i.e. find the lines in the
// 	 equation system which gives the most stable 2x2 equation
// 	 system when finding the 2 last derivatives when the 2
// 	 first are set.  */

//       sdet[0] = sdet[1] = (double)0.0;
      
//       for (ki=0; ki<ieq; ki++)
// 	for (kj=ki+1; kj<ieq; kj++)
// 	  for (kk=0,kh=0; kh<2; kk+=2,kh++)
// 	    {
// 	      tdet = ea[4*ki+kk]*ea[4*kj+kk+1] - ea[4*ki+kk+1]*ea[4*kj+kk];
// 	      if (fabs(tdet) > fabs(sdet[kh]))
// 		{
// 		  sdet[kh] = tdet;
// 		  l[kh][0] = ki, l[kh][1] = kj;
// 		}
// 	    }
//       for (kh=0; kh<2; kh++)
// 	if (fabs(sdet[kh])<epsco) 
// 	  sdet[kh] = (double)1.0;
      
//       kcount = 0;
      
//       /* Set the two first derivatives equal to those of a corresponding
// 	 linear blending function.   */

//       ki = 0, kk = 2;
//       sxx[ki] = sx[ki];
//       sxx[ki+1] = sx[ki+1];
      
//       while (1)
// 	{
// 	  /* Iterate until the conditions on the derivatives is satisfied
// 	     and the difference from the linear case is eqal for all 
// 	     derivatives or the number of iterations is equal to 4.   */

// 	  /* Find right side of local equation system for finding the 
// 	     remaining two derivatives.   */

// 	  for (kj=0; kj<ieq; kj++)
// 	    sb[kj] = eb[kj] - ea[4*kj+ki]*sxx[ki] - ea[4*kj+ki+1]*sxx[ki+1];
	  
// 	  /* Compute the remaining two derivatives.  */

// 	  kj = kk/2;
// 	  sxx[kk] = (-ea[4*l[kj][0]+kk+1]*sb[l[kj][1]] 
// 		       + ea[4*l[kj][1]+kk+1]*sb[l[kj][0]])/sdet[kj];
// 	  sxx[kk+1] = (ea[4*l[kj][0]+kk]*sb[l[kj][1]] 
// 		       - ea[4*l[kj][1]+kk]*sb[l[kj][0]])/sdet[kj];

	  
// 	  /* Check if the difference between current derivatives and 
// 	     derivatives of a linear blend is small.  */

// 	  tdiff = fabs(sxx[0] - sx[0]);
// 	  for (kh=1; kh<4; kh++)
// 	    if (fabs(tdiff-fabs(sxx[kh]-sx[kh])) > teps) break;
	  
// 	  if (kh == 4 || kcount == 4) break;   /* Stop iteration.  */
	  
// 	  /* Initiate next iteration step.  */

// 	  sxx[kk] = (sx[kk] + sxx[kk])/(double)2.0;
// 	  sxx[kk+1] = (sx[kk+1] + sxx[kk+1])/(double)2.0;
	  
// 	  ki = (ki + 2) % 4;
// 	  kk = (kk + 2) % 4;

// 	  kcount++;
// 	}
//     }
//   else if (irang == 3)
//     {

//       Matrix smat3(3,3);
//       ColumnVector sb3(3);

//       /* Rang of equation system equal to 3. */
//       ki = 0;
//       kcount = 0;
      
//       while (1)
// 	{
// 	  /* Iterate until the conditions on the derivatives is satisfied
// 	     and the difference from the linear case is eqal for all 
// 	     derivatives or the number of iterations is equal to 4.   */

// 	  /* Set the current derivative equal to those of a corresponding
// 	     linear blending function.   */

// 	  sxx[ki] = sx[ki];
	  
// 	  /* Set up local equation system.  */

// 	  for (kj=0; kj<ieq; kj++)
// 	    {
// 	      sb3.element(kj) = eb[kj] - ea[kj*4+ki]*sxx[ki];
// 	      for (kk=0,kh=0; kk<ix; kk++)
// 		{
// 		  if (kk == ki) continue;
// 		  smat3.element(kj, kh) = ea[kj*ix+kk];
// 		  kh++;
// 		}
// 	    }
	  
// 	  /* Compute the remaining three derivatives.  */

// // 	  s6lufacp(smat,lpiv,ieq,&kstat);
// // 	  if (kstat < 0) 
// // 	    {
// // 	      *jstat = kstat;
// // 	      return;
// // 	    }
	  
// // 	  s6lusolp(smat,sb,lpiv,ieq,&kstat);
// // 	  if (kstat < 0) 
// // 	    {
// // 	      *jstat = kstat;
// // 	      return;
// // 	    }
// 	  CroutMatrix smat3LUfact = smat3;
// 	  ALWAYS_ERROR_IF(smat3LUfact.IsSingular(),
// 			  "System couldn't be solved; should not happen!");

// 	  ColumnVector x = smat3LUfact.i() * sb3;

	  
// 	  for (kk=0,kh=0; kk<ix; kk++)
// 	    {
// 	      if (kk == ki) continue;
// 	      sxx[kk] = x.element(kh++);
// 	    }
	  
// 	  /* Check if the difference between current derivatives and 
// 	     derivatives of a linear blend is small.  */

// 	  tdiff = fabs(sxx[0] - sx[0]);
// 	  for (kh=1; kh<4; kh++)
// 	    if (fabs(tdiff-fabs(sxx[kh]-sx[kh])) > teps) break;
	  
// 	  if (kh == 4 || kcount == 4) break;   /* Stop iteration.  */
	  
// 	  /* Initiate next iteration step.  */

// 	  kj = (ki + 1) % 4;
// 	  sxx[kj] = (sx[kj] + sxx[kj])/(double)2.0;
	  
// 	  ki = kj;
// 	  kcount++;
// 	}
      
//     }
  
//   /* Copy derivatives into output array.  */ 

//   ealfai[1] = sxx[0];
//   ebetai[1] = sxx[1];
//   ealfaj[nmbcoef2-2] = sxx[2];
//   ebetaj[nmbcoef2-2] = sxx[3];
  
//   /* Test result.  */

//   for (ki=0; ki<ieq; ki++)
//     {
//       tdum = (double)0.0;
//       for (kj=0; kj<ix; kj++)
// 	tdum += ea[ki*ix+kj]*sxx[kj];
//       if (fabs(tdum - eb[ki]) > REL_PAR_RES) break;
//     }
// //    if (ki < ieq) 
// //      /* Twist requirement not satisfied completely.  */
// //        //    *jstat = 1;
// //        MESSAGE("Twist requirement not satisfied completely.");

//   return;
// }


//===========================================================================
void
checkContinuity(const SplineCurve* pos,
		const SplineCurve* der1, 
		const SplineCurve* der2,
		const SplineCurve* cross,
		double& minAng, double& medAng, 
		double& maxAng)
//===========================================================================
{
  double accuracy = 0.01;
  double pihalf = 3.14159/2.0;


  // Initiate output
  minAng = 100000000.0;
  maxAng = medAng = 0.0;

  // Estimate curve length (total distance between spline coefs).
  double crv_length = 0;
  double length;
  int dim = pos->dimension();
  int i, j;
  for (i = 0; i < pos->numCoefs() - 1; ++i) {
      length = 0.0;
      for (j = 0; j < dim; ++j)
	  length += pos->coefs_begin()[(i + 1) * dim + j] -
	      pos->coefs_begin()[i * dim + j];
      crv_length += length;
  }


  int nmbStep = (int)(crv_length/accuracy);
  nmbStep = std::max(nmbStep, 2*pos->numCoefs());
  nmbStep = std::min(nmbStep, 500);
  
  double startPar = pos->startparam();
  double incre = (pos->endparam() - startPar)/(double)(nmbStep - 1);

  // March along the edge curve, checking continuity between adjacent
  // surfaces.
  int ki;
  double tpar;
  double ang;
  Point crosspt, derpt1, derpt2, norm;
  for (ki=0, tpar=startPar; ki<nmbStep; ki++, tpar+=incre)
    {
      // Evaluate the curves
      der1->point(derpt1, tpar);
      der2->point(derpt2, tpar);
      cross->point(crosspt, tpar);
      norm = derpt1 % derpt2;

      // Returning euclidean norm.
      ang = norm.angle(crosspt);
      ang = fabs(pihalf-ang);

      minAng = std::min(ang, minAng);
      maxAng = std::max(ang, maxAng);
      medAng += ang;
    }

  // Compute medium values
  medAng /= nmbStep;


  //  return GR_OK;
}


}; // end anonymous namespace
