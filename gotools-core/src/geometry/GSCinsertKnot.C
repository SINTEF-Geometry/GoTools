
//===========================================================================
//                                                                           
// File: GSCinsertKnot.C                                                     
//                                                                           
// Created: Mon Jan 22 13:17:19 2001                                         
//                                                                           
// Author: Atgeirr F Rasmussen <atgeirr@sintef.no>
//                                                                           
// Revision: $Id: GSCinsertKnot.C,v 1.14 2005-06-27 15:32:58 oan Exp $
//                                                                           
// Description:
//                                                                           
//===========================================================================


#include "GoTools/geometry/SplineCurve.h"
#include "GoTools/geometry/SplineUtils.h"
using namespace Go;

//===========================================================================
void SplineCurve::insertKnot(double apar)
//===========================================================================

/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : Insert a given knot into the description of
*              a B-spline curve.
* NOTE       : When the curve is periodic, the input parameter value
*              must lie in the HALFOPEN [et[kk-1], et[kn), the function
*              will automatically update the extra knots and
*              coeffisients.
*              rcnew->in is still eq to pc->in + 1!
*
*
* INPUT      : pc        - SISLCurve to be refined.
*              apar      - Parameter values of knot to be s1017ed.
*
*
*
* OUTPUT     : rc        - The new, refined curve.
*              jstat     - status messages
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
*-
* CALLS      : newCurve  - Allocate space for a new curve-object.
*              freeCurve - Free space occupied by given curve-object.
*              S1701.C   - Making the knot-s1017en-transformation matrix.
*              s1017knots in periodic case.
* WRITTEN BY : Arne Laksaa, SI, 88-11.
* CHANGED BY : Ulf J. Krystad, SI, 92-01
*              Treatment of periodic curves.
* Revised by : Paal Fugelli, SINTEF, Oslo, Norway, Nov. 1994.  Undefined
*              symbol 's1017knots()' changed to 's1018()'. Closed ('cuopen'
*              flag=0) and 'ikind' flag are now handled correctly. Fixed memory
*              problem for rationals.  Cleaned up a bit.
*
**********************************************************************/
{

    //  int kstat;			/* Local status variable.                     */
    //  int kpos = 0;			/* Position of error.                         */
  int kmy;			/* An index to the knot-vector.               */
  int kpl, kfi, kla;		/* To posisjon elements in trans.-matrix.     */
  int kk = order();		/* Order of the input curve.                  */
  int kn = numCoefs();		/* Number of the vertices in input curve.     */
  int kdim = dimension();	/* Dimension of the space in whice curve lies.*/
  if (rational_)
      ++kdim;
  int kn1;			/* Number of vertices in the new curve.       */
  int kch;			/* First vertice to be changes.               */
  int knum = 0;			/* Number of knots less and equal than
                                   the intersection point.                    */
  int ki;			/* Control variable in loop.                  */
  int kj, kj1, kj2;		/* Control variable in loop.                  */
  double *s1;           	/* Pointers used in loop.                     */
//    double *st = NULL;		/* The first new knot-vector.                 */
//    double *salfa = NULL;	/* A line of the trans.-matrix.               */
//    double *scoef = NULL;		/* The first new vertice.                     */
//    double *ecoef;                /* Pointer to input curve ecoef/rcoef vector. */
//    SISLCurve *qc = NULL;	/* Pointer to new curve-object.               */



  /* Make sure the returned pointer is valid. */

//    *rc = NULL;


  /* Check that we have a curve. */

//    if (!pc)
//      goto err150;


  /* Find the number of vertices in the new curve. */

  kn1 = kn + 1;

//    if ( kn1 <= 0 )
//      goto err150;



  /* Periodicity treatment -------------------------- */
//    if (pc->cuopen == SISL_CRV_PERIODIC)
//      {
//        s1018(pc, &apar, 1, rc, &kstat);
//        if (kstat < 0)
//  	goto err153;
//        goto out;
//      }



  /* Check that the intersection point is an interior point. */

//    if (apar < *(pc->et) || apar > *(pc->et + kn + kk - 1))
//      goto err158;
  // May change apar, if it is very close to a knot.
  (void)basis_.knotIntervalFuzzy(apar);
  std::vector<double>::iterator orig_knot = basis_.begin();
  if (apar < orig_knot[0] || apar > orig_knot[kn + kk - 1]) {
      THROW("New knot outside knot vector interval.");
  }


  /* Check if the curve is rational.  */

  /* Allocate space for the kk elements which may not be zero in eache
     line of the basic transformation matrix and for a help array of
     kk elements.*/

  std::vector<double> salfa(2*kk);


  /* Find the number of the knots which is smaller or like
     the intersection point.*/

  s1 = &orig_knot[0];

  if (apar > orig_knot[0] && apar < orig_knot[kn + kk - 1])
    {
      /* Using binary search*/
      kj1 = 0;
      kj2 = kk + kn - 1;
      knum = (kj1 + kj2) / 2;
      while (knum != kj1)
	{
	  if (s1[knum] < apar)
	    kj1 = knum;
	  else
	    kj2 = knum;
	  knum = (kj1 + kj2) / 2;
	}
      knum++;			/* The smaller knots. */

      while (s1[knum] == apar)
	/* The knots that are equal to the intersection point. */
	knum++;
    }
  else if (apar == orig_knot[0])
    {
      knum = 0;
      while (s1[knum] == apar)
	/* The knots that are equal to the intersection point. */
	knum++;
    }
  else if (apar == orig_knot[kn + kk - 1])
    {
      knum = 0;
      while (s1[knum] < apar)
	/* The knots that are less than or equal to the intersection point. */
	knum++;
    }



  /* Allocating the new arrays to the new curve. */

  std::vector<double> scoef(kn1 * kdim);
  std::vector<double> st(kn1 + kk);


  /* Copying the knotvectors, all but the intersection point from
     the old curve to the new curves */

  std::copy(orig_knot, orig_knot + knum, st.begin());
  st[knum] = apar;
  if (knum < kn + kk)
      std::copy (orig_knot + knum, orig_knot + kn + kk, st.begin() + knum + 1);


  /* Copying the coefisientvector to the new curve.  (Here 'ecoef' points
     to 'pc->rcoef' or 'pc->ecoef' depening on if 'pc' is rational or not). */

  kch = knum - kk + 1;
  std::vector<double>& co = rational_ ? rcoefs_ : coefs_;

  if (kch > 0)
      std::copy(co.begin(), co.begin() + kdim*kch, scoef.begin());
  if (knum < kn1)
      std::copy(co.begin() + kdim*(knum-1),
		co.begin() + kdim*(kn1-1),
		scoef.begin() + kdim*knum);


  /* Updating the coefisientvectors the new curve.*/

  /* Updating the first curve. */

  for (ki = std::max(0, kch), kmy = 0, s1 = &scoef[0] + ki * kdim;
       ki<std::min(knum + 1, kn1);
       ki++)
    {
      /* Initialising:
           ki = kch,        Index of the vertices we are going to
                             change. Starting with kch, but if
                             kch is negativ we start at zero.
           s1=scoef1+ki*kdim,Pointer at the first vertice to
                             change. */


      /* Using the Oslo-algorithm to make a transformation-vector
         from the old vertices to one new vertice. */

      while (kmy < kn + kk && orig_knot[kmy] <= st[ki])
	kmy++;

      Go::osloalg (ki, kmy - 1, kk, kn, &kpl, &kfi, &kla,
	     &st[0], &orig_knot[0], &salfa[0]);

      /* Compute the kdim vertices with the same "index". */

      for (kj = 0; kj < kdim; kj++, s1++)
	for (*s1 = 0, kj1 = kfi, kj2 = kfi + kpl; kj1 <= kla; kj1++, kj2++)
	  *s1 += salfa[kj2] * co[kj1 * kdim + kj];
    }



  /* Allocating new curve-objects.*/


    co.swap(scoef);
    if (rational_) {
	updateCoefsFromRcoefs();
    }
    basis_.insertKnot(apar);
  
}


//===========================================================================
void SplineCurve::insertKnot(const std::vector<double>& new_knots)
//===========================================================================
{
    // @@ This could be optimized a lot!
    for (size_t i = 0; i < new_knots.size(); ++i) {
	insertKnot(new_knots[i]);
    }
}
