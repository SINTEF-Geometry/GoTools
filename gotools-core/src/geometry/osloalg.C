//===========================================================================
//                                                                           
// File: osloalg.C                                                           
//                                                                           
// Created: Fri Apr  6 14:23:36 2001                                         
//                                                                           
// Author: Atgeirr F Rasmussen <atgeirr@sintef.no>
//                                                                           
// Revision: $Id: osloalg.C,v 1.7 2003-05-08 15:37:14 afr Exp $
//                                                                           
// Description:
//                                                                           
//===========================================================================

#include "GoTools/geometry/SplineUtils.h"
#include <vector>
//using namespace Go;

void 
Go::osloalg(int ij,int imy,int ik,int in,int *jpl,int *jfi,int *jla,
	    double *et,double *etau,double *galfa)

/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : To compute in a compact format a line in the discrete
*              B-spline matrix converting between an orginal basis
*              "etau" and a new basis "et".
*
*
*
* INPUT    : ij     - The index of the new vertice
*            imy    - An index on etau, where the input value are to be
*                     etau(imy) <= et(ij) < etau(imy + 1).
*            ik     - The order of the B-spline.
*            in     - The number of the orginal vertices.
*            et     - The new knot vector.
*            etau   - The old knot vector.
*            ep     - An array ep(ik) to local use. Such that we
*                     do not need to allocate the array locally after
*                     each call.
*
*
*
* OUTPUT   : jpl    - The negativ difference between the index in galfa
*                     and the real knot inserten matrix.
*            jfi    - The index of the first element in the line j in the
*                     the real knot inserten matrix whice is not zero.
*                     The element with the index (jfi+jpl) in galfa
*                     is the same as the element with index jfi in
*                     the real line j in the knot inserten matrix.
*            jla    - The index of the last element in the line j in the
*                     real knot inserten matrix whice is not zero.
*                     The element with the index (jla+jpl) in galfa
*                     is the same as the element with index jla in
*                     the real line j in the knot inserten matrix.
*            galfa  - A compressed line in the knot inserten matrix.
*            jstat  - status messages
*                                         > 0      : warning
*                                         = 0      : ok
*                                         < 0      : error
*
*
* METHOD     : Using the Oslo-algorithm
*
*
* REFERENCES : Making The Oslo algorithm more efficient.
*              by  T.Lyche and K.Moerken.
*              SIAM J.NUMER.ANAL  Vol. 23, No. 3, June 1986.
*
*-
* CALLS      :
*
* WRITTEN BY : Arne Laksaa, SI, 88-11.
* REVISED BY : Atgeirr F Rasmussen, Sintef, 06/04/2001. Go-ified.
*
**********************************************************************/
{
  int kp;                  /* Control variable in loop.     */
  int kv,kkv;              /* Help variables.               */
  double *ah;              /* Help pointer to galfa.        */
  double tbeta,tbeta1;     /* Help variables.               */
  double td1,td2;          /* Help variables.               */
  double *tu;              /* Pointer to the knot vector.   */

  std::vector<double> epvec(ik);
  double* ep = &epvec[0];
  
  /* Correction of imy to be sure that the the old knot etau(imy)
     is not passing the new knot et(ij). */
  
  kp=ij+1; kkv=ij+ik; in--;
  while ((et[kp] == etau[imy]) && kp<kkv) {kp++; imy--;}
  
  
  /* Counting the old and the new knot and copying the new knots
     in the area between et(ij) and et(ij+ik) in the array ep. */
  
  for (kp=imy+1,kv=0,ij++; ij<kkv; ij++)
    if (et[ij] == etau[kp]) kp++;
    else      ep[kv++] = et[ij];
  
  
  /* Compute the negativ difference between the index in galfa and
     the real knot inserten matrix. */
  
  *jpl=ik-imy-1;
  
  
  /* Changing the galfa so we may use the index in the real matrix. */
  
  galfa += *jpl;
  
  
  /* Initialise the last element. */
  
  galfa[imy] = 1;
  
  
  /* Here we go one time for each new knot from et(j+1)
     until et(j+k) we insert. */
  
  for (kp=0,kkv=ik-kv,ij=in+kv-1,in +=ik; kp<kv; kp++,kkv++,ep++)
    {
      /* The initialising:  The two first are not changing.
	 kkv = ik-kv , the nuber of old knots in the field.
	 This variabel is counting up to ik
	 (the order) during the loops.
	 ij = in+kv-1, minus the maximum of kp it
	 gives the index of the last
	 orginal vertices.
	 in = in+kv-1, the index of the last element in et. */
      
      
      /* Here we note the special case where we are at the
	 start of the matrix and we does not have a k-touple
	 knot at this end. */
      
      if (kp>=imy) tbeta1=(*ep - *etau)*(*galfa)/(etau[kkv] - *etau);
      else         tbeta1=0;
      
      *jfi=(1>imy-kp) ? 1 : imy-kp;
      *jla=(imy < ij-kp) ? imy : ij-kp;
      
      
      /* For details about this loop look in the reference. */
      
      for (et=etau+ *jfi,tu=etau+ *jla,ah=galfa+ *jfi; et<=tu; et++,ah++)
	{
	  td1 = *ep - *et;
	  td2 = et[kkv] - *ep;
	  tbeta = *ah/(td1+td2);
	  *(ah-1) = td2*tbeta + tbeta1;
	  tbeta1 = td1*tbeta;
	}
      
      
      /* Here we note the special case where we are at the
	 end of the matrix and we does not have a k-touple
	 knot at this end. */
      
      if (*jla<imy)
	{
	  et = etau + in;
	  *(ah-1) = tbeta1+(*et-*ep)*(*ah)/(*et - *(tu+1));
	} else  *(ah-1) = tbeta1;
    }
  
  
  /* Adjusting the index of first and last in galfa. */
  
  if (kv) (*jfi)--;
  else   *jfi = *jla = imy;
  
  if ((*jfi)<0) *jfi = 0;
  if ((*jla)>in-ik) *jla = in-ik;
  
}
