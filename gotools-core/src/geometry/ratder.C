//===========================================================================
//                                                                           
// File: ratder.C                                                            
//                                                                           
// Created: Fri Apr  6 13:45:59 2001                                         
//                                                                           
// Author: Atgeirr F Rasmussen <atgeirr@sintef.no>
//                                                                           
// Revision: $Id: ratder.C,v 1.7 2005-06-24 07:50:27 oan Exp $
//                                                                           
// Description:
//                                                                           
//===========================================================================


#include "GoTools/geometry/SplineUtils.h"
#include <vector>

namespace Go
{

void curve_ratder(double const eder[],int idim,int ider,double gder[])
/*
*********************************************************************
*
* PURPOSE    : To calculate the ider derivatives of a rational

*              point described in homogenous coordinates
*
* INPUT      : eder    - The derivatives in homogenous coordinates
*                        In sequence:
*                         Position (x,y,...h)
*                         1st der (x,y,...h)
*                         2nd der (x,y,...h)
*                         etc.
*              idim    - The dimension of the non homogenous space
*              ider    - The number of input derivatives
*
*
* OUTPUT     : gder    - The derivatives in the nonhomogenous space
*
*
* METHOD     :  The curve P(u) can be written as the quotient
*               P(u) = T(u) / w(u) where T and w are ordinary splines.
*               The dimensions of T and w are idim and 1
*               respectively. The array eder contains position
*               and derivatives of the idim+1 dimensional curve
*               (T(u),w(u)).
*
*               Now, since wP = T, we find, by the Leibnitz formula,
*
*                 k
*                         k!     (k-i) (i)         (k)
*                sum   -------- w     P       =   T    .
*                      i!(k-i)!
*                i=0
*
*               Therefore
*
*
*                   --         k-1                      --
*             (k)   |   (k)             k!     (k-i) (i) |
*            P    = |  T    -  sum   -------- w     P    | / w .
*                   |                i!(k-i)!            |
*                   --         i=0                      --
*
*               This formula is applied recursively to evaluate P's derivatives.
*
*                                                          MF.
*
* CALLS      :
*
* WRITTEN BY : Tor Dokken, SI, Oslo, Norway. 1988-des-1988
* REVISED BY : Michael Floater, SI, 30/9/91 Removed division by t=1.
* REWRITTEN BY : Michael Floater, SI, 16/12/91. New algorithm.
* REWRITTEN BY : Michael Floater, SI, 25/8/92. Extend to arbitrary
*                   number of derivatives (by Leibnitz). Finally!
* REVISED BY : Paal Fugelli, SINTEF, 07/07-94. Added free'ing of binom and
*              initiation to SISL_NULL to avoid memory leakage.
* REVISED BY : Atgeirr F Rasmussen, SINTEF, 06/04/2001. Go and C++-ified.
* 
*********************************************************************
*/
{
  double w0;           /* The denominator.                       */
  int ki;              /* Count through dimensions.              */
  int id;              /* Count through derivatives.             */
  std::vector<double> binom;
  double sum;          /* Binomial (Leibnitz) expansion.         */
  int idimp1;          /* idim + 1.                              */
  int iw;              /* Pointer to a weight.                   */
  int igder;           /* Pointer to already calculated derivs.  */
  int i,j,k;           /* Counters.                              */
  int iwfix;           /* Initial value of iw in Leibnitz loop.  */

  ALWAYS_ERROR_IF(ider<0, "Less than zero derivatives ?!?");
  ALWAYS_ERROR_IF(idim<1, "Less than zero derivatives ?!?");

  idimp1 = idim + 1;

  /* Find denominator. */

  w0 = eder[idim];
  if (fabs(w0)<1e-13) w0 = (double)1.0; // Maybe we should throw instead?

  /* Set up initial binomial coefficient (1). */

  binom.resize(ider+1);

  binom[0] = 1;

  /* Calculate position first. */

  for(ki=0; ki<idim; ki++)
  {
      gder[ki] = eder[ki] / w0;
  }



  /* Then derivatives if there are any. */

  for(id=1,j=idim,k=idimp1; id<=ider; id++,k++)
  {
      /* Calculate the new row of binomial coefficients. */

      binom[id] = 1;

      for(i=id-1; i>=1; i--)
      {
	  binom[i] += binom[i-1];
      }


      /* Run through the idim dimensions, calculating each
	 coefficient of the id'th derivative of
	 the rational curve (in gder). */

      iwfix = k + idim;

      for(ki=0; ki<idim; ki++,j++,k++)
      {
	  /* Calculate the Leibnitz sum (the binomial
	     coefficient in the first term is always 1). */

	  sum = eder[iwfix] * gder[ki];

          for(i=1,igder=idim+ki,iw=iwfix-idimp1;
	  i<id;
	  i++,igder+=idim,iw-=idimp1)
          {
	      sum += (double)binom[i] * eder[iw] * gder[igder];
          }

	  gder[j] = (eder[k] - sum) / w0;

      }

  }

  /* Done. */

}


void surface_ratder(double const eder[],int idim,int ider,double gder[])
/*
*********************************************************************
*                                                                   
* PURPOSE    : To calculate the value and ider*ider derivatives of 
*              a rational B-spline surface.
*
* INPUT      : eder    - Double array of dimenson [(ider+1)*(ider+2)*(idim+1)/2]
*                        containing the position and the derivative vectors
*                        of the homogeneous surface at the point with parameter value
*                        (epar[0],epar[1]).
*                        (idim+1 is the number of components of each B-spline
*                        coefficient, i.e. the dimension of the homogemous
*                        space in which the surface lies.)
*                        These vectors are stored in the following order:
*                        First the idim+1 components of the position vector,
*                        then the idim+1 components of the D(1,0) vector,
*                        then the idim+1 components of the D(0,1) vector,
*                        then the idim+1 components of the D(2,0) vector,
*                        followed by D(1,1), D(0,2)
*                        and so on up to the idim+1 components of the D(0,ider).
*              idim    - The dimension of the non homogenous space
*              ider    - The number of input derivatives with respect
*                        to both parameter directions.
*
*
* OUTPUT     : jstat   - Status message
*                                        >0      : Warning
*                                        =0      : ok
*                                        <0      : Error
*              gder    - Double array of dimension [(ider+1)*(ider+2)*idim/2]
*                        containing the position and the derivative vectors
*                        of the surface at the point with parameter value
*                        (epar[0],epar[1]).
*                        (idim is the number of components of each B-spline
*                        coefficient, i.e. the dimension of the Euclidean
*                        space in which the surface lies.)
*                        These vectors are stored in the following order:
*                        First the idim components of the position vector,
*                        then the idim components of the D(1,0) vector,
*                        then the idim components of the D(0,1) vector,
*                        then the idim components of the D(2,0) vector,
*                        followed by D(1,1), D(0,2)
*                        and so on up to the idim components of the D(0,ider).
*
*
* METHOD     :  The surface P(u,v) can be written as the quotient
*               P(u,v) = T(u,v) / w(u,v) where T and w are ordinary splines.
*               The dimensions of T and w are idim and 1
*               respectively. The array eder contains position
*               and derivatives of the idim+1 dimensional surface
*               (T(u,v),w(u,v)).
*
*               Now, since wP = T, we find, by the Leibnitz formula,
*
*      k   l
*                  k!       l!     (k-i,l-j) (i,j)         (k,l)
*     sum sum   -------- -------- w         P         =   T       .
*               i!(k-i)! j!(l-j)!
*     i=0 j=0
*
*               Therefore
*               
*
*              --            k   l                                     --
*      (k,l)   |   (k,l)                k!       l!     (k-i,l-j) (i,j) |    
*     P      = |  T      -  sum sum  -------- -------- w         P      | / w .
*              |                     i!(k-i)! j!(l-j)!                  |
*              --           i=0 j=0                                    --
*                               i+j<k+l
*
*               This formula is applied recursively to evaluate P's derivatives.
*
*                                                          MF.
*
*
*
* CALLS      :
*
* WRITTEN BY : Michael Floater, SI, 3.9.92.
*                Essentially the same as s6sratder
*                except that we work with triangular matrices
*                ((0,0), (1,0), (0,1), (2,0), (1,1), ...)
*                instead of rectangular ones.
* Revised by : Christophe Rene Birkeland, SINTEF Oslo, May 1993.
*              Error message corrected
* Revised by : Atgeirr F Rasmussen, SINTEF, 06/04/2001. Go and C++-ified.
*
*********************************************************************
*/
{
  double w0;           /* The denominator.                       */
  int ki;              /* Count through dimensions.              */
  int idu;             /* Count through derivatives in u.        */
  int idv;             /* Count through derivatives in v.        */
  int* binom;
  std::vector<int> binomvec;
  int *binomu=0;    /* Pointer to binomial coefficients in u. */
  int *binomv=0;    /* Pointer to binomial coefficients in v. */
  double *sum1=0;   /* Leibnitz expansion in u                */
  double *sum2=0;   /* Leibnitz expansion in u and v.         */
  std::vector<double> sum1v;
  std::vector<double> sum2v;
  double sumdum1[4];   /* Fixed space for sum1.                  */
  double sumdum2[4];   /* Fixed space for sum2.                  */
  int idimp1;          /* idim + 1.                              */
  int iw;              /* Pointer to a weight.                   */
  int igder;           /* Pointer to already calculated derivs.  */
  int i,iu,iv,j,k;     /* Counters.                              */
  int iderp1;          /* ider + 1.                              */
  int igrow;           /* (ider+1) * idim.                       */  
  int iwrow;           /* (ider+1) * idimp1.                     */  
  int iutemp,ivtemp;   /* Used to find next weight in the sum.   */
  int tot,temp1;       /* Temporary variables.                   */
  int bidum[10];       /* Array for storing binomial coeffs.     */
  double temp;         /* Temporary multiple.                    */
  
  ALWAYS_ERROR_IF(ider<0, "Less than zero derivatives ?!?");
  ALWAYS_ERROR_IF(idim<1, "Less than zero derivatives ?!?");

  /* Find denominator. */ 
  
  w0 = eder[idim];
  if (fabs(w0)<1e-13) w0 = (double)1.0; // Maybe we should throw instead?

  /* If we're only asked for position, we'll do it
     now and exit for the sake of speed. */

  if(ider == 0)
  {
    for(ki=0; ki<idim; ki++)
      gder[ki] = eder[ki] / w0;

    return;
  }

  /* Set up some constants. */

  idimp1  = idim + 1;
  iderp1 = ider + 1;
  igrow   = iderp1 * idim;
  iwrow   = igrow + iderp1;  /* = iderp1 * idimp1 */

  /* Set up  binomial coefficients.
     Use new array only when ider > 3. */

  if (ider > 3)
  { 
    binomvec.resize((iderp1*(iderp1+1))/2);
    binom = &binomvec[0];
  }
  else
  { 
    binom = bidum;
  }

  for(j=0,k=0; j<=ider; j++,k+=j)
  {
      /* Calculate the new row of binomial coefficients. */
  
      binom[k] = 1;
  
      for(i=k+1; i<k+j; i++)
      {
          binom[i] = binom[i-j-1] + binom[i-j];
      }

      binom[k+j] = 1;
  }

  /* Set up space for sum1 and sum2 if necessary.
     Use new arrays only when idim > 4. */

  if (idim > 4)
  { 
    sum1v.resize(idim);
    sum1 = &sum1v[0];
    sum2v.resize(idim);
    sum2 = &sum2v[0];
  }
  else
  { 
    sum1=sumdum1;
    sum2=sumdum2;
  }
  
  /* Loop through derivatives in u and v. */

  for(idv=0,binomv=binom; idv<=ider; idv++,binomv+=idv)
  {
    for(idu=0,binomu=binom; idu<=ider-idv; idu++,binomu+=idu)
    {
      if(idu == 0 && idv == 0)
      {
          /* Position is a special case. */
    
          for(ki=0; ki<idim; ki++)
            gder[ki] = eder[ki] / w0;
      }
      else
      {
          /* Calculate indices in eder and gder. */
    
          tot = idu + idv;
          temp1 = ((tot * (tot+1)) >> 1) + idv;
    
          j = temp1 * idim;
          k = j + temp1;

          /* Calculating each coefficient of the (idu,idv)'th
	     derivative of the rational surface (in gder).
        
  	     This requires calculating the Liebnitz sum from
  	     the subarray of gder (0,..,idu, 0,...,idv) and
             the subarray of eder (0,..,idu, 0,...,idv). */

          /* Calculate the Leibnitz sum. */

          for(ki=0; ki<idim; ki++)
            sum2[ki] = (double)0.0;        

          for(iv=0; iv<=idv; iv++)
          {
            for(ki=0; ki<idim; ki++)
               sum1[ki] = (double)0.0;	               
            ivtemp = idv-iv;

            for(iu=0; iu<=idu; iu++)
            {
                tot = iu + iv;
                temp1 = ((tot * (tot+1)) >> 1) + iv;

	        igder = temp1 * idim;
                iutemp = idu-iu;

                tot = iutemp + ivtemp;
	        temp1 = ((tot * (tot+1)) >> 1) + ivtemp;

                iw   = temp1 * idimp1 + idim;

     	      /* Add the next Leibnitz term unless we
       		 have reached the last one (the unknown). */
  
                if(iu<idu || iv<idv)
  	        {
  	       	  /* If iu=0 or iu=idu, the u binomial
  	       	     coefficient is 1 so don't multiply. */
  
  	            if(iu>0 && iu<idu)
  	       	    {
  		      temp = (double)binomu[iu] * eder[iw];
                      for(ki=0; ki<idim; ki++,igder++)
  	       	         sum1[ki] += temp * gder[igder];
  		     }
  		     else
                       for(ki=0; ki<idim; ki++,igder++)
  		         sum1[ki] += eder[iw] * gder[igder];
                }
            }
  
  	    /* If iv=0 or iv=idv, the v binomial
  	       coefficient is 1 so don't multiply. */
  
  	    if(iv>0 && iv<idv)
              for(ki=0; ki<idim; ki++)
  	          sum2[ki] += (double)binomv[iv] * sum1[ki];
  	    else
              for(ki=0; ki<idim; ki++)
		 sum2[ki] += sum1[ki];		    
          }
          for(ki=0; ki<idim; ki++,j++,k++)
            gder[j] = (eder[k] - sum2[ki]) / w0;
      }
    }
  }  

  return;
}


} // namespace Go

