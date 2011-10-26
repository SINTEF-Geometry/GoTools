//===========================================================================
//
// File : volume_ratder.C
//
// Created: Fri Nov 14 12:15:00 2008
//
// Author: Kjell Fredrik Pettersen
//
// Revision: $Id: volume_ratder.C,v 1.2 2009-01-13 12:12:00 kfp Exp $
//
// Description:
//
//===========================================================================

#include "GoTools/geometry/SplineUtils.h"
#include <vector>

namespace Go
{

//===========================================================================
void volume_ratder(double const eder[],int idim,int ider,double gder[])
//===========================================================================
{
  double w0;           /* The reciprocal of the denominator.     */

  int* binom;       /* Array for storing binomial coeffs.     */
  std::vector<int> binomvec;
  int bidum[10];    /* Fixed space for binom */
  int *binomu=0;    /* Pointer to binomial coefficients in u. */
  int *binomv=0;    /* Pointer to binomial coefficients in v. */
  int *binomw=0;    /* Pointer to binomial coefficients in w. */

  double *sum1=0;   /* Leibnitz expansion in u                */
  double *sum2=0;   /* Leibnitz expansion in u and v.         */
  double *sum3=0;   /* Leibnitz expansion in u, v and w.      */
  std::vector<double> sum1v;
  std::vector<double> sum2v;
  std::vector<double> sum3v;
  double sumdum1[4];   /* Fixed space for sum1.                  */
  double sumdum2[4];   /* Fixed space for sum2.                  */
  double sumdum3[4];   /* Fixed space for sum2.                  */

  int idu;             /* Count through derivatives in u.        */
  int idv;             /* Count through derivatives in v.        */
  int idw;             /* Count through derivatives in w.        */

  int tot,temp1;       /* Temporary variables.                   */

  ALWAYS_ERROR_IF(ider<0, "Less than zero derivatives ?!?");
  ALWAYS_ERROR_IF(idim<1, "Less than zero derivatives ?!?");

  /* Find 1/denominator. */ 

  if (fabs(eder[idim])<1e-13)
    w0 = (double)1.0; // Maybe we should throw instead?
  else
    w0 = 1.0 / eder[idim];

  /* If we're only asked for position, we'll do it
     now and exit for the sake of speed. */

  if (ider == 0)
    {
      for (int i = 0; i < idim; ++i)
	gder[i] = eder[i] * w0;

      return;
    }

  /* Set up  binomial coefficients.
     Use new array only when ider > 3. */

  if (ider > 3)
    { 
      binomvec.resize(((ider+1)*(ider+2))>>1);
      binom = &binomvec[0];
    }
  else
    binom = bidum;

  for( int j = 0, k = 0; j <= ider; ++j, k += j)
    {
      /* Calculate the new row of binomial coefficients. */
  
      binom[k] = binom[k+j] = 1;
      
      for(int i = k+1; i < k+j; ++i)
	binom[i] = binom[i-j-1] + binom[i-j];
    }

  /* Set up space for Leibnitz expansions.
     Use new arrays only when idim > 4. */

  if (idim > 4)
  { 
    sum1v.resize(idim);
    sum1 = &sum1v[0];
    sum2v.resize(idim);
    sum2 = &sum2v[0];
    sum3v.resize(idim);
    sum3 = &sum3v[0];
  }
  else
  { 
    sum1=sumdum1;
    sum2=sumdum2;
    sum3=sumdum3;
  }

  /* Loop through derivatives in u, v and w. */

  for(idw=0,binomw=binom; idw<=ider; idw++,binomw+=idw)
    for(idv=0,binomv=binom; idv<=ider-idw; idv++,binomv+=idv)
      for(idu=0,binomu=binom; idu<=ider-idv-idw; idu++,binomu+=idu)
	{
	  if(idu+idv+idw == 0)
	    {
	      /* Position is a special case. */
	      for(int i = 0; i < idim; ++i)
		gder[i] = eder[i] * w0;
	      continue;
	    }

          /* Calculating each coefficient of the (idu,idv,idw)'th
	     derivative of the rational surface (in gder).
        
  	     This requires calculating the Liebnitz sum from
  	     the subarray of gder and eder */

          /* Calculate the Leibnitz sum. */

          for(int i = 0; i < idim; ++i)
            sum3[i] = 0.0;        

          for(int iw = 0; iw <= idw; ++iw)
	    {
	      for(int i = 0; i < idim; ++i)
		sum2[i] = 0.0;	               
	      int iwtemp = idw-iw;
 
	      for(int iv = 0; iv <= idv; ++iv)
		{
		  for(int i = 0; i < idim; ++i)
		    sum1[i] = 0.0;	               

		  int ivtemp = idv-iv;
		  for(int iu = 0; iu <= idu; ++iu)
		    {
		      int iutemp = idu-iu;

		      tot = iv + iw;
		      temp1 = ((tot * (tot+1)) >> 1) + iw;
		      tot += iu;
		      temp1 += (tot * (tot+1) * (tot+2)) / 6;

		      int igder = temp1 * idim;

		      tot = ivtemp + iwtemp;
		      temp1 = ((tot * (tot+1)) >> 1) + iwtemp;
		      tot += iutemp;
		      temp1 += (tot * (tot+1) * (tot+2)) / 6;

		      int ieder   = temp1 * (idim+1) + idim;

		      /* Add the next Leibnitz term unless we
			 have reached the last one (the unknown). */
  
		      if(iu<idu || iv<idv || iw<idw)
			{
			  /* If iu=0 or iu=idu, the u binomial
			     coefficient is 1 so don't multiply. */
  
			  if(iu>0 && iu<idu)
			    {
			      double temp = (double)binomu[iu] * eder[ieder];
			      for(int i = 0; i < idim; ++i, ++igder)
				sum1[i] += temp * gder[igder];
			    }
			  else
			    for(int i = 0; i < idim; ++i, ++igder)
			      sum1[i] += eder[ieder] * gder[igder];
			}
		    }
  
		  /* If iv=0 or iv=idv, the v binomial
		     coefficient is 1 so don't multiply. */
  
		  if(iv > 0 && iv < idv)
		    for(int i = 0; i < idim; ++i)
		      sum2[i] += (double)binomv[iv] * sum1[i];
		  else
		    for(int i = 0; i < idim; ++i)
		      sum2[i] += sum1[i];		    
		}


	      /* If iw=0 or iw=idw, the w binomial
		 coefficient is 1 so don't multiply. */
  
	      if(iw > 0 && iw < idw)
		for(int i = 0; i < idim; ++i)
  	          sum3[i] += (double)binomw[iw] * sum2[i];
	      else
		for(int i = 0; i < idim; ++i)
		  sum3[i] += sum2[i];		    
	    }

          /* Calculate indices in eder and gder. */
    
          tot = idv + idw;
          temp1 = ((tot * (tot+1)) >> 1) + idw;
	  tot += idu;
	  temp1 += (tot * (tot+1) * (tot+2)) / 6;

          int pos_gder = temp1 * idim;
          int pos_eder = pos_gder + temp1;

          for(int i = 0; i < idim; ++i, ++pos_gder, ++pos_eder)
            gder[pos_gder] = (eder[pos_eder] - sum3[i]) * w0;
	}

  return;
}

}
