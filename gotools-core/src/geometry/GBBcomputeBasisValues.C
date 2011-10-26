#include "GoTools/geometry/BsplineBasis.h"
#include <algorithm>
#include <math.h>

using namespace Go;

//-----------------------------------------------------------------------------
std::vector<double>
BsplineBasis::computeBasisValues(double tval, int derivs ) const
//-----------------------------------------------------------------------------
{
    std::vector<double> resultvec(order_*(derivs+1));
    computeBasisValues(tval, &resultvec[0], derivs);
    return resultvec;
}

//-----------------------------------------------------------------------------
void BsplineBasis::computeBasisValues(const double tval, 
				      double* basisvals_start,
				      int derivs ,
				      double resolution) const
//-----------------------------------------------------------------------------
/*
*********************************************************************
*
*********************************************************************
*                                                                   
* PURPOSE    : To compute the value and ider first derivatives of the
*              ik (possibly) nonzero B-splines associated with the knot
*              vector et at the point ax.
*
*
*
* INPUT      : et     - Double array of dimension [in+ik] containing
*                       the knot vector.
*              ik     - The polynomial order of the B-splines associated
*                       with et.
*              in     - The dimension of the spline space associated with
*                       the knot vector et.
*              ax     - The point at which the B-spline values and derivatives
*                       are to be computed.
*              ider   - The number of derivatives to be computed.
*                       < 0 : Error.
*                       = 0 : Compute position.
*                       = 1 : Compute position and first derivative.
*                       etc.
*
*                
*
* INPUT/OUTPUT : ileft - Pointer to the interval in the knot vector
*                        where ax is located. The relation
*                          
*                          et[ileft] <= ax < et[ileft+1]
* 
*                        should hold. (If ax == et[in] then ileft should
*                        be in-1.)
*                        If ileft does not have the right value upon
*                        entry to the routine, its value will be changed
*                        to the value satisfying the above condition.
*
*
*
* OUTPUT     : ebder  - Double array of dimension [ik*(ider+1)] containing
*                       values of the ik nonzero B-splines and their
*                       derivatives at the point ax. These numbers are stored
*                       in the following order:
*                       First the 0-ider derivatives of the first nonzero
*                       B-splines at ax and so on. In other words, if ebder
*                       is considered as a two dimensional array its
*                       declaration in C would be ebder[ik,ider+1].
*              jstat  - Status messages  
*                                         > 0      : Warning.
*                                         = 0      : Ok.
*                                         < 0      : Error.
*
*
* METHOD     : It is well known that at a given point there are at most
*              ik nonzero B-splines. If ileft satisfies the condition
*              above, these B-splines, at the point ax, are
*
*                    B(ileft-ik+1),B(ileft-ik+2),...,B(ileft).
*
*              The number ileft is therefore determined at the beginning
*              of the routine.
*              To compute the values and derivatives of these B-splines,
*              the fundamental recurrence relations (k=ik, t=et, x=ax)
*
*                 B(i,k)(x)=(x-t(i))*Q(i,k-1)(x) + (t(i+k)-x)*Q(i+1,k-1)(x) (1)
*
*              and
*
*                DB(i,k)(x)=(k-1)*(Q(i,k-1)(x) - Q(i+1,k-1)(x))             (2)
*
*              are used. Here
*
*                             / B(i,k)(x)/(t(i+k)-t(i)),    if t(i) < t(i+k);
*                  Q(i,k)(x)=/
*                            \
*                             \           0,                otherwise;
*
*
*                            / 1,           if t(i) <= x < t(i+k);
*                 B(i,1)(x)=/
*                           \
*                            \ 0,           otherwise;
*
*              and DB(i,k)(x) denotes differentiation with respect to x.
*              (This makes the B-splines continuous from the right; note that
*              at the right hand end of the parameter value, at et[in],
*              ileft is still in-1 so that at that point the B-splines
*              are left continuous which is common practice.)
*              From (2) one obtains
*
*               D(r)B(i,k)(x)=(k-1)*(D(r-1)Q(i,k-1)(x) - D(r-1)Q(i+1,k-1)(x)),
*
*              (D(r) denotes r-fold differentation). From this it follows that
*              the r'th derivative of a B-spline is a linear combination
*              of B-splines of order k-r for r<k and zero for larger values
*              of r. Therefore, to compute the r'th derivative of all the
*              nonzero B-splines at x one starts with the one B-spline of
*              order 1 that is nonzero at x, namely B(ileft,1).
*              By using (1) the two B-splines of order 2 that are nonzero
*              at x are obtained, B(ileft-1,2) and B(ileft,2).
*              Continuing in this way one can compute the value of the k-r
*              B-splines of order k-r that are nonzero at x.
*              At this point one switches to (2) and obtains first the first
*              derivative of all the nonzero B-splines of order k-r+1 at x,
*              B(ileft-k+r,k-r+1),...,B(ileft,k-r+1), then the second
*              derivative of all the nonzero B-splines of order k-r+2 at x.
*              Applying (2) repeatedly in this fashion one eventually
*              obtains the r'th derivative of the k nonzero B-splines of order
*              k at x,
*
*                D(r)B(ileft-k+1,k),D(r)B(ileft-k+2,k),....,D(r)B(ileft,k).
*
*              Now the objective in this routine is not only to compute
*              the r'th derivative of all the nonzero B-splines for
*              one particular value of r, but for all r with 0 <= r <= ider.
*              With the above in mind, a natural way to do this is to start
*              with the one nonzero B-spline of order 1 and then apply (1)
*              k-ider-1 times to obtain all the nonzero B-splines of order
*              k-ider. During the next iteration both (1) and (2) are applied
*              to these B-spline values to obtain the values and first
*              derivatives of the nonzero B-splines of order k-ider+1.
*              In the next iteration (2) is applied to both the 0'th and
*              1'st derivatives to obtain the 1'st and 2'nd derivatives of
*              all nonzero derivatives of order k-ider+2, and (1) is applied
*              to the 0'th derivative to obtain the value of all nonzero
*              B-splines of order k-ider+2.
*              Continuing this process ider-2 more times one ends up with
*              the desired derivatives.
*
*              All these computations are performed in ebder and the author
*              finds it easiest to think of (the one dimensional array) ebder
*              as a two dimensional array of dimension [ider+1,ik]
*              (ider+1 rows and ik columns).
*          NB! If one conforms with C's view of two dimensional arrays, then
*              ebder would be of dimension [ik,ider+1] (first come the
*              ider+1 derivatives of the first nonzero B-spline and so on).
*
*              One then starts with the one nonzero B-spline of order 1 in
*              the lower right corner of ebder and then while applying (1)
*              fills in the last row of ebder from the right.
*              After k-ider-1 iterations the last k-ider elements of the last
*              row will be filled in with values of the k-ider nonzero
*              B-splines of order k-ider. These values are then copied up
*              to the second but last row and during the next iteration
*              (1) is applied to this row and (2) to the bottom row.
*              Then the second but last row is copied up to the third but
*              last row and (1) is applied to this row and (2) to the two
*              bottom rows. This is then repeated until ebder is filled
*              with the B-spline values in the first row, the first derivatives
*              in the second row and so on.
*
*              Note that (1) may be rewritten as
*
*             B(i,k)(x) = w(i)(x)*B(i,k-1)(x) + (1-w(i+1)(x))*B(i+1,k-1)(x) (3)
*
*              where
*
*                           / (x-t(i))/(t(i+k-1)-t(i)), if t(i) < t(i+k);
*                  w(i)(x)=/
*                          \
*                           \            0,             otherwise.
*
*              Note also that if et[ileft] <= x < et[ileft+1], then for
*              the nonzero B-splines B(ileft-k+1),...,B(ileft) the denominator
*              in the definition of w(i) is always nonzero. Since this
*              denominator is the same number as the one that divides
*              B(i,k-1) in (2), it is evident that with a correct knot vector,
*              no division by zero can occur.
*
* REFERENCES : Any work on elementary B-spline theory.
*
*-
* CALLS      : s1219 - Determine the value of ileft for the given ax.
*
* WRITTEN BY : Knut Moerken, University of Oslo, August 1988.
* REWISED BY : Vibeke Skytt, SI, 11.92.  Avoid memory error with order 1.
*
*********************************************************************
*/                                     
{
//  void s1220(et,ik,in,ileft,ax,ider,ebder,jstat)
//       double *et;
//       int    ik;
//       int    in;
//       int    *ileft;
//       double ax;
//       int    ider;
//       double ebder[];
//       int    *jstat;


  const double* et = &knots_[0];
  int ik = order_;
  int ider = derivs;
  double* ebder = basisvals_start;

  //  int kstat=0;        /* Local status variable.                          */
  int kdeg;           /* Convenience variable which is set to ik-1.      */
  int kleft;          /* Local version of ileft in order to avoid
			 the pointer.                                    */
  int ki,kj,ks;       /* Control variables in for loops and for stepping
			 through arrays.                                 */
  int ki1,ki2,kjh;    /* Control variables in for loops and for stepping
			 through arrays.                                 */
  int kder;           /* Local version of ider. All derivatives of order
			 higher than ik-1 are zero so
			 kder=min(ik-1,ider).                            */
  double td1,td2;     /* These variables are used to store the inverse of
			 the number that divides B(i,k-1) and B(i+1,k-1)
			 in (2) above.                                   */
  double tw1,tw2;     /* These variabels are used to store the factors
			 multiplying B(i,k-1) and B(i+1,k-1) in (1)
			 above.                                          */
  double ts1,ts2;     /* These variables are similar to td1 and td2 except
			 that they will also contain the appropriate
			 integer factor stemming from (2).               */
  ts2 = 0.0;
  double tt,tth;      /* Auxiliary variables used to avoid unnecessary
			 look ups in the knot vector.                    */
  /* Check the input. */                                       
  
  ALWAYS_ERROR_IF(ider < 0, "Number of derivatives must be >= 0.");
  
  // Find the right value of ileft and check input.
  // knotInterval may throw, in which case we have nothing delete
  // or release, so we let any exceptions propagate
  double val = tval;
  kleft = knotIntervalFuzzy(val, resolution);
  
  
  /* Initialize. */
  
  kdeg = ik - 1;
// #ifdef _MSC_VER
//   kder = (ider < ik-1) ? ider : ik-1;
// #else
  kder = std::min(ik-1,ider);
// #endif

  /* The fact that kder can be less than ider causes some problems.
     If kder < ider we know that ebder in the end should be a (kder+1)xik
     matrix augmented with ider-kder rows of zeros at the bottom.
     Since we store the matrix column by column, care must be taken to
     access the entries correctly.
     In the comments below ebder will usually be considered to be a
     (kder+1)xik matrix.
     ki2 is set to point to the last element of ebder (the lower right entry
     of the matrix, cf. above), and this element is set to one
     (this is the lower right corner of the (kder+1)*ik part of ebder). */
  
  ki2 = (ik-1)*(ider+1) + kder;
  ebder[ki2] = (double)1.0;
  
  if (ik == 1)
  {
     /* VSK. Constant. Task is done.  */
     
     return;
  }
  
  /* Get ready for the main iteration loop where (1) and/or (2) are applied
     each time. In these iterations ki1 will run through the entries of ebder
     to be computed, starting at the upper left corner and running down the
     rows.
     ki2 will follow ki1 but be one row ahead of ki1 to take care of the second
     term in (1) and (2).
     Note that when accessing the entries of ebder the convention is
     used that the pointer to the entry has to be incremented first.
     ki1 is therefore initialized to point to the element before the first
     element to be computed in the first iteration (the second but 
     last diagonal element of ebder).
     If (ider == ik-1) then we first have to copy up the last row and
     initialize ki1 to point to the entry above the second but last diagonal
     element of ebder.                                                       */
  
  ki1 = ki2 - ider - 2;
  if (kder == kdeg)
    {
      ebder[ki2-1] = (double)1.0;
      ki1 -= 1;
    }
  
  /* ki2 should be one row ahead of ki1 and each row has ider+1 entries. */
  
  ki2 = ki1 + ider + 1;
  
  /* Iterate and apply (1) and/or (2) each time. ks counts the degree of the
     B-splines whose values and derivaties are to be computed in this
     iteration.                                                          */
  
  for(ks=1; ks<ik; ks++)
    {
      
      /* In (1) and (2) the denominators that divide B(i,k-1) and B(i+1,k-1)
	 are on the form (t(i+k-1)-t(i)). Below kj is used as the (i) index
	 and kjh as the (i+k-1) index. It is the alternative form (3) of
	 (1) that is used below and tw2 is used as 1-w(i+1) and tw1 as
	 w(i).
	 For the first nonzero B-spline of degree ks, the first term in
	 (1) and (2) is zero.
	 kj is initialized to point to the first knot that gives a contribution
	 during this iteration and kjh to point to `t(kj+ks)=t(kleft+1)'.
	 If (t(kjh)-t(kj)) <= 0.0 there must be an error in the knot vector. */
      
      kj = kleft - ks + 1;
      kjh = kleft + 1;
      tt = et[kjh++];
      tth = tt - et[kj];

      // @@@ Temporary hack to avoid program crash. @jbt
//       if (tth <= (double)0.0)
// 	GO_ERROR("Error in knot vector.", CorruptData());
//       td2 = (double)1.0/tth;
      td2 = (double)1.0 / (tth + 1.0e-25); 

      tw2 = (tt-tval)*td2;
      
      ebder[++ki1] = tw2*ebder[++ki2];
      
      /* Check to see if there is either copying or differentiation to do. */
      
      if (ks >= kdeg-kder && kder > 0)
	{
	  /* Copy the first element of the row up to the previous
	     unless it is the last iteration.                     */
	  
	  if (ks < kdeg)
	    ebder[ki1-1] = ebder[ki1];
	  
	  /* Apply (2) to the rest of this column. Remember that this
	     is the first nonzero column of ebder so the first term in (2)
	     is zero.                                                     */
	  
	  ts2 = ks*td2;
	  for (ki=0; ki<ks-kdeg+kder; ki++)
	    ebder[++ki1] = -ts2*ebder[++ki2];
	  
	  /* Step to the top of the next column (the last time ebder
	     is full and there is no stepping to do (unless kder <ider). */
	  
	  ki1 += ider - kder + kdeg - ks;
	  ki2 = ki1 + ider + 1;
	}
      else
	{
	  
	  /* If there was no copying or differentiation to be done
	     we just step to the top of the next column. This step is zero
	     if there are no derivatives to be computed.                  */
	  
	  ki1 += ider;
	  ki2 += ider;
	}
      
      /*  Loop through the ks-1 middle columns of ebder. */
      
      for (kj=kleft-ks+2; kj<=kleft; kj++)
	{
	  
	  /* Compute the denominators and weights (w(i+1)) to be used.
	     See the comments above for more details.                  */
	  
	  tt = et[kjh++];
	  tth = tt - et[kj];
	  ALWAYS_ERROR_IF(tth <= (double)0.0, "Error in knot vector.");

	  td1 = td2; td2 = (double)1.0/tth;
	  tw1 = (double)1.0 - tw2; tw2 = (tt-tval)*td2;
	  
	  ki1 += 1;
	  ebder[ki1] = tw1*ebder[ki1] + tw2*ebder[++ki2];
	  
	  /* Check if there is copying and differentiation to be done. */
	  
	  if (ks >= kdeg-kder && kder > 0)
	    {
	      
	      /* Copy unless it is the last iteration. */
	      
	      if (ks < kdeg)
		ebder[ki1-1] = ebder[ki1];
	      
	      /* Do the differentiation. */
	      
	      ts1 = ts2; ts2 = ks*td2;
	      for (ki=0; ki<ks-kdeg+kder; ki++)
		{
		  ki1 += 1;
		  ebder[ki1] = ts1*ebder[ki1] - ts2*ebder[++ki2];
		}
	      
	      /* Jump to the next column. */
	      
	      ki1 += ider - kder + kdeg - ks;
	      ki2 = ki1 + ider + 1;
	    }
	  else
	    {
	      
	      /* Jump to the next column. */
	      
	      ki1 += ider;
	      ki2 += ider;
	    }
	}
      
      /* Compute the last column of ebder. Remember that now the last term
	 in (1) and (2) is zero, so there is no new td2 or tw2 to compute. */
      
      td1 = td2;
      tw1 = (double)1.0 - tw2;
      
      ki1 += 1;
      ebder[ki1] = tw1*ebder[ki1];
      
      /* Check if there is copying or differentiation to do. */
      
      if (ks >= kdeg-kder && kder > 0)
	{
	  
	  /* Copy. */
	  
	  if (ks < kdeg)
	    ebder[ki1-1] = ebder[ki1];
	  
	  /* Differentiate. */
	  
	  ts1 = ts2;
	  for (ki=0; ki<ks-kdeg+kder; ki++)
	    {
	      ki1 += 1;
	      ebder[ki1] = ts1*ebder[ki1];
	    }
	  
	  /* Move ki1 back to the first nonzero element of the last column.
	     Each column now has ks-kdeg+kder+2 nonzero elements.         */
	  
	  ki1 -= ks - kdeg + kder + 1;
	}
      
      /* Move ki1 from the first element of the last column to the element
	 prior to the first nonzero element in the first column.
	 ki2 is as usual one column ahead of ki1.                          */
      
      ki1 -= (ks+1)*(ider+1) + 1;
      ki2 = ki1 + ider + 1;
    }
  
  /* Set the remaining derivatives to zero. */
  
  for (ki=kder+1; ki<=ider; ki++)
    {
      ki1 = ki;
      for (kj=0; kj<ik; kj++)
	{
	  ebder[ki1] = (double)0.0;
	  ki1 += ider + 1;
	}
    }
  
  /* Successful computations.  */
  
  return;
}


//-----------------------------------------------------------------------------
void
BsplineBasis::computeBasisValues(const double* parvals_start,
				   const double* parvals_end,
				   double* basisvals_start,
				   int* knotinter_start,
				   int derivs) const
//-----------------------------------------------------------------------------
{
    for (; parvals_start < parvals_end; ++parvals_start) {
	computeBasisValues(*parvals_start, basisvals_start, derivs);
	*knotinter_start = lastKnotInterval();
	++knotinter_start;
	basisvals_start += order()*(derivs+1);
    }
}




//-----------------------------------------------------------------------------
std::vector<double>
BsplineBasis::computeBasisValuesLeft(double tval, int derivs ) const
//-----------------------------------------------------------------------------
{
    std::vector<double> resultvec(order_*(derivs+1));
    computeBasisValuesLeft(tval, &resultvec[0], derivs);
    return resultvec;
}

//-----------------------------------------------------------------------------
void
BsplineBasis::computeBasisValuesLeft(double tval, 
				     double*      basisvals_start,
				     int          derivs,
				     double       resolution) const
//-----------------------------------------------------------------------------
{
    // Method taken from s1227. If tval is a knot, make new basis ending in tval.

    // We locate the interval in which tval belongs.
    int left = knotIntervalFuzzy(tval, resolution);

    // Adjust knot interval for numerical noice
    if (left < num_coefs_-1 && knots_[left+1]-tval <= resolution)
	left++;
    while (left < num_coefs_-1 && knots_[left] == knots_[left+1])
	left++;

    // If tval is not a knot, left evaluation is exactly the same as right eval.
    if (fabs(tval-startparam()) <= resolution ||  
	fabs(knots_[left]-tval) > resolution) {
      computeBasisValues(tval, basisvals_start, derivs);
      return;
    }

    /* To force the derivative to be taken from the left we artificially
       shorten the curve if ax==st[kleft]  */

    int mult = knotMultiplicity(tval);
    --last_knot_interval_;

    // Copy the knots in the basis.
    int new_num_coefs = left - mult + 1;
    std::vector<double> new_knots(new_num_coefs + order_);
    std::copy(knots_.begin(), knots_.begin() + new_knots.size(), new_knots.begin());

    BsplineBasis new_basis(new_num_coefs, order_, new_knots.begin());
    new_basis.computeBasisValues(tval, basisvals_start, derivs);
    last_knot_interval_ = left - mult;
    if (last_knot_interval_ < order_-1)
	last_knot_interval_ = order_ - 1;
}

//-----------------------------------------------------------------------------
void
BsplineBasis::computeBasisValuesLeft(const double* parvals_start,
				       const double* parvals_end,
				       double* basisvals_start,
				       int* knotinter_start,
				       int derivs) const
//-----------------------------------------------------------------------------
{
    for (; parvals_start < parvals_end; ++parvals_start) {
	computeBasisValuesLeft(*parvals_start, basisvals_start, derivs);
	*knotinter_start = lastKnotInterval();
	++knotinter_start;
	basisvals_start += order()*(derivs+1);
    }
}

