/*
 * Copyright (C) 1998, 2000-2007, 2010, 2011, 2012, 2013 SINTEF ICT,
 * Applied Mathematics, Norway.
 *
 * Contact information: E-mail: tor.dokken@sintef.no                      
 * SINTEF ICT, Department of Applied Mathematics,                         
 * P.O. Box 124 Blindern,                                                 
 * 0314 Oslo, Norway.                                                     
 *
 * This file is part of GoTools.
 *
 * GoTools is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Affero General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version. 
 *
 * GoTools is distributed in the hope that it will be useful,        
 * but WITHOUT ANY WARRANTY; without even the implied warranty of         
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the          
 * GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public
 * License along with GoTools. If not, see
 * <http://www.gnu.org/licenses/>.
 *
 * In accordance with Section 7(b) of the GNU Affero General Public
 * License, a covered work must retain the producer line in every data
 * file that is created or manipulated using GoTools.
 *
 * Other Usage
 * You can be released from the requirements of the license by purchasing
 * a commercial license. Buying such a license is mandatory as soon as you
 * develop commercial activities involving the GoTools library without
 * disclosing the source code of your own applications.
 *
 * This file may be used in accordance with the terms contained in a
 * written agreement between you and SINTEF ICT. 
 */

/* #include "constant.h" */
/* #include "boxdir.dcl" */
/* #include "surf.dcl" */
#include "GoTools/creators/CreatorsOffsetUtils.h"
#include <vector>

using std::vector;

namespace Go
{

void OffsetUtils::blend_s1421(const SplineSurface* ps, double aoffset, int ider,
                              const Point& epar, int& ilfs, int& ilft,
                              std::vector<Point>& eoffpnt, 
                              std::vector<Point>& epnt, int* jstat)
    /*
*********************************************************************
*
*********************************************************************
*                                                                   
* PURPOSE    : Evaluate the surface and the offset surface in the 
*              distance aoffset at the parameter pair value epar.
*
*
*
* INPUT      : ps     - The input B-spline surface.
*              aoffset- The offset distance.
*                       If idim=2 a positive signe on this value put the
*                       offset on the side of the positive normal vector,
*                       and a negative sign puts the offset on the sign
*                       of the negative normal vector.
*                       If idim=3 the offset is determined by the cross
*                       product of the tangent vector and the anorm vector.
*                       The offset distance is multiplied by this vector.
*              epar   - Parameter pair value at which to compute position
*                       and derivatives.
*              ider   - The number of derivatives to compute.
*                       < 0: Error
*                       = 0: Compute position
*                       = 1: Compute position and first derivative
*                       etc.
* INPUT/OUTPUT:
*              ilfs   - Pointer to the interval in the knot vector in 1.
*                       parameter direction.
*              ilft   - Pointer to the interval in the knot vector in 2.
*                       parameter direction.
*
* OUTPUT     : 
*              jstat  - status messages  
*                                         = 2      : Surface is degenerate
*                                                    at the point, normal
*                                                    has zero length
*                                         = 1      : Surface is close to
*                                                    degenerate at the point
*                                                    Angle between tangents,
*                                                    less than angular tolerance
*                                         = 0      : ok
*                                         < 0      : error
*              eoffpnt -  Array of dimension idim*(1+2+...+(ider+1)).
*                       The sequence is position, first derivative in first
*                       parameter direction, first derivative in second parameter
*                       direction, (2,0) derivative, (1,1) derivative, (0,2) 
*                       derivative etc. The offset surface is evaluated. Higher 
*                       order derivatives than second order is only allowed
*                       when aoffset = 0. And finally we include normal in point,
*                       given that ider >= 1. The normal is not normalized.
*              epnt   - Array of dimension idim*(1+2+...+(ider+1)). 
*                       The sequence is position, first derivative in first
*                       parameter direction, first derivative in second parameter
*                       direction, (2,0) derivative, (1,1) derivative, (0,2) 
*                       derivative etc. The input surface is evaluated. Higher 
*                       order derivatives than second order is only allowed
*                       when aoffset = 0. And finally we include normal in point,
*                       given that ider >= 1. The normal is not normalized.
*
* METHOD     : Points, derivatives, cross derivatives and normals are calcu-
*              lated at the selected point on the surface. The selected point
*              corresponds to the input parameter-pair (epar).
*              The offset values are calculated in the following way:
*
*              O(u,v)   = P(u,v)   + aoffset (N(u,v)/n(u,v))
*
*              O(u,v)   = P(u,v)   + aoffset (N(u,v)/n(u,v))
*                    u          u                           u
*
*              O(u,v)   = P(u,v)   + aoffset (N(u,v)/n(u,v))
*                    v          v                           v
*
*              O(u,v)   = P(u,v)   + aoffset (N(u,v)/n(u,v))
*                    uu         uu                          uu
*
*              O(u,v)   = P(u,v)   + aoffset (N(u,v)/n(u,v))
*                    uv         uv                          uv
*
*              O(u,v)   = P(u,v)   + aoffset (N(u,v)/n(u,v))
*                    vv         vv                          vv
*
*
* EXAMPLE OF USE:
*
*
* REFERENCES :
*
*
*-
* CALLS      : s1421     - Evaluate the surface at the parameter pair value
*                          (epar).
*              s6length, s6crss, s6scpr, s6ang, s6err
*
* WRITTEN BY : Per Evensen,  SI, 89-4.
* REWISED BY : Vibeke Skytt, SI, 89-10.
*
*********************************************************************
*/
{

    double ANGULAR_TOLERANCE = 0.001; // @@sbr

    int kdim;          /* Dimension of the space in which the surface lies.*/
    int kder = ider+1; /* Number of derivatives to calculate.              */
    int knmb;          /* Number of array-elements.                        */
    int knmb2;         /* Number of elements in output arrays.             */
    /* double snorm[3];   /\* Pointer to array containing the normal.          *\/ */
    Point snorm(3);
    // double sdum[30];   /* Array used instead of allocation.                */
    // double *sder=NULL; /* Pointer to array containing the derivatives.     */
    vector<Point> sder;


    double tnorm,tang=(double)0.0;

    double tl3,tl5;    /* tl3=tnorm**3. tl5=tnorm**5.                */
    double tder1,tder2;/* Length of tangents.                        */
    Point snoru1(3);  /* The first derivative of the surface normal
			   in first parameter direction, 1. part.     */
    Point snoru2(3);  /* The first derivative of the surface normal
			   in first parameter direction, 2. part.     */
    Point snoru(3);   /* The first derivative of the surface normal
			   in first parameter direction.              */
    Point snorv1(3);  /* The first derivative of the surface normal
			   in second parameter direction, 1. part.    */
    Point snorv2(3);  /* The first derivative of the surface normal
			   in second parameter direction, 2. part.    */
    Point snorv(3);   /* The first derivative of the surface normal
			   in second parameter direction.             */
    Point snoruu1(3); /* The second derivative of the surface normal
			   in first parameter direction. 1. part.     */
    Point snoruu2(3); /* The second derivative of the surface normal
			   in first parameter direction. 2. part.     */
    Point snoruu3(3); /* The second derivative of the surface normal
			   in first parameter direction. 3. part.     */
    Point snoruu(3);  /* The second derivative of the surface normal
			   in first parameter direction.              */
    Point snoruv1(3); /* The cross derivative of the surface normal,
			   1. part.                                   */
    Point snoruv2(3); /* The cross derivative of the surface normal,
			   2. part.                                   */
    Point snoruv3(3); /* The cross derivative of the surface normal,
			   3. part.                                   */
    Point snoruv(3);  /* The cross derivative of the surface normal.*/
    Point snorvv1(3); /* The second derivative of the surface normal
			   in second parameter direction. 1. part.    */
    Point snorvv2(3); /* The second derivative of the surface normal
			   in second parameter direction. 2. part.    */
    Point snorvv3(3); /* The second derivative of the surface normal
			   in second parameter direction. 3. part.    */
    Point snorvv(3);  /* The second derivative of the surface normal
			   in second parameter direction.             */
    double tnun;       /* The scalar product; snorm snoru.           */
    double tnvn;       /* The scalar product; snorm snorv.           */
    double tunvn;      /* The scalar product; snoru snorv.           */
    double tnuun;      /* The scalar product: snorm snoruu.          */
    double tnuvn;      /* The scalar product; snorm snoruv.          */
    double tnvvn;      /* The scalar product: snorm snorvv.          */
    double tunun;      /* The scalar product: snoru snoru.           */
    double tvnvn;      /* The scalar product: snorv snorv.           */
    Point snorus(3);  /* The first derivative of the surface normal
			   divided by the length of the normal
			   in first parameter direction.              */
    Point snorvs(3);  /* The first derivative of the surface normal
			   divided by the length of the normal
			   in second parameter direction.             */
    Point snoruus(3); /* The second derivative of the surface normal
			   in the first parameter direction divided
			   by the length of the normal.               */
    Point snoruvs(3); /* The cross derivative of the surface normal
			   divided by the length of the normal.       */
    Point snorvvs(3); /* The second derivative of the surface normal
			   in second parameter direction divided by
			   the length of the normal.                  */





    kdim = ps->dimension();
    /* if (kdim !=3) goto err105; */
    ALWAYS_ERROR_IF(kdim != 3,
		"Dimension of input surface must be 3!");

    /* if (ider > 2 && !DEQUAL(aoffset,DNULL)) goto err105; */
    ALWAYS_ERROR_IF((ider > 2) && (aoffset != 0.0),
		"We do not handle derivatives larger than 2 if offset dist > 0.0.");


    /* Only allocate space if sdum is too small. */

    int ki;
    for (knmb=0,ki=1; ki <= kder+1; knmb+=ki,++ki);
    if (knmb>10)
	//     sder = newarray(knmb,DOUBLE);
	sder.resize(knmb);
    else
	//     sder = &sdum[0];
	sder.resize(10);

    /* Get number of elements in output array.  */

    for (knmb2=0,ki=1; ki <= kder; knmb2+=ki,++ki);

    /* if (DEQUAL(aoffset,DNULL)) */
    if (aoffset == 0.0) {
	/* Evaluate the surface (ps) at the parameter value spar.  */
	ps->point(sder, epar[0], epar[1], kder);
	copy(sder.begin(), sder.begin() + knmb2, epnt.begin());
	ilfs = ps->basis_u().knotInterval(epar[0]);
	ilft = ps->basis_v().knotInterval(epar[1]);
	ps->normal(epnt[knmb2], epar[0], epar[1]);

	/* The offset surface is equal to the input surface in this case.
	   Return the result of the evaluation of the offset surface.     */
	eoffpnt = epnt;
   
    } else {



	/* Evaluate the surface (ps) at the parameter value spar.  */
	ps->point(sder, epar[0], epar[1], kder);
	copy(sder.begin(), sder.begin() + knmb2, epnt.begin());
	ilfs = ps->basis_u().knotInterval(epar[0]);
	ilft = ps->basis_v().knotInterval(epar[1]);
	// We next must compute normal in the same point.
	ps->normal(epnt[knmb2], epar[0], epar[1]);

	/*  Copy to output array.  */
	//     memcopy(eder,sder,knmb2,DOUBLE);
	//     memcopy(enorm,snorm,kdim,DOUBLE);
    
	/*  Calculate angle between tangents */
	tang = sder[1].angle(sder[2]);
	if (tang == 0.0)
	    *jstat = 2;
	else if (tang <= ANGULAR_TOLERANCE)
	    *jstat = 1;   
	else
	    *jstat = 0;
    
	/* Calculate normal of surface and length of normal.   */
	snorm = sder[1]%sder[2];
	tnorm = snorm.length();
    
	/*
	 *   The offset point, O(s,t,o), where o is the offset distance, is
	 *   found by the expression:
	 *
	 *   O = P + o*(N/n)    
	 */
    
	/* Calculate position of offset point. */
	eoffpnt[0] = sder[0] + aoffset*snorm/tnorm;
    
	if (ider > 0) {
	    /* Calculate length of tangents. */
	    tder1 = sder[1].length();
	    tder2 = sder[2].length();	

	    /* The tangent length might be very different from 1. Scale it and
	       higher order derivatives to give tangent length one. */

	    /* --------------------------------------------------------------------------------*/
	    /* Scaling to avoid over(under) flow                                               */
	    /* --------------------------------------------------------------------------------*/

	    sder[1] /= tder1;
	    sder[2] /= tder2;
	    sder[3] /= (tder1*tder1);
	    sder[4] /= (tder1*tder2);
	    sder[5] /= (tder2*tder2);
	    sder[6] /= (tder1*tder1*tder1);
	    sder[7] /= (tder1*tder1*tder2);
	    sder[8] /= (tder1*tder2*tder2);
	    sder[9] /= (tder2*tder2*tder2);

	    /* Calculate normal of offset surface in point and length of normal.  */

	    snorm = sder[1]%sder[2];
	    tnorm = snorm.length();
	    tl3 = tnorm*tnorm*tnorm;
	    tl5 = tl3*tnorm*tnorm;

	    /* --------------------------------------------------------------------------------*/
	    /* Evaluation of the derivatives of the (scaled) normal vector.                    */
	    /* --------------------------------------------------------------------------------*/


	    /*
	     *   The first derivative of the surface normal in first parameter direction
	     *   is calculated by the expression:
	     *
	     *
	     *   N  = P  x P   + P x P
	     *    u    uu   v     u   uv
	     */
	
	    /* Calculate first derivative of the surface normal in first
	       parameter direction. */

	    snoru1 = sder[3]%sder[2];
	    snoru2 = sder[1]%sder[4];
	
	    for (ki = 0;ki < kdim;++ki)
		snoru[ki]=snoru1[ki]+snoru2[ki];
	
	    /*
	     *   The first derivative of the surface normal in second parameter direction
	     *   is calculated by the expression:
	     *
	     *
	     *   N  = P  x P   + P x P
	     *    v    uv   v     u   vv
	     */
	
	    /* Calculate first derivative of the surface normal in second
	       parameter direction. */
	
	    snorv1 = sder[4]%sder[2];
	    snorv2 = sder[1]%sder[5];
	
	    snorv = snorv1 + snorv2;
	
    

	    /* -------------------------------------------------------------------------------*/
	    /* Evaluation of the derivatives of the derivatives of the unit normal.           */
	    /* -------------------------------------------------------------------------------*/


	    /*
	     *   The first derivative of the surface normal divided by the length of
	     *   the normal in first parameter direction is calculated by the 
	     *   expression:
	     *
	     *                             3
	     *   (N/n)  = N /n - (N * N )/n  N
	     *        u    u           u
	     */
	
	    /* Calculate first derivative of the surface normal divided
	       by the length of the normal in first parameter direction. */

	    tnun = snorm*snoru;

	    snorus = snoru/tnorm - tnun*snorm/tl3;

	    /*
	     *   The first derivative of the surface normal divided by the length of
	     *   the normal in second parameter direction is calculated by the 
	     *   expression:
	     *
	     *                             3
	     *   (N/n)  = N /n - (N * N )/n  N
	     *        v    v           v
	     */
	
	    /* Calculate first derivative of the surface normal divided
	       by the length of the normal in second parameter direction. */
	
	    tnvn = snorm*snorv;

	    snorvs = snorv/tnorm - tnvn*snorm/tl3;


	    /* -------------------------------------------------------------------------------*/
	    /* Evaluation of the derivatives of the offset surface.                           */
	    /* -------------------------------------------------------------------------------*/

	    /*
	     *   The first derivative in first parameter direction offset point, 
	     *   where o is the offset distance, is found by the expression:
	     *
	     *           
	     *   O  = P  + o*(N/n)    
	     *    u    u          u
	     */
	
	    /* Calculate position of first derivative in first parameter direction
	       offset point. */
	    eoffpnt[1] = sder[1] + aoffset*snorus;

	    /*
	     *   The first derivative in second parameter direction offset point, 
	     *   where o is the offset distance, is found by the expression:
	     *
	     *
	     *   O  = P  + o*(N/n)    
	     *    v    v          v
	     */

	    /* Calculate position of first derivative in second parameter direction
	       offset point. */
	    eoffpnt[2] = sder[2] + aoffset*snorvs;

	    /* Calculate normal of offset surface.   */
	    eoffpnt[knmb2] = eoffpnt[1]%eoffpnt[2];
	
	    if (ider > 1) {

		/* ---------------------------------------------------------------------------*/
		/* Evaluation of the derivatives of the (scaled) normal vector.               */
		/* ---------------------------------------------------------------------------*/

		/*  The second derivative of the surface normal in the first parameter
		 *   direction is calculated by the expression:
		 *
		 *
		 *   N   = P   x P   + 2*(P  x P )  + P  x P
		 *    uu    uuu   v       uu   uv     u    uuv
		 */
	
		/* Calculate second derivative of surface normal in first direction. */
	
		// 	    s6crss(sder+6*kdim,sder+2*kdim,snoruu1);
		snoruu1 = sder[6]%sder[2];
		// 	    s6crss(sder+3*kdim,sder+4*kdim,snoruu2);
		snoruu2 = sder[3]%sder[4];
		// 	    s6crss(sder+kdim,sder+7*kdim,snoruu3);
		snoruu3 = sder[1]%sder[7];
		snoruu = snoruu1 + 2.0*snoruu2 + snoruu3;
		
		/*
		 *   The cross derivative of the surface normal is calculated by the 
		 *   expression:
		 *
		 *
		 *   N   = P   x P   + P  x P   + P  x P
		 *    uv    uuv   v     uu  vv     u    uvv
		 */
	    
		/* Calculate cross derivative of the surface normal. */
		snoruv1 = sder[7]%sder[2];
		snoruv2 = sder[3]%sder[5];
		snoruv3 = sder[1]%sder[8];

		snoruv = snoruv1 + snoruv2 + snoruv3;
	
	
		/*  The second derivative of the surface normal in the second parameter
		 *   direction is calculated by the expression:
		 *
		 *
		 *   N   = P   x P  + 2*(P  x P  ) + P x P
		 *    vv    uvv   v       uv   vv     u   vvv
		 */
	
		/* Calculate second derivative of normal in second direction.  */
		snorvv1 = sder[7]%sder[2];
		snorvv2 = sder[4]%sder[5];
		snorvv3 = sder[1]%sder[9];
	
		snorvv = snorvv1 + 2.0* snorvv2 + snorvv3;
	
    		/*  The second derivative of the surface normal in the first parameter
		 *   direction divided by the length of the normal is calculated by the
		 *   expression :
		 *
		 *                                 3               3              3
		 *  (N/n)  = N  /n - 2N  (N * N )/n - N (N * N  )/n - N (N * N )/n +
		 *      uu    uu        u      u              uu          u   u
		 *
		 *                      2   5
		 *           3N (N * N ) / n
		 *                    u
		 */
	
		/* Calculate second derivative of surface normal in first parameter
		   direction divided by the length of the normal.                   */
		tunun = snoru*snoru;
		tnuun = snorm*snoruu;
	
		snoruus = snoruu/tnorm - 2.0*snoru*tnun/tl3 -
		    snorm*((tnuun + tunun)/tl3 - 3*tnun*tnun/tl5);

	
		/*
		 *   The cross derivative of the surface normal divided by the length of
		 *   the normal is calculated by expression:
		 *
		 *                               3                3
		 *   (N/n)   = N  /n - (N * N )/n  N  - (N * N )/n  N  -  ==>
		 *        uv    uv           v      u         u      v
		 *
		 *                        3                3                        5
		 *             (N  * N )/n  N - (N * N  )/n  N + 3(N * N )(N * N )/n  N
		 *               u    v               uv                u       v
		 */
	
		/* Calculate cross derivative of the surface normal divided
		   by the length of the normal. */
		tunvn = snoru*snorv;
		tnuvn = snorm*snoruv;
		snoruvs = snoruv/tnorm - tnvn*snoru/tl3 - tnun*snorv/tl3 -
		    tunvn*snorm/tl3 - tnuvn*snorm/tl3 + 3.0*tnun*tnvn*snorm/tl5;
	    
		/*   The second derivative of the surface normal in the second parameter
		 *   direction divided by the length of the normal is calculated by the
		 *   expression :
		 *
		 *                                 3               3              3
		 *  (N/n)  = N  /n - 2N  (N * N )/n - N (N * N  )/n - N (N * N )/n +
		 *      vv    vv        v      v              vv          v   v
		 *
		 *                      2   5
		 *           3N (N * N ) / n
		 *                    v
		 */
	
		/* Calculate second derivative of surface normal in second parameter
		   direction divided by the length of the normal.                   */
		tvnvn = snorv*snorv;
		tnvvn = snorm*snorvv;
		
		snorvvs = snorvv/tnorm - 2.0*snorv*tnvn/tl3 -
		    snorm*((tnvvn + tvnvn)/tl3 - 3.0*tnvn*tnvn/tl5);

		/* ------------------------------------------------------------------------------*/
		/* Evaluation of the derivatives of the offset surface.                          */
		/* ------------------------------------------------------------------------------*/
    
		/*  The second derivative in first parameter direction offset point,
		 *  where o is the offset distance is found by the expression:
		 *
		 *
		 *  O  = P   + o*(N/n)
		 *   uu   uu          uu
		 */
	    
		/* Calculate position of second derivative in first parameter direction
		   offset point.  */
		eoffpnt[3] = sder[3] + aoffset*snoruus;
	
		/*
		 *   The cross derivative offset point, 
		 *   where o is the offset distance, is found by the expression:
		 *
		 *
		 *   O   = P   + o*(N/n)    
		 *    uv    uv          uv
		 */
	    
		/* Calculate position of cross derivative offset point. */
		eoffpnt[4] = sder[4] + aoffset*snoruvs;
	
	
		/*  The second derivative in second parameter direction offset point,
		 *  where o is the offset distance is found by the expression:
		 *
		 *
		 *  O  = P   + o*(N/n)
		 *   vv   vv          vv
		 */
	
		/* Calculate position of second derivative in second parameter direction
		   offset point.  */
		eoffpnt[5] = sder[5] + aoffset*snorvvs;
	    }
	
	    /* --------------------------------------------------------------------------------*/
	    /* Scale the derivatives to match the original parametrization.                    */
	    /* --------------------------------------------------------------------------------*/

	    eoffpnt[1] *= tder1;

	    eoffpnt[2] *= tder2;

	    eoffpnt[3] *= (tder1*tder1);
	     
	    eoffpnt[4] *= (tder1*tder2);

	    eoffpnt[5] *= (tder2*tder2);

	}
    }

    /* Point and derivatives calculated, with a possible offset distance. */

    *jstat = 0;
}

} // end namespace Go
