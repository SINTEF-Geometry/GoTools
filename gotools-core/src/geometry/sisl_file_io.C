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

#include "GoTools/geometry/sisl_file_io.h"
#include "sislP.h"

#include <vector>
#include <memory>

#ifndef FALSE
#define FALSE 0
#endif

#ifdef __BORLANDC__
using std::FILE;
using std::sscanf;
#endif

void read_non_comment(FILE* fp, char* string)
{
  char c;    /* char to find end of line and end of file */

  fscanf(fp,"%s",string);

  /* read until first string is a letter */
  while (string[0] == '$')
    {
      /* find end of line */
	while ((c = (char)getc(fp)) != '\n')
	{
	  if (c == EOF)
	    {
	      string[0] = c;
	      string[1] = '\0';
	      return;
	    }
	}
      /* get data */
      fscanf(fp,"%s",string);
    }

}


void file_to_obj(FILE* fp, SISLObject** wo, int* jstat)

/*
*********************************************************************************
*
* Purpose : To read a B-spline communication file.
*           File format hplabs.
*           Naming conventions SI/HP.
*
* Input   : from_file   - File name.
*
* Output  : wo          - Pointer to a SISLObject.
*           jstat       - status flag
*                       >= 0 - OK.
*                       <  0 - Error
*
* Written by : Bjoern Olav Hoset, SI, 06-89.
*
********************************************************************************
*/
{
  SISLCurve *qc=NULL;        /* Pointer to struct SISLCurve */
  SISLSurf *qs=NULL;         /* Pointer to struct SISLSurf  */
  SISLObject *qo=NULL;       /* Pointer to struct SISLObject */
//   FILE *fp=NULL;/* FILE pointer to read from */
  char string[80];       /* character to read $ which identifies a comment */
  int value;             /* type of B-Spline and numbers */
  int value2;            /* second type of B-Spline and numbers */
  int kcopy=1;           /* Flag */
  int kcuopen;           /* Open/closed flag for curves */
  int kcuopen1;          /* Open/closed flag for surf  */
  int kcuopen2;          /* Open/closed flag for surf  */
  int kdim;              /* dimension of target space */
  int kdim1;             /* dimension + 1             */
  int kind;              /* identify point, curve, surface */
  int kk1;               /* Order of curve or in 1st parameter direction of surface */
  int kk2;               /* Order in 2nd parameter direction of surface */
  int kn1;               /* # of vertices of curve or in 1st parameter direction of surf */
  int kn2=0;               /* # of vertices in 2nd parameter direction of surface */
  double *st1=NULL;      /* Knot vector of curve or in 1st parameter direction of surface */
  double *st2=NULL;      /* Knot vector in 2nd parameter direction of surface */
  double *scoef=NULL;    /* Vertices of curve or surface */
  double dummy;          /* trashbin for fourth value in vertices */
  double help[4];        /* Help array in computing transformated vertic. */
  double matrix[16];     /* instance matrix */
  double *sptr=NULL;     /* Dummy pointer */
  int ki, kj;            /* Counters */
  int kkind;
  int *ind_hex=0;
  int hex=FALSE;
  int *ind_2 = (int*)&dummy;

  //vector__double__ seg_u;
  //vector__double__ seg_v;
  //vector__int__ prio_u;
  //vector__int__ prio_v;
  //vector__double__ pbox;

  /*
   * Open file with name from_file
   * -----------------------------
   */

  /*
   * Read data from file.
   * --------------------
   */

  read_non_comment(fp,string);

  /*
   * Test if end of file is reached
   * ------------------------------
   */

  if (string[0] == EOF) goto err001;
  else
  {
     sscanf(string,"%d",&value);

     if (value >=10)
     {
	hex=true;
	value -=10;
     }
     switch (value)
     {
	case 0: /* object is a curve or a surface */

	  read_non_comment(fp,string);
	  sscanf(string,"%d",&value);
	  for (ki = 0; (string[ki] = (char)getc(fp)) != '\n'; ki++);
	  string[ki] = '\0';

	  /*
	   * if a second numerical value is specified in this line,
	   * a surface is assumed
	   * ------------------------------------------------------
	   */

	  kind = sscanf(string,"%d",&value2);

	  /* now back to reading the object */

	  if (kind !=1 )
	    {
	      /*
	       * object is a curve
	       * -----------------
	       */

	      /*
	       * order
	       * ------
	       */

	      kk1 = value;

	      /*
	       * number of vertices
	       * ------------------
	       */

	      read_non_comment(fp,string);
	      sscanf(string,"%d",&kn1);

	      /*
	       * dimension
	       * ---------
	       */

	      read_non_comment(fp,string);
	      sscanf(string,"%d",&kdim);

	      /*
	       * curve open/closed
	       * ------------------
	       */

	      read_non_comment(fp,string);
	      sscanf(string,"%d",&kcuopen);

	      /*
	       * rational B-Spline or not
	       * ------------------------
	       */

	      read_non_comment(fp,string);
	      sscanf(string,"%d",&value);
	      /* if (value != 0) goto err002; */

	      if (value != 0)
	      {
		 kdim1 = kdim+1;
		 kkind = 2;
	      }
	      else
	      {
		 kdim1 = kdim;
		 kkind = 1;
	      }

	      /*
	       * Allocate space for knot vector
	       * ------------------------------
	       */

	      st1 = new double[kk1+kn1];
	      if (!st1) goto err101;

	      /*
	       * Read in knot vector
	       * -------------------
	       */

	      read_non_comment(fp,string);
	      if (hex)
	      {
		 ind_hex = (int *)st1;
		 sscanf(string,"%x",ind_hex);
		 fscanf(fp,"%x",ind_hex+1);
		 for (ind_hex=(int *)(st1+1),ki=1;ki < kk1+kn1;ki++,ind_hex+=2)
		    fscanf(fp,"%x %x",ind_hex,ind_hex+1);
	      }
	      else
	      {
		 sscanf(string,"%lf",st1);
		 for (sptr=st1+1,ki=1;ki < kk1+kn1;ki++,sptr++)
		    fscanf(fp,"%lf",sptr);
	      }

	      /*
	       * Allocate space for control vertices
	       * -----------------------------------
	       */

	      scoef = new double[kn1*kdim1];
	      if (!scoef) goto err101;

	      /*
	       * Read in control vertices
	       * ------------------------
	       */

	      read_non_comment(fp,string);
	      if (hex)
	      {
		 ind_hex = (int *)scoef;
		 sscanf(string,"%x",ind_hex);
		 fscanf(fp,"%x",ind_hex+1);
		 if (kdim > 1)
		 {
		    ind_hex +=2;
		    fscanf(fp,"%x %x",ind_hex,ind_hex+1);
		    if (kdim > 2)
		    {
		       ind_hex +=2;
		       fscanf(fp,"%x %x",ind_hex,ind_hex+1);
		    }
		 }
		 fscanf(fp,"%x %x",ind_2,ind_2+1);
		 if (kkind == 2)
		 {
		    ind_hex +=2;
		    scoef[kdim] = dummy;
		 }
		 for (ki = kdim1; ki < kn1*kn2*kdim1; ki += kdim1)
		 {
		    ind_hex +=2;
		    fscanf(fp,"%x %x",ind_hex,ind_hex+1);
		    if (kdim > 1)
		    {
		       ind_hex +=2;
		       fscanf(fp,"%x %x",ind_hex,ind_hex+1);
		       if (kdim > 2)
		       {
			  ind_hex +=2;
			  fscanf(fp,"%x %x",ind_hex,ind_hex+1);
		       }
		    }
		    fscanf(fp,"%x %x",ind_2,ind_2+1);
		    if (kkind == 2)
		    {
		       ind_hex +=2;
		       scoef[ki+kdim] = dummy;
		    }
		 }

	      }
	      else
	      {
		 sscanf(string,"%lf",scoef);
		 if (kdim > 1) fscanf(fp,"%lf",&scoef[1]);
		 if (kdim > 2) fscanf(fp,"%lf",&scoef[2]);
		 fscanf(fp,"%lf",&dummy);
		 if (kkind == 2) scoef[kdim] = dummy;
		 for (ki = kdim1; ki < kn1*kdim1; ki += kdim1)
		 {
		    fscanf(fp,"%lf",&scoef[ki]);
		    if (kdim > 1) fscanf(fp,"%lf",&scoef[ki+1]);
		    if (kdim > 2) fscanf(fp,"%lf",&scoef[ki+2]);
		    fscanf(fp,"%lf",&dummy);
		    if (kkind == 2) scoef[ki+kdim] = dummy;
		 }
	      }
	      /*
	       * transformation data
	       * -------------------
	       */

	      read_non_comment(fp,string);
	      sscanf(string,"%lf",matrix);
	      for (ki = 1; ki < (kdim+1)*(kdim+1); ki ++)
		  fscanf(fp,"%lf",&matrix[ki]);

	      /*
	       * transform vertices with instance matrix
	       * ---------------------------------------
	       */


	      /*
	       * test whether matrix is unit matrix
	       * ----------------------------------
	       */

	      if (kdim      == 1           &&
		  matrix[0] == (double)1.0 &&
		  matrix[1] == (double)0.0 &&
		  matrix[2] == (double)0.0 &&
		  matrix[3] == (double)1.0 )
	      {
	      }
	      else if (kdim      == 2           &&
		       matrix[0] == (double)1.0 &&
		       matrix[1] == (double)0.0 &&
		       matrix[2] == (double)0.0 &&
		       matrix[3] == (double)0.0 &&
		       matrix[4] == (double)1.0 &&
		       matrix[5] == (double)0.0 &&
		       matrix[6] == (double)0.0 &&
		       matrix[7] == (double)0.0 &&
		       matrix[8] == (double)1.0)
	      {
	      }
	      else if (kdim       == 3           &&
		       matrix[0]  == (double)1.0 &&
		       matrix[1]  == (double)0.0 &&
		       matrix[2]  == (double)0.0 &&
		       matrix[3]  == (double)0.0 &&
		       matrix[4]  == (double)0.0 &&
		       matrix[5]  == (double)1.0 &&
		       matrix[6]  == (double)0.0 &&
		       matrix[7]  == (double)0.0 &&
		       matrix[8]  == (double)0.0 &&
		       matrix[9]  == (double)0.0 &&
		       matrix[10] == (double)1.0 &&
		       matrix[11] == (double)0.0 &&
		       matrix[12] == (double)0.0 &&
		       matrix[13] == (double)0.0 &&
		       matrix[14] == (double)0.0 &&
		       matrix[15] == (double)1.0)
	      {
	      }

	      else
	      {
		 if (kdim > 3) goto err005;
		 for (ki = 0; ki < kn1*kdim; ki += kdim)
		 {
		    kdim1 = kdim + 1;
		    for (kj = 0; kj < kdim1; kj++)
		       help[kj]  = s6scpr(&scoef[ki],
					   &matrix[kj*kdim1],kdim)
			            + matrix[kj*kdim1 + kdim];

		    if (help[kdim] == 0.0) goto err004;

		    for (kj = 0; kj < kdim; kj++)
		       scoef[ki + kj] = help[kj]/help[kdim];
		 }
	      }



	      /*
	      * Create an instance of SISLCurve
	      * ---------------------------
	      */

	      qc = newCurve(kn1,kk1,st1,scoef,kkind,kdim,kcopy);
	      if (!qc) goto err101;
	      else
	      {
		  /*
		   * Fill in open/closed flag
		   * ------------------------
		   */

		  qc->cuopen = kcuopen;

		  /*
		   * Create an instace of SISLObject and connect curve to the SISLObject.
		   * ------------------------------------------------------------
		   */

		  qo = newObject(SISLCURVE);
		  if (!qo) goto err101;
		  else
		    {
		      qo->c1 = qc;
		      *wo    = qo;
		    }
		}

	      /*
	       * Remember to free allocated space.
	       * ---------------------------------
	       */

	      if (st1) delete [] st1;
	      if (scoef) delete [] scoef;

	    }
	  else
	    /* object is a surface */
	    {
	      /*
	       * orders of parameter directions
	       * ------------------------------
	       */

	      kk1 = value;
	      kk2 = value2;

	      /*
	       * number of vertices
	       * ------------------
	       */

	      read_non_comment(fp,string);
	      sscanf(string,"%d",&kn1);
	      fscanf(fp,"%d",&kn2);


	      /*
	       * dimension
	       * ---------
	       */

	      read_non_comment(fp,string);
	      sscanf(string,"%d",&kdim);
	      kdim1 = kdim;

	      /*
	       * surface open/closed
	       * ------------------
	       */

	      read_non_comment(fp,string);
	      sscanf(string,"%d",&kcuopen1);
	      fscanf(fp,"%d",&kcuopen2);

	      /*
	       * rational B-Spline or not
	       * ------------------------
	       */

	      read_non_comment(fp,string);
	      sscanf(string,"%d",&value);
	      /* if (value != 0) goto err002;

	      fscanf(fp,"%d",&value); */
	      if (value) kdim1++;

	      /*
	       * Allocate space for first knot vector
	       * ------------------------------------
	       */

	      st1 = new double[kk1+kn1];
	      if (!st1) goto err101;

	      /*
	       * Read in first knot vector.
	       * --------------------------
	       */

	      read_non_comment(fp,string);
	      if (hex)
	      {
		 ind_hex = (int *)st1;
		 sscanf(string,"%x",ind_hex);
		 fscanf(fp,"%x",ind_hex+1);
		 for (ind_hex=(int*)(st1+1),ki=1;ki < kk1+kn1;ki++,ind_hex+=2)
		    fscanf(fp,"%x %x",ind_hex,ind_hex+1);
	      }
	      else
	      {

		 sscanf(string,"%lf",st1);
		 for (sptr = st1+1,ki = 1; ki < kn1+kk1; sptr++,ki++)
		    fscanf(fp,"%lf",sptr);
	      }
	      /*
	       * Allocate space for second knot vector
	       * ------------------------------------
	       */

	      st2 = new double[kk2+kn2];
	      if (!st2) goto err101;

	      /*
	       * Read in second knot vector.
	       * --------------------------
	       */

	      read_non_comment(fp,string);
	      if (hex)
	      {
		 ind_hex = (int *)st2;
		 sscanf(string,"%x",ind_hex);
		 fscanf(fp,"%x",ind_hex+1);
		 for (ind_hex=(int*)(st2+1),ki=1;ki < kk2+kn2;ki++,ind_hex+=2)
		    fscanf(fp,"%x %x",ind_hex,ind_hex+1);
	      }
	      else
	      {

		 sscanf(string,"%lf",st2);
		 for (sptr = st2+1,ki = 1; ki < kn2+kk2; sptr++,ki++)
		    fscanf(fp,"%lf",sptr);
	      }

	      /*
	       * Allocate space for control vertices
	       * -----------------------------------
	       */

	      scoef = new double[kn1*kn2*kdim1];
	      if (!scoef) goto err101;

	      /*
	       * Read in control vertices.
	       * -------------------------
	       */

	      read_non_comment(fp,string);
	      if (hex)
	      {
		 ind_hex = (int *)scoef;
		 sscanf(string,"%x",ind_hex);
		 fscanf(fp,"%x",ind_hex+1);
		 if (kdim > 1)
		 {
		    ind_hex +=2;
		    fscanf(fp,"%x %x",ind_hex,ind_hex+1);
		    if (kdim > 2)
		    {
		       ind_hex +=2;
		       fscanf(fp,"%x %x",ind_hex,ind_hex+1);
		    }
		 }
		 fscanf(fp,"%x %x",ind_2,ind_2+1);
		 if (value)
		 {
		    ind_hex +=2;
		    scoef[kdim] = dummy;
		 }
		 for (ki = kdim1; ki < kn1*kn2*kdim1; ki += kdim1)
		 {
		    ind_hex +=2;
		    fscanf(fp,"%x %x",ind_hex,ind_hex+1);
		    if (kdim > 1)
		    {
		       ind_hex +=2;
		       fscanf(fp,"%x %x",ind_hex,ind_hex+1);
		       if (kdim > 2)
		       {
			  ind_hex +=2;
			  fscanf(fp,"%x %x",ind_hex,ind_hex+1);
		       }
		    }
		    fscanf(fp,"%x %x",ind_2,ind_2+1);
		    if (value)
		    {
		       ind_hex +=2;
		       scoef[ki+kdim] = dummy;
		    }
		 }

	      }
	      else
	      {
		 sscanf(string,"%lf",scoef);
		 if (kdim > 1) fscanf(fp,"%lf",&scoef[1]);
		 if (kdim > 2) fscanf(fp,"%lf",&scoef[2]);
		 fscanf(fp,"%lf",&dummy);
		 if (value) scoef[kdim] = dummy;
		 for (ki = kdim1; ki < kn1*kn2*kdim1; ki += kdim1)
		 {
		    fscanf(fp,"%lf",&scoef[ki]);
		    if (kdim > 1) fscanf(fp,"%lf",&scoef[ki+1]);
		    if (kdim > 2) fscanf(fp,"%lf",&scoef[ki+2]);
		    fscanf(fp,"%lf",&dummy);
		    if (value) scoef[ki+kdim] = dummy;
		 }

	      }
	      /*
	       * transformation data
	       * -------------------
	       */

	      read_non_comment(fp,string);
	      sscanf(string,"%lf",matrix);
	      for (ki = 1; ki < (kdim+1)*(kdim+1); ki ++)
		  fscanf(fp,"%lf",&matrix[ki]);

	      /*
	       * transform vertices with instance matrix
	       * ---------------------------------------
	       */


	      /*
	       * test whether matrix is unit matrix
	       * ----------------------------------
	       */

	      if (kdim      == 1           &&
		  matrix[0] == (double)1.0 &&
		  matrix[1] == (double)0.0 &&
		  matrix[2] == (double)0.0 &&
		  matrix[3] == (double)1.0 )
	      {
	      }
	      else if (kdim      == 2           &&
		       matrix[0] == (double)1.0 &&
		       matrix[1] == (double)0.0 &&
		       matrix[2] == (double)0.0 &&
		       matrix[3] == (double)0.0 &&
		       matrix[4] == (double)1.0 &&
		       matrix[5] == (double)0.0 &&
		       matrix[6] == (double)0.0 &&
		       matrix[7] == (double)0.0 &&
		       matrix[8] == (double)1.0)
	      {
	      }
	      else if (kdim       == 3           &&
		       matrix[0]  == (double)1.0 &&
		       matrix[1]  == (double)0.0 &&
		       matrix[2]  == (double)0.0 &&
		       matrix[3]  == (double)0.0 &&
		       matrix[4]  == (double)0.0 &&
		       matrix[5]  == (double)1.0 &&
		       matrix[6]  == (double)0.0 &&
		       matrix[7]  == (double)0.0 &&
		       matrix[8]  == (double)0.0 &&
		       matrix[9]  == (double)0.0 &&
		       matrix[10] == (double)1.0 &&
		       matrix[11] == (double)0.0 &&
		       matrix[12] == (double)0.0 &&
		       matrix[13] == (double)0.0 &&
		       matrix[14] == (double)0.0 &&
		       matrix[15] == (double)1.0)
	      {
	      }
              else if (1)
              {
              }

	      else
		{
		  if (kdim > 3 || kdim < 1) goto err005;
		  if (!value)
		    {
		      for (ki = 0; ki < kn2*kn1*kdim; ki += kdim)
			{
			  kdim1 = kdim + 1;
			  for (kj = 0; kj < kdim1; kj++)
			    help[kj]  = s6scpr(&scoef[ki],
					       &matrix[kj*kdim1],kdim)
			      + matrix[kj*kdim1 + kdim];
			  
			  if (help[kdim] == 0.0) goto err004;
			  
			  for (kj = 0; kj < kdim; kj++)
			    scoef[ki + kj] = help[kj]/help[kdim];
			}
		    }
		  else 
		    goto err002;
		}
	      
	      // Read segmentation and parameter box for surface
	      // Nope I do not
	      //stat = readSegments(fp,seg_u,seg_v,prio_u,prio_v,pbox);
	      //if (stat < 0)
	      //goto error;

	      /*
	       * Create an instance of SISLSurf.
	       * ---------------------------
	       */

	      qs = newSurf(kn1,kn2,kk1,kk2,st1,st2,scoef,(value)?2:1,kdim,kcopy);
	      if (!qs) goto err101;
	      else
		{
		  // Set segmentation values and parameter box
                  // Nope this is ignored here
		  //if (pbox.size() == 4)
		  //SI_setparbox(qs, pbox.begin());

		  //if (seg_u.size() > 0 || seg_v.size() > 0)
		  //{
		  //stat = SI_set_seg(qs, seg_u.size(), seg_u.begin(),
		  //	      prio_u.begin(),seg_v.size(),
		  //	      seg_v.begin(), prio_v.begin());
		  //if (stat < 0)
		  //  goto error;

		  //		  }
		  /*
		   * Create an instace of Object and connect surface to the Object.
		   * --------------------------------------------------------------
		   */

		  qs->cuopen_1 = kcuopen1;
		  qs->cuopen_2 = kcuopen2;

		  qo = newObject(SISLSURFACE);
		  if (!qo) goto err101;
		  else
		    {
		      qo->s1 = qs;
		      *wo    = qo;
		    }
		}

	      /*
	       * Remember to free allocated space.
	       * ---------------------------------
	       */

	      if (st1) delete [] st1;
	      if (st2) delete [] st2;
	      if (scoef) delete[] scoef;


	    } /* end of if curve or surface */
	  break;

	default:
	  goto err003;
	}
    }/* end else */

  /*
   * Object read from file.
   * ----------------------
   */

  *jstat = 0;
  goto out;

  /*
   * Error in opening/reading from file.
   * -----------------------------------
   */

  err001: *jstat = -1;
          s6err("file_to_obj",*jstat,0);
          goto out;

  /*
   * Error, object in file was of type rational.
   * -------------------------------------------
   */

 err002: *jstat = -2;
          s6err("file_to_obj",*jstat,0);
          goto out;
	  
  /*
   * Error, object in file was of unknown type.
   * ------------------------------------------
   */

  err003: *jstat = -3;
          s6err("file_to_obj",*jstat,0);
          goto out;

  /*
   * Error, Homeogenous division by zero.
   * ------------------------------------------
   */

  err004: *jstat = -4;
          s6err("file_to_obj",*jstat,0);
          goto out;

  /*
   * Error, Cannot treat this dimension.
   * ------------------------------------------
   */

  err005: *jstat = -5;
          s6err("file_to_obj",*jstat,0);
          goto out;
  /*
   * Error in memory allocation.
   * ---------------------------
   */

  err101: *jstat = -101;
          s6err("file_to_obj",*jstat,0);
          goto out;

  /*
   * Exit file_to_obj, remember to close file if opened.
   * ---------------------------------------------------
   */

   out:

      return;
}


void curve_to_file(FILE *f,
		   struct SISLCurve *c1)
/*
*********************************************************************************
*
* Purpose : To write a B-spline curve to file.
*           File format hplabs.
*           Naming conventions SI/HP.
*
* Input   : to_file   - File name.
*           c1        - Pointer to curve.
*
* Output  : *           jstat       - status flag
*                       >= 0 - OK.
*                       <  0 - Error
*
* Written by :
*
********************************************************************************
*/
{
  int i,j;
  int linenum;

  fprintf(f,"$ This is a B-Spline curve\n");
  fprintf(f,"$ type: 0 is usual, 5 point, 6 analytic\n");
  fprintf(f,"0\n");

  /* order */
  fprintf(f,"$ order ik\n");
  fprintf(f,"%d\n",c1->ik);

  /* number of control vertices */
  fprintf(f,"$ number of control vertices in\n");
  fprintf(f,"%d\n",c1->in);

  /* dimension of geometry space  */
  fprintf(f,"$ dimension\n");
  fprintf(f,"%d\n",c1->idim);

  /* curve open/closed */
  fprintf(f,"$ curve open/closed\n");
  fprintf(f,"%d\n",c1->cuopen);

  /* nonrational, i.e. polynomial */
  fprintf(f,"$ rational or not\n");
  fprintf(f,"%d\n",(c1->ikind==2 || c1->ikind==4));

  /* knot vector */
  linenum = (c1->in + c1->ik)/4;
  fprintf(f,"$ knot vector\n");
  for (j=0; j < linenum; j++){
    for (i=0; i < 4; i++)
      fprintf(f,"%20.16g ",c1->et[j * 4 + i]);
    fprintf(f,"\n");
  }
  for (i = linenum * 4; i < (c1->in + c1->ik); i++)
    fprintf(f,"%20.16g ",c1->et[i]);
  fprintf(f,"\n");

  /* control vertices */
  fprintf(f,"$ control vertices\n");

  if (c1->ikind == 1 || c1->ikind == 3)
  {
     for ( i = 0; i < c1->in; i++ )
     {
	for ( j = 0; j < c1->idim; j++ ) fprintf(f,"%20.16g ",
						 c1->ecoef[i*c1->idim+j]);
	fprintf(f,"%g\n",(double)1.0);
     }
  }
  else
  {
     for ( i = 0; i < c1->in; i++ )
     {
	for ( j = 0; j < c1->idim+1; j++ ) fprintf(f,"%20.16g ",
						 c1->rcoef[i*(c1->idim+1)+j]);
	fprintf(f,"\n");
     }
  }

  /* instance matrix */
  fprintf(f,"$ instance matrix\n");
  if (c1->idim == 3)
  {
     fprintf(f,"%f %f %f %f\n",1.0,0.0,0.0,0.0);
     fprintf(f,"%f %f %f %f\n",0.0,1.0,0.0,0.0);
     fprintf(f,"%f %f %f %f\n",0.0,0.0,1.0,0.0);
     fprintf(f,"%f %f %f %f\n",0.0,0.0,0.0,1.0);
  }
  else if (c1->idim == 2)
  {
     fprintf(f,"%f %f %f\n",1.0,0.0,0.0);
     fprintf(f,"%f %f %f\n",0.0,1.0,0.0);
     fprintf(f,"%f %f %f\n",0.0,0.0,1.0);
  }
  else
  {
     fprintf(f,"%f %f\n",1.0,0.0);
     fprintf(f,"%f %f\n",0.0,1.0);
  }

  fprintf(f,"$ ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^\n");

  return;
}


void surface_to_file(FILE *f,
		     struct SISLSurf *surf)
/*
*********************************************************************************
*
* Purpose : To write a B-spline surface to file.
*           File format hplabs.
*           Naming conventions SI/HP.
*
* Input   : to_file   - File name.
*           surf      - Pointer to surface.
*
* Output  : *           jstat       - status flag
*                       >= 0 - OK.
*                       <  0 - Error
*
* Written by :
*
********************************************************************************
*/
{
  int i,j;
  int linenum;

  fprintf(f,"$ This is a B-Spline surface\n");
  fprintf(f,"$ type: 0 is usual, 5 point, 6 analytic\n");
  fprintf(f,"0\n");

  /* orders */
  fprintf(f,"$ orders in the two directions ik1, ik2\n");
  fprintf(f,"%d %d\n",surf->ik1,surf->ik2);

  /* numbers of control vertices */
  fprintf(f,"$ numbers of control vertices in the two directions in1, in2\n");
  fprintf(f,"%d %d\n",surf->in1,surf->in2);

  /* dimension of geometry space  */
  fprintf(f,"$ dimension\n");
  fprintf(f,"%d\n",surf->idim);

  /* surface open/closed */
  fprintf(f,"$ surface open/closed\n");
  fprintf(f,"%d %d\n",surf->cuopen_1,surf->cuopen_2);

  /* nonrational, i.e. polynomial in the two directions */
  fprintf(f,"$ rational or not\n");
  fprintf(f,"%d\n",(surf->ikind==2 || surf->ikind==4));

  /* knot vectors */
  linenum = (surf->in1 + surf->ik1)/4;
  fprintf(f,"$ knot vector in the first direction et1[]\n");
  for (j=0; j < linenum; j++){
    for (i=0; i < 4; i++)
      fprintf(f,"%20.16g ",surf->et1[j * 4 + i]);
    fprintf(f,"\n");
  }
  for (i = linenum * 4; i < (surf->in1 + surf->ik1); i++)
    fprintf(f,"%20.16g ",surf->et1[i]);
  fprintf(f,"\n");

  linenum = (surf->in2 + surf->ik2)/4;
  fprintf(f,"$ knot vector in the second direction et2[]\n");
  for (j=0; j < linenum; j++){
    for (i=0; i < 4; i++)
      fprintf(f,"%20.16g ",surf->et2[j * 4 + i]);
    fprintf(f,"\n");
  }
  for (i = linenum * 4; i < (surf->in2 + surf->ik2); i++)
    fprintf(f,"%20.16g ",surf->et2[i]);
  fprintf(f,"\n");

  /* control vertices */
  fprintf(f,"$ control vertices in the two directions \n");

  if (surf->ikind == 1 || surf->ikind == 3)
  {
     for ( i = 0; i < surf->in1 * surf->in2; i++ )
     {
	for ( j = 0; j < surf->idim; j++ )
	   fprintf(f,"%20.16g ",surf->ecoef[i*surf->idim+j]);
	fprintf(f,"%g\n",(double)1.0);
     }
  }
  else
  {
     for ( i = 0; i < surf->in1 * surf->in2; i++ )
     {
	for ( j = 0; j < surf->idim+1; j++ )
	   fprintf(f,"%20.16g ",surf->rcoef[i*(surf->idim+1)+j]);
	fprintf(f,"\n");
     }
  }

  /* instance matrix */
  fprintf(f,"$ instance matrix\n");
  if (surf->idim == 3)
  {
     fprintf(f,"%f %f %f %f\n",1.0,0.0,0.0,0.0);
     fprintf(f,"%f %f %f %f\n",0.0,1.0,0.0,0.0);
     fprintf(f,"%f %f %f %f\n",0.0,0.0,1.0,0.0);
     fprintf(f,"%f %f %f %f\n",0.0,0.0,0.0,1.0);
  }
  else if (surf->idim == 2)
  {
     fprintf(f,"%f %f %f\n",1.0,0.0,0.0);
     fprintf(f,"%f %f %f\n",0.0,1.0,0.0);
     fprintf(f,"%f %f %f\n",0.0,0.0,1.0);
  }
  else
  {
     fprintf(f,"%f %f\n",1.0,0.0);
     fprintf(f,"%f %f\n",0.0,1.0);
  }
  fprintf(f,"$ ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^\n");

  return;
}



int get_next_surface(FILE       *fp,
		     SISLSurf **qc)
{
  SISLSurf *qs=NULL;         /* Pointer to struct SISLSurf  */
  SISLObject *qo=NULL;       /* Pointer to struct SISLObject */
  //  SISLObject *qo=NULL;       /* Pointer to struct SISLObject */
  char string[80];       /* character to read $ which identifies a comment */
  int value;             /* type of B-Spline and numbers */
  int kcopy=1;           /* Flag */
  int kcuopen1;          /* Open/closed flag for surf  */
  int kcuopen2;          /* Open/closed flag for surf  */
  int kdim;              /* dimension of target space */
  int kdim1;             /* dimension + 1             */
  int kk1;               /* Order of curve or in 1st parameter direction of surface */
  int kk2;               /* Order in 2nd parameter direction of surface */
  int kn1;               /* # of vertices of curve or in 1st parameter direction of surf */
  int kn2;               /* # of vertices in 2nd parameter direction of surface */
  double *st1=NULL;      /* Knot vector of curve or in 1st parameter direction of surface */
  double *st2=NULL;      /* Knot vector in 2nd parameter direction of surface */
  double *scoef=NULL;    /* Vertices of curve or surface */
  double dummy;          /* trashbin for fourth value in vertices */
  double help[4];        /* Help array in computing transformated vertic. */
  double matrix[16];     /* instance matrix */
  double *sptr=NULL;     /* Dummy pointer */
  int ki, kj;            /* Counters */
  int *ind_hex=0;
  int hex=FALSE;
  int *ind_2 = (int*)&dummy;
  std::vector<double> seg_u;
  std::vector<double> seg_v;
  std::vector<int> prio_u;
  std::vector<int> prio_v;
  std::vector<double> pbox;
  int stat = 0;


  /*
   * Read data from file.
   * --------------------
   */
  read_non_comment(fp,string);
  if (string[0] == EOF) goto err001;
  sscanf(string,"%d",&value);

  read_non_comment(fp,string);
  sscanf(string,"%d",&kk1);
  fscanf(fp,"%d",&kk2);

  read_non_comment(fp,string);
  sscanf(string,"%d",&kn1);
  fscanf(fp,"%d",&kn2);


  /*
   * dimension
   * ---------
   */

  read_non_comment(fp,string);
  sscanf(string,"%d",&kdim);
  kdim1 = kdim;

  /*
   * surface open/closed
   * ------------------
   */

  read_non_comment(fp,string);
  sscanf(string,"%d",&kcuopen1);
  fscanf(fp,"%d",&kcuopen2);

  /*
   * rational B-Spline or not
   * ------------------------
   */

  read_non_comment(fp,string);
  sscanf(string,"%d",&value);
  /* if (value != 0) goto err002;

  fscanf(fp,"%d",&value); */
  if (value) kdim1++;

  /*
   * Allocate space for first knot vector
   * ------------------------------------
   */

//   st1 = newarray(kk1+kn1,DOUBLE);
  st1 = new double[kk1+kn1];
  if (!st1) goto err101;

  /*
   * Read in first knot vector.
   * --------------------------
   */

  read_non_comment(fp,string);
  if (hex)
      {
	  ind_hex = (int *)st1;
	  sscanf(string,"%x",ind_hex);
	  fscanf(fp,"%x",ind_hex+1);
	  for (ind_hex=(int*)(st1+1),ki=1;ki < kk1+kn1;ki++,ind_hex+=2)
	      fscanf(fp,"%x %x",ind_hex,ind_hex+1);
      }
  else
      {

	  sscanf(string,"%lf",st1);
	  for (sptr = st1+1,ki = 1; ki < kn1+kk1; sptr++,ki++)
	      fscanf(fp,"%lf",sptr);
      }
  /*
   * Allocate space for second knot vector
   * ------------------------------------
   */

//   st2 = newarray(kk2+kn2,DOUBLE);
  st2 = new double[kk2+kn2];
  if (!st2) goto err101;

  /*
   * Read in second knot vector.
   * --------------------------
   */

  read_non_comment(fp,string);
  if (hex)
      {
	  ind_hex = (int *)st2;
	  sscanf(string,"%x",ind_hex);
	  fscanf(fp,"%x",ind_hex+1);
	  for (ind_hex=(int*)(st2+1),ki=1;ki < kk2+kn2;ki++,ind_hex+=2)
	      fscanf(fp,"%x %x",ind_hex,ind_hex+1);
      }
  else
      {

	  sscanf(string,"%lf",st2);
	  for (sptr = st2+1,ki = 1; ki < kn2+kk2; sptr++,ki++)
	      fscanf(fp,"%lf",sptr);
      }

  /*
   * Allocate space for control vertices
   * -----------------------------------
   */

//   scoef = newarray(kn1*kn2*kdim1,DOUBLE);
  scoef = new double[kn1*kn2*kdim1];
  if (!scoef) goto err101;

  /*
   * Read in control vertices.
   * -------------------------
   */

  read_non_comment(fp,string);
  if (hex)
      {
	  ind_hex = (int *)scoef;
	  sscanf(string,"%x",ind_hex);
	  fscanf(fp,"%x",ind_hex+1);
	  if (kdim > 1)
	      {
		  ind_hex +=2;
		  fscanf(fp,"%x %x",ind_hex,ind_hex+1);
		  if (kdim > 2)
		      {
			  ind_hex +=2;
			  fscanf(fp,"%x %x",ind_hex,ind_hex+1);
		      }
	      }
	  fscanf(fp,"%x %x",ind_2,ind_2+1);
	  if (value)
	      {
		  ind_hex +=2;
		  scoef[kdim] = dummy;
	      }
	  for (ki = kdim1; ki < kn1*kn2*kdim1; ki += kdim1)
	      {
		  ind_hex +=2;
		  fscanf(fp,"%x %x",ind_hex,ind_hex+1);
		  if (kdim > 1)
		      {
			  ind_hex +=2;
			  fscanf(fp,"%x %x",ind_hex,ind_hex+1);
			  if (kdim > 2)
			      {
				  ind_hex +=2;
				  fscanf(fp,"%x %x",ind_hex,ind_hex+1);
			      }
		      }
		  fscanf(fp,"%x %x",ind_2,ind_2+1);
		  if (value)
		      {
			  ind_hex +=2;
			  scoef[ki+kdim] = dummy;
		      }
	      }

      }
  else
      {
	  sscanf(string,"%lf",scoef);
	  if (kdim > 1) fscanf(fp,"%lf",&scoef[1]);
	  if (kdim > 2) fscanf(fp,"%lf",&scoef[2]);
	  fscanf(fp,"%lf",&dummy);
	  if (value) scoef[kdim] = dummy;
	  for (ki = kdim1; ki < kn1*kn2*kdim1; ki += kdim1)
	      {
		  fscanf(fp,"%lf",&scoef[ki]);
		  if (kdim > 1) fscanf(fp,"%lf",&scoef[ki+1]);
		  if (kdim > 2) fscanf(fp,"%lf",&scoef[ki+2]);
		  fscanf(fp,"%lf",&dummy);
		  if (value) scoef[ki+kdim] = dummy;
	      }

      }
  /*
   * transformation data
   * -------------------
   */

  read_non_comment(fp,string);
  sscanf(string,"%lf",matrix);
  for (ki = 1; ki < (kdim+1)*(kdim+1); ki ++)
      fscanf(fp,"%lf",&matrix[ki]);

  /*
   * transform vertices with instance matrix
   * ---------------------------------------
   */


  /*
   * test whether matrix is unit matrix
   * ----------------------------------
   */

  if (kdim      == 1           &&
      matrix[0] == (double)1.0 &&
      matrix[1] == (double)0.0 &&
      matrix[2] == (double)0.0 &&
      matrix[3] == (double)1.0 )
      {
      }
  else if (kdim      == 2           &&
	   matrix[0] == (double)1.0 &&
	   matrix[1] == (double)0.0 &&
	   matrix[2] == (double)0.0 &&
	   matrix[3] == (double)0.0 &&
	   matrix[4] == (double)1.0 &&
	   matrix[5] == (double)0.0 &&
	   matrix[6] == (double)0.0 &&
	   matrix[7] == (double)0.0 &&
	   matrix[8] == (double)1.0)
      {
      }
  else if (kdim       == 3           &&
	   matrix[0]  == (double)1.0 &&
	   matrix[1]  == (double)0.0 &&
	   matrix[2]  == (double)0.0 &&
	   matrix[3]  == (double)0.0 &&
	   matrix[4]  == (double)0.0 &&
	   matrix[5]  == (double)1.0 &&
	   matrix[6]  == (double)0.0 &&
	   matrix[7]  == (double)0.0 &&
	   matrix[8]  == (double)0.0 &&
	   matrix[9]  == (double)0.0 &&
	   matrix[10] == (double)1.0 &&
	   matrix[11] == (double)0.0 &&
	   matrix[12] == (double)0.0 &&
	   matrix[13] == (double)0.0 &&
	   matrix[14] == (double)0.0 &&
	   matrix[15] == (double)1.0)
      {
      }
  else if (1)
      {
      }

  else
      {
	  if (kdim > 3 || kdim < 1) goto err005;
	  if (!value)
	      {
		  for (ki = 0; ki < kn2*kn1*kdim; ki += kdim)
		      {
			  kdim1 = kdim + 1;
			  for (kj = 0; kj < kdim1; kj++)
			      help[kj]  = s6scpr(&scoef[ki],
						 &matrix[kj*kdim1],kdim)
				  + matrix[kj*kdim1 + kdim];

			  if (help[kdim] == 0.0) goto err004;

			  for (kj = 0; kj < kdim; kj++)
			      scoef[ki + kj] = help[kj]/help[kdim];
		      }
	      }
	  else goto err002;
      }

  // Read segmentation and parameter box for surface

  // 			  stat = readSegments(fp,seg_u,seg_v,prio_u,prio_v,pbox);
  // 			  if (stat < 0)
  // 			      goto error;

  /*
   * Create an instance of SISLSurf.
   * ---------------------------
   */

  qs = newSurf(kn1,kn2,kk1,kk2,st1,st2,scoef,(value)?2:1,kdim,kcopy);
  if (!qs) goto err101;
  else
      {
// 	  // Set segmentation values and parameter box

// 	  if (pbox.size() == 4)
// 	      SI_setparbox(qs, pbox.begin());

// 	  if (seg_u.size() > 0 || seg_v.size() > 0)
// 	      {
// 		  stat = SI_set_seg(qs, seg_u.size(), seg_u.begin(),
// 				    prio_u.begin(),seg_v.size(),
// 				    seg_v.begin(), prio_v.begin());
// 		  if (stat < 0)
// 		      goto error;

// 	      }
	  /*
	   * Create an instace of Object and connect surface to the Object.
	   * --------------------------------------------------------------
	   */

	  qs->cuopen_1 = kcuopen1;
	  qs->cuopen_2 = kcuopen2;

	  qo = newObject(SISLSURFACE);
	  if (!qo) goto err101;
	  else
	      {
		  qo->s1 = qs;
		  *qc = qs; //    = qo;
	      }
      }

  /*
   * Remember to free allocated space.
   * ---------------------------------
   */

  stat = 0;
  goto out;

 err001: stat = -1;
  s6err("file_to_obj",stat,__LINE__);
  goto out;

  /*
   * Error, object in file was of type rational.
   * -------------------------------------------
   */

 err002: stat = -2;
  s6err("file_to_obj",stat,__LINE__);
  goto out;

//   /*
//    * Error, object in file was of unknown type.
//    * ------------------------------------------
//    */

//  err003: stat = -3;
//   s6err("file_to_obj",stat,__LINE__);
//   goto out;

  /*
   * Error, Homeogenous division by zero.
   * ------------------------------------------
   */

 err004: stat = -4;
  s6err("file_to_obj",stat,__LINE__);
  goto out;

  /*
   * Error, Cannot treat this dimension.
   * ------------------------------------------
   */

 err005: stat = -5;
  s6err("file_to_obj",stat,__LINE__);
  goto out;
  /*
   * Error in memory allocation.
   * ---------------------------
   */

 err101: stat = -101;
  s6err("file_to_obj",stat,__LINE__);
  goto out;

  /*
   * Exit file_to_obj, remember to close file if opened.
   * ---------------------------------------------------
   */

 out:
  if (st1) delete [] st1;
  if (st2) delete [] st2;
  if (scoef) delete [] scoef;

  return stat;
}


// Read sisl sfs from file.
// If returned value is negative something went wrong.
int get_sisl_surfaces(FILE* fp,
		      std::vector<shared_ptr<SISLSurf> >& sisl_sfs)
{
//   read_non_comment(fp,string);
  int kstat = 0;
  while (kstat == 0)
    {
      SISLSurf* sisl_sf = NULL;
      kstat = get_next_surface(fp, &sisl_sf);
      if (kstat < 0)
	{
	  return kstat;
	}
      else
	sisl_sfs.push_back(shared_ptr<SISLSurf>(sisl_sf));


      char string[80]; /* character to read $ which identifies a comment */
      read_non_comment(fp,string);
      if (string[0] == EOF)
	break; // Hmm, but what about next surface ... Will we not
	       // miss the first comment?
    }

  return 0;
}
