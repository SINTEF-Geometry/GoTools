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

#ifndef _SISL_CODE_H
#define _SISL_CODE_H

//===========================================================================
// SISL DEFINEs
//===========================================================================
/// \cond SISL
/* Name of geometry objects. Used in branching. */
#define SISLPOINT    0
#define SISLCURVE    1
#define SISLSURFACE  2


//===========================================================================
// SISL structs
//===========================================================================

//===========================================================================
typedef struct SISLdir
//===========================================================================
{
  int igtpi;			/* 0 - The direction of the surface or curve
			               is not greater than pi in any
			               parameter direction.
			           1 - The direction of the surface or curve
			               is greater than pi in the first
			               parameter direction.
			           2 - The direction of the surface is greater
			               than pi in the second parameter
			               direction. 			     */
  double *ecoef;		/* The coordinates to the center of the cone.*/
  double aang;			/* The angle from the center whice describe the
			           cone.				     */
  double *esmooth;		/* Coordinates of object after smoothing.    */
} SISLdir;
 /* The following structure contains 3 different boxes. The
    first box is the plain box given by the coefficients of the
    object. The second box is expanded with the half of a given
    tolerance. The third box is expanded with half the tolerance
    in the inner and for the vertices at the edges/endpoints
    a distance of half the tolerance is removed. The minimum and
    maximum values of the boxes are given by the arrays
    e2min[0] - e2min[2] and e2max[0] - e2max[2]. The tolerances used
    when making the boxes are stored in etol[0] - etol[2]. etol[0]
    will always be equal to zero. If a box is made, the pointers
    belonging to this box points to arrays, otherwise they point
    to SISL_NULL.                                                       */

//===========================================================================
typedef struct SISLbox
//===========================================================================
{
  double *emax;			/* The minimum values to the boxes.	     */
  double *emin;			/* The maximum values to the boxes.	     */
  int imin;			/* The index of the min coeff (one-dim case) */
  int imax;			/* The index of the max coeff (one-dim case) */

  double *e2max[3];		/* The minimum values dependant on tolerance.*/
  double *e2min[3];		/* The maximum values dependant on tolerance.*/
  double etol[3];		/* Tolerances of the boxes.                  */
} SISLbox;

//===========================================================================
typedef struct SISLCurve
//===========================================================================
{
  int ik;			/* Order of curve.                           */
  int in;			/* Number of vertices.                       */
  double *et;			/* Pointer to the knotvector.                */
  double *ecoef;		/* Pointer to the array containing vertices. */
  double *rcoef;		/*Pointer to the array of scaled vertices if
				  rational.  */
  int ikind;			/* Kind of curve
	                           = 1 : Polynomial B-spline curve.
	                           = 2 : Rational B-spline curve.
	                           = 3 : Polynomial Bezier curve.
	                           = 4 : Rational Bezier curve.             */
  int idim;			/* Dimension of the space in which the curve
				   lies.      */
  int icopy;			/* Indicates whether the arrays of the curve
				   are copied or referenced by creation of the
				   curve.
	                           = 0 : Pointer set to input arrays.
			           = 1 : Copied.
	                           = 2 : Pointer set to input arrays,
				         but are to be treated as copied.   */
  SISLdir *pdir;		/* Pointer to a structur to store curve
				   direction.      */
  SISLbox *pbox;		/* Pointer to a structur to store the
				   surrounded boxes. */
  int cuopen;			/* Open/closed flag.                         */
} SISLCurve;

//===========================================================================
typedef struct SISLIntcurve
//===========================================================================
{
  int ipoint;			/* Number of points defining the curve.      */
  int ipar1;			/* Number of parameter directions of first
				   object.                                   */
  int ipar2;			/* Number of parameter directions of second
				 * object.                                   */
  double *epar1;		/* Pointer to the parameter-values of the
				   points
	                           in the first object.                      */
  double *epar2;		/* Pointer to the parameter-values of the
				   points
	                           in the second object. If one of the objects
	                           is an analytic curve or surface epar2 points
	                           to nothing.                               */
  SISLCurve *pgeom;		/* Pointer to the intersection curve in the
				   geometry space. If the curve is not
				   computed, pgeom points to nothing.       */
  SISLCurve *ppar1;		/* Pointer to the intersection curve in the
				   parameter plane of the first object. If
				   the curve is not computed, ppar1 points
				   to nothing.                              */
  SISLCurve *ppar2;		/* Pointer to the intersection curve in the
				   parameter plane of the second object. If
				   the curve is not computed, ppar2 points
				   to nothing.                              */
  int itype;			/* Kind of curve.
	                           = 1 : Straight line.
	                           = 2 : Closed loop. No singularities.
	                           = 3 : Closed loop. One singularity.
				         Not used.
	                           = 4 : Open curve. No singularity.
	                           = 5 : Open curve. Singularity at the
	                                 beginning of the curve.
	                           = 6 : Open curve. Singularity at the end
	                                 of the curve.
	                           = 7 : Open curve. Singularity at the
				         beginning  and end of the curve.
	                           = 8 : An isolated singularity. Not used.
				   = 9 : The curve is exact, pgeom and either
				   	 ppar1 or ppar2 is set.      */

  int pretop[4];		/* Pretopology */
} SISLIntcurve;

//===========================================================================
typedef struct SISLIntsurf
//===========================================================================
{
  int ipoint;			/* Number of points defining the curve.    */
  int ipar;			/* Number of parameter directions of       */
  double *epar;		        /* Pointer to the parameter-values of the
				   points, dimension: ipoint*ipar          */
  int *const_par;               /* Constant parameter direction between
				   two points in epar.                     */
} SISLIntsurf;

//===========================================================================
typedef struct SISLIntpt
//===========================================================================
{
  int ipar;			/* Number of parameter directions in
				 * intersection problem.                     */
  double *epar;			/* Parametervalues of point, possibly in two
				 * objects.                                  */
  double adist;			/* Distance between the objects in this point.
				 * tdist is used in closest point problems.  */
  struct SISLIntpt *pcurve;	/* Not used, kept for compatibility with old
			           version on the structure.*/
  int iinter;			/* = 1 ORDINARY MAIN POINT
				   = 2 SINGULAR MAIN POINT
				   = 3 TRIM MAIN POINT
		                   = -1 ORDINARY HELP POINT
				   = -2 SINGULAR HELP POINT
				   = -3 TRIM HELP POINT */
  struct SISLIntpt **pnext;	/* Pointers to next points in each curve
				 * chain.                                   */
  int *curve_dir;		/* An array of curve directions + from - to
				 * this point.                              */
  int no_of_curves;		/* Number of curves containing this point.  */
  int no_of_curves_alloc;	/* The size of the arrays allocated         */
  int *left_obj_1;		/* Pretopology information, one for each
				 * curve.                                   */
  int *left_obj_2;		/* Pretopology information, one for each
				 * curve.                                   */
  int *right_obj_1;		/* Pretopology information, one for each
				   curve.                                   */
  int *right_obj_2;		/* Pretopology information, one for each
				   curve.                                   */
  int size_1;			/* Size of geo_data_1                       */
  int size_2;			/* Size of geo_data_2                       */
  double *geo_data_1;		/* Containing geometric info first object   */
  double *geo_data_2;		/* Containing geometric info second object  */
  /*  double  geo_aux[3]; Containing auxiliary geo info, see sh6idput*/
  double geo_track_3d[10];	/* To store intersection curve info */
  double geo_track_2d_1[7];
  double geo_track_2d_2[7];
  int edge_1;                   /* Edge flag for topology           */
  int edge_2;
  int marker;                   /* Help attribute when creating lists  */
  int evaluated;                /* Help attribute when creating tracks */
  struct SISLTrimpar *trim[2];          /* Used if pt is in trim curve. */
  int iside_1;			/* Left/right evaluator flag.  -1,0+ */
  int iside_2;			/* Left/right evaluator flag.  -1,0+*/
} SISLIntpt;

//===========================================================================
typedef struct SISLIntlist
//===========================================================================
{
  SISLIntpt *pfirst;		/* Pointer to first point in list. */
  SISLIntpt *plast;		/* Pointer to last point in list.  */
  int ind_first;		/* Index pointer in pfirst         */
  int ind_last;			/* Index pointer in plast          */
  int itype;			/* Status of curve-segment.
                               = 0 : open curve, no singularities.
                               = 1 : closed curve, no singularities.
                               = 2 : more than two curves meet at start point.
                               = 3 : more than two curves meet at end point.
                               = 4 : more than two curves meet at start
                                     and end point.
                               = 5 : isolated singularity.
                               = 6 : touching area of surface.            */
  int inumb;			/* Number of points in the list.          */
  int pretop[4];		/* Pretopology */
} SISLIntlist;



//===========================================================================
typedef struct SISLIntdat
//===========================================================================
{
  SISLIntpt **vpoint;
  int ipoint;
  int ipmax;
  SISLIntlist **vlist;
  int ilist;
  int ilmax;
} SISLIntdat;



//===========================================================================
typedef struct SISLSurf
//===========================================================================
{
  int ik1;			/* Order of surface in first parameter
				   direction.       */
  int ik2;			/* Order of surface in second parameter
				   direction.      */
  int in1;			/* Number of vertices in first parameter
				   direction.     */
  int in2;			/* Number of vertices in second parameter
				   direction.    */
  double *et1;			/* Pointer to knotvector in first parameter
				   direction.  */
  double *et2;			/* Pointer to knotvector in second parameter
				   direction. */
  double *ecoef;		/* Pointer to array of vertices of surface. */
  double *rcoef;		/* Pointer to the array of scaled vertices
				   if surface is rational. */
  int ikind;			/* Kind of surface
	                           = 1 : Polynomial B-spline tensor-product
				         surface.
	                           = 2 : Rational B-spline tensor-product
				         surface.
	                           = 3 : Polynomial Bezier tensor-product
				         surface.
	                           = 4 : Rational Bezier tensor-product
				         surface.                           */
  int idim;			/* Dimension of the space in which the surface
				   lies.    */
  int icopy;			/* Indicates whether the arrays of the surface
				   are copied or referenced by creation of
				   the surface.
	                           = 0 : Pointer set to input arrays.
			           = 1 : Copied.
	                           = 2 : Pointer set to input arrays,
				         but are to be treated as copied.               */
  SISLdir *pdir;		/* Pointer to a structur to store surface
				   direction.    */
  SISLbox *pbox;		/* Pointer to a structur to store the
				   surrounded boxes. */
  int use_count;                /* use count so that several tracks can share
				   surfaces, no internal use */
 int cuopen_1;                  /* Open/closed flag, 1. par directiion */
 int cuopen_2;                  /* Open/closed flag. 2. par direction  */
} SISLSurf;

typedef struct SISLTrack
{
  SISLSurf *psurf_1;		/* Pointer to first surface in intersection */
  SISLSurf *psurf_2;		/* Pointer to second surface in intersection */
  SISLCurve *pcurve_3d;		/* Pointer to 3D support curve. */
  SISLCurve *pcurve_2d_1;	/* Pointer to 2D support curve in first
				   parameter space. */
  SISLCurve *pcurve_2d_2;	/* Pointer to 2D support curve in second
				   parameter space. */
  int ideg;			/* Type of track.
				    = 0, Bspline vs Bspline
				= 1, Bspline vs Plane
				= 2, Bspline vs Quadric surface
				= 1001 Bspline vs Torus surface
				= 1003 Bspline silhouette line, parallel
				  projection
				= 1004 Bspline silhouette line, perspective
				  projection
				= 1005 Bspline silhouette line, circular
			          projection */

  double eimpli[16];		/* Description of the implicit surface */
  int turned;			/* Connection between the direction of the
				   support curve and the cross product
				   between the two surface normals.
				= 0, same direction
				= 1, oposite direction */
  int exact;                    /* Flag if curve is exact */
  int pretop[4];		/* Pretopology */
  int sing_start;               /* Singular start end point markers */
  int sing_end;
} SISLTrack;

//===========================================================================
typedef struct SISLPoint
//===========================================================================
{
  double ec[3];
  int idim;			/* The dimension the point lies in           */
  double *ecoef;		/* Pointer to the array containing the
				   coordinates */
  int icopy;			/* Indicates whether the arrays of the point
				   are copied or referenced by creation of
				   the point.
				   = 0 : Pointer set to input arrays.
				   = 1 : Copied.
				   = 2 : Pointer set to input arrays,
				         but are to be treated as copied.    */

  SISLbox *pbox;		/*Pointer to a structur to store the boxes.  */
} SISLPoint;

//===========================================================================
typedef struct SISLObject
//===========================================================================
{
  int iobj;			/* Integer indicates which kind of geometric
				   object is contained in a particular
				   instance of the structure.
		                   = 1 (SISLCurve) - curve.
		                   = 2 (SURFACE) - tensor-product surface.   */
  SISLPoint *p1;		/* Pointer to a point (instance of Point).   */
  SISLCurve *c1;		/* Pointer to a curve
				 * (instance of SISLCurve).                  */
  SISLSurf *s1;			/* Pointer to a surface
				 * (instance of SISLSurf).                   */
  struct SISLObject *o1;	/* Pointer to parent object
				 * (instance of Object).                     */
  struct SISLObject *edg[4];	/* Pointer to objects edges
				 * (instance of Object).                     */
  struct SISLObject *psimple;	/* Indicates if object/object intersection
				 * is simple case. */
} SISLObject;

//===========================================================================
typedef struct SISLPtedge
//===========================================================================
{
  SISLIntpt *ppt;		/* Pointer to intersection points.      */
  struct SISLPtedge *pnext;	/* Pointer to next element in the list. */
} SISLPtedge;

//===========================================================================
typedef struct SISLEdge
//===========================================================================
{
  int iedge;			/* Number of edges/endpoints of object.      */
  int ipoint;			/* Number of intersection points found on
				 * the edges.                                */
  SISLPtedge **prpt;		/* Array containing lists of pointers to the
			         * intersections at the edges.               */
} SISLEdge;



//===========================================================================
// SISL constructors/destructors
//===========================================================================
SISLObject   *newObject(int);
SISLCurve    *newCurve(int,int,double *,double *,int,int,int);
void freeCurve(SISLCurve *);
void freeObject(SISLObject *);
void freeIntdat(SISLIntdat *pintdat);
void sh1761 (SISLObject * po1, SISLObject * po2, double aepsge, 
	     SISLIntdat ** pintdat, int *jstat);
SISLPoint *newPoint (double *ecoef, int idim, int icopy);
SISLSurf     *newSurf(int,int,int,int,double *,double *,double *,int,int,int);
void freeIntcrvlist(SISLIntcurve **,int);
void freeIntcurve(SISLIntcurve *pintc);
void freeSurf(SISLSurf *);

//===========================================================================
// SISL functions directly used by 'sisl_dependent'
//===========================================================================

void s1220(double *et,int ik,int in,int *ileft,
	   double ax,int ider,double ebder[],int *jstat);
void s1310(SISLSurf *,SISLSurf *,SISLIntcurve *,double,double,int,int,int *);
void s1314(SISLSurf *,double *,double *,int,double,double,double,
	   SISLIntcurve *,int,int,int *);
void s1316(SISLSurf *ps1,double *epoint,double *edirec,double aradiu,
	   int idim,double aepsco,double aepsge,double amax,
	   SISLIntcurve *pintcr,int icur,int igraph,int *jstat);
void s1421(SISLSurf *,int,double [],int *,int *,double [],double [],int *);
void s1770(SISLCurve *,SISLCurve *,double,double,double,double,double,
	   double,double,double *,double *,int *);
void s1785(SISLCurve *pcurve,SISLSurf *psurf,double aepsge,
	   double epar1[],double epar2[],int icur,int *jstat);
void s1851(SISLSurf *,double [],double [],int,double,double,
	   int *,double **,int *,SISLIntcurve ***,int *);
void s1853(SISLSurf *ps1,double epoint[],double edirec[],double aradius,
	   int idim,double aepsco,double aepsge,int *jpt,double **gpar,
	   int *jcrv,SISLIntcurve ***wcurve,int *jstat);
void s1856(SISLSurf *ps1,double epoint[],double edir[],int idim,
	   double aepsco,double aepsge,int *jpt,double **gpar,
	   int *jcrv,SISLIntcurve ***wcurve,int *jstat);
void s1859(SISLSurf *,SISLSurf *,double,double,
	   int *,double **,double **,int *,SISLIntcurve ***,int *);

void sh1857(SISLCurve *,SISLCurve *,double,double,int,int *,SISLTrack ***,
	    int *,double **,double **,int **,int *,SISLIntcurve ***,int *);

void sh1858(SISLSurf *,SISLCurve *,double,double,int,int *,SISLTrack ***,
	    int *,double **,double **,int **,int *,SISLIntcurve ***,int *);


void
   s1871(SISLCurve *pc1, double *pt1, int idim, double aepsge,
	 int *jpt,double **gpar1,int *jcrv,SISLIntcurve ***wcurve,int *jstat);
double s6scpr(double e1[],double e2[],int idim);
void s6err(const char *,int,int);

/// \endcond SISL

#endif // _SISL_CODE_H

