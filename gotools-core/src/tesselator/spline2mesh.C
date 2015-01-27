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

#include "GoTools/tesselator/spline2mesh.h"
#include "GoTools/utils/errormacros.h"
#include "GoTools/tesselator/2dpoly_for_s2m.h"
#include "GoTools/geometry/SplineCurve.h"


#include <fstream>



//======================================================================
//
// 040607: This is very similar to "marching squares", so maybe 
//         the routines should be based on some common code...
//
//         NB! Vertices added are only the new ones! We assume that
//         the regular non-trimmed net of vertices are already
//         present.
// 040614: No, we do this in the function.
//
// 040614: Triangles are oriented ccw.
//         The param 'n' denotes number of "quads", not corner nodes.
//
//======================================================================


#include <functional>

#ifdef __BORLANDC__
#include <algorithm>

using std::fabs;
using std::lower_bound;
using std::printf;
using std::puts;
using std::sqrt;
#endif

using namespace Go;
using std::vector;




 

#define S2M_VERBOSEnot		// May be helpful during debugging.

//#define DBG			// Produces all sorts of output, including various matlab scripts.
				// 100214: Defined by the build system instead.

#define CLC_DBGq			// Corner list construction. Note: c? are still created if DBG is set.

#define CORNER_SPLITTING	// Initial splitting of quads (and secondary triangles) containing corners.
#define SPLIT_NON_CORNER_QUADS	//
#define SPLIT_TRIANGLES_AT_OUTER_CURVE
#define SPLIT_TRIANGLES_AT_INNER_CURVES

#ifdef DBG			// 100222: A construction for combining debugging facilities, readability and
				//         optimization in the release-version.
#  define DBG_FLAG ,dbg
#else
#  define DBG_FLAG
#endif

#define SPLIT_LARGE_TRIANGLESqq	// 100225: Experimental

#define EVAL_SRF_NOT_INTERP	// 100225: Experimental. Evaluating surface instead of interpolating triangle corners.






namespace Go
{

  const bool s2m_with_boundary=false; // true;	// 090117: Enabling this makes "point-inside-contour"-tests true for
  //         points on the contour. At the moment this is done by
  //         testing the point against all segments of a discretized
  //         contour. This kind of negates the speed-up gained by the
  //         current fast implementation of 'point_inside_contour'. This
  //         can be fixed by doing the same trick to eliminate
  //         irrelevant contour segments (without having to loop over
  //         them all as is now done.)
  // 	   Will the same exercise have to be done to 'split_triangle'?!






  //==============================================================================================================
  //
  // 090204: Now also testing for non-degeneracy before actually adding it.
  //
  //==============================================================================================================
  
  void add_triangle(const vector<Vector2D> &vert_p,
		    vector<int> &mesh,
		    const int c1, const int c2, const int c3
#ifdef DBG
		    , const bool dbg = false
#endif
    )
  {
#ifndef DBG
    const bool dbg = false;
#endif

    if (!degenerate_triangle(vert_p[c1], vert_p[c2], vert_p[c3] DBG_FLAG))
      {
	mesh.push_back(c1);
	mesh.push_back(c2);
	mesh.push_back(c3);
      }
   else
     if (dbg)
       {
	 printf("    add_triangle(): Not adding triangle. Degenerate.\n");
	 printf("\n");
	 printf("    hold on; plot(%f, %f, 'ro', 'markersize', 10, 'linewidth', 2);\n", vert_p[c1][0], vert_p[c1][1]);
	 printf("    plot(%f, %f, 'go', 'markersize', 10, 'linewidth', 2);\n", vert_p[c2][0], vert_p[c2][1]);
	 printf("    plot(%f, %f, 'bo', 'markersize', 10, 'linewidth', 2); hold off\n", vert_p[c3][0], vert_p[c3][1]);
	 printf("\n");
	 // const bool tmp = degenerate_triangle(vert_p[c1], vert_p[c2], vert_p[c3] DBG_FLAG);
       }
 }






#ifdef DBG

  //==============================================================================================================
  //
  // 090216: For debugging.
  //
  //==============================================================================================================
  
  bool same_point(const Vector2D &c1, const Vector2D &d1, const double tol)
  {
    return (fabs(c1[0]-d1[0])<tol) && (fabs(c1[1]-d1[1])<tol);
  }

  bool same_triangle(const Vector2D &c1, const Vector2D &c2, const Vector2D &c3, 
		     const Vector2D &d1, const Vector2D &d2, const Vector2D &d3,
		     const double tol)
  {
    return ( ((same_point(c1, d1, tol)) && (same_point(c2, d2, tol)) && (same_point(c3, d3, tol))) ||
	     ((same_point(c1, d1, tol)) && (same_point(c2, d3, tol)) && (same_point(c3, d2, tol))) ||
	     ((same_point(c1, d2, tol)) && (same_point(c2, d1, tol)) && (same_point(c3, d3, tol))) ||
	     ((same_point(c1, d2, tol)) && (same_point(c2, d3, tol)) && (same_point(c3, d1, tol))) ||
	     ((same_point(c1, d3, tol)) && (same_point(c2, d1, tol)) && (same_point(c3, d2, tol))) ||
	     ((same_point(c1, d3, tol)) && (same_point(c2, d2, tol)) && (same_point(c3, d1, tol)))    );
  }

#endif






  //--------------------------------------------------------------------------------------------------------------
  //
  // This subroutine assumes that either one, two or three points are inside, and the rest outside, so that two
  // intersections (between "grid lines" and trim curve) are to be found, and either one, two or three triangles
  // are to be added, respectively.
  //
  // 081208: Note that 'srf' is not evaluated to create new vertices, the old ones are simply interpolated! This
  //         could be improved!
  //
  // 081208: Hmm... Can this interpolation result in non-normalized normals?! Think so! But is it specified that
  //         they should be normalized?! Safest thing would probably be to normalize them... Not doing it for
  //         split_quad, but for split_triangle...
  //
  // 081208: Adding 'vert_p' for use by later 'split_triangle'-calls. This will hold parameter pairs.
  //
  // 090117: Adding 'with_boundary' so that the caller may decide whether or not the boundary (with an
  //         "epsilon-tube" around it) should be defined as "inside".
  //
  //--------------------------------------------------------------------------------------------------------------

  void split_quad(const int i0, const int j0,
		  const int i1, const int j1,	// May not be used, if so, -1 is used.
		  const int i2, const int j2,	// May not be used, if so, -1 is used.
		  const vector<int> &si0, const vector<int> &sj0,
		  const vector<int> &si1, const vector<int> &sj1,

		  // 090115: All the vectors above have two elements each, s?0[i] and s?1[i] are the two end
		  //         points of two of the sides of some quad, that the contour is assumed to intersect.
		  //         (The contour is given in 'contour', and 'vert', together. I think.)
		  //         (i?, j?) (actually, it makes more sense to say (j?, i?)...) specifies corners of
		  //         the quad. Which corners, is dependent on the configuration in question... (sigh.)

		  shared_ptr<ParamSurface> srf,
		  vector< Vector3D > &vert, vector<Vector2D> &vert_p,
		  vector< int > &bd, vector< Vector3D > &norm,
		  //vector< Vector3D > &col,
		  vector<int> &mesh, const int dn, const int dm,
		  //vector< Vector3D > &extra_v,
		  const vector< Vector3D > &trim_curve_p,
		  const vector<int> &contour,
		  const bool with_boundary, const bool dbg=false)
  {
    //
    // This one always get two existing points and two segments between two pairs of points in, and intersect
    // the trimming curve twice against these segments for the two remaining points. The triangles are not
    // always oriented the same way, but that shouldn't matter too much. [Could easily be fixed by a switch
    // among the arguments...]
    //

    int i;

    //
    // Do the intersections.
    //

    const RectDomain dom = srf->containingDomain();
    const double u0 = dom.umin();
    const double u1 = dom.umax();
    const double v0 = dom.vmin();
    const double v1 = dom.vmax();
    for (i=0; i<2; i++)
      {
	double x, y, s;
	if (dbg)
	  printf("%f %f     %f %f\n", u0*(1.0-sj0[i]/double(dn)) + u1*sj0[i]/double(dn),
		 v0*(1.0-si0[i]/double(dm)) + v1*si0[i]/double(dm),
		 u0*(1.0-sj1[i]/double(dn)) + u1*sj1[i]/double(dn),
		 v0*(1.0-si1[i]/double(dm)) + v1*si1[i]/double(dm));
	bool inters_found = segment_contour_intersection_for_s2m(u0*(1.0-sj0[i]/double(dn)) + u1*sj0[i]/double(dn),
								 v0*(1.0-si0[i]/double(dm)) + v1*si0[i]/double(dm),
								 u0*(1.0-sj1[i]/double(dn)) + u1*sj1[i]/double(dn),
								 v0*(1.0-si1[i]/double(dm)) + v1*si1[i]/double(dm),
								 &trim_curve_p[0][0], contour, 
								 x, y, s,
								 false // 100222: New 'snap_ends' flag.
								 DBG_FLAG);
	if (dbg)
	  printf("=1===> i=%d, inters=%d\n", i, inters_found?1:0);

	if ((!inters_found) && (with_boundary))
	  {
	    // 090117: Didn't find intersection, even though one should be here. Checking if the segment's ends
	    //         are on the contour.
	    if (is_on_contour(trim_curve_p, contour,
			      u0*(1.0-sj0[i]/double(dn)) + u1*sj0[i]/double(dn),
			      v0*(1.0-si0[i]/double(dm)) + v1*si0[i]/double(dm)))
	      {
		x=u0*(1.0-sj0[i]/double(dn)) + u1*sj0[i]/double(dn);
		y=v0*(1.0-si0[i]/double(dm)) + v1*si0[i]/double(dm);
		s=0.0;
		inters_found=true;
		if (dbg)
		  printf("=2===> i=%d, inters=%d\n", i, inters_found?1:0);
	      }
	    else
	      if (is_on_contour(trim_curve_p, contour,
				u0*(1.0-sj1[i]/double(dn)) + u1*sj1[i]/double(dn),
				v0*(1.0-si1[i]/double(dm)) + v1*si1[i]/double(dm)))
		{
		  x=u0*(1.0-sj1[i]/double(dn)) + u1*sj1[i]/double(dn);
		  y=v0*(1.0-si1[i]/double(dm)) + v1*si1[i]/double(dm);
		  s=1.0;
		  inters_found=true;
		  if (dbg)
		    printf("=3===> i=%d, inters=%d\n", i, inters_found?1:0);
		}
	  }

	if (!inters_found)
	  {
	    MESSAGE("Failed finding intersection. Expected one.");
	    fflush(stderr);
	  } 
	else 
	  {
	    const double eps=1e-13; // 090115: Used for zero-tests for distances in the parameter domain.
	    //         Hmm... these should *really*, *really* be taken from some global
	    //         variable or something
	    // 090117: Anyway... s will never be out of range here...
	    ASSERT2((s>=-eps) && (s<=1.0+eps), printf("Huh?! s=%f\n", s));
	  }
	// eval surf i x og y for aa faa vert og norm, evt. returnere param
	// for skjaering mellom punktene s.a. vi igjen kan danne lin komb.
	if (inters_found)
	  {
	    vert.push_back((1.0-s)*vert[si0[i]*(dn+1)+sj0[i]] +
			   s*vert[si1[i]*(dn+1)+sj1[i]]   );
	    vert_p.push_back((1.0-s)*vert_p[si0[i]*(dn+1)+sj0[i]] +
			     s*vert_p[si1[i]*(dn+1)+sj1[i]]   );
	    if (dbg)
	      printf("=4===> inters: %f %f\n", vert_p[vert_p.size()-1][0], vert_p[vert_p.size()-1][1]);
	    bd.push_back(1);
	    Vector3D new_norm = (1.0-s)*norm[si0[i]*(dn+1)+sj0[i]] +
	      s*norm[si1[i]*(dn+1)+sj1[i]];
	    new_norm.normalize();
	    norm.push_back(new_norm);
	  }
	//extra_v.push_back(vert[vert.size()-1]); // for debugging
      }

    //
    // Add triangles to the mesh.
    //
    // Adding three triangles: (A, B, 0), (B, 1, 0), (0, 2, A).
    //
    //   (i1, j1)  B
    //      o------x----o
    //      |        \  |
    //      |         \ x A
    //      |           |
    //      o-----------o
    //   (i0, j0)    (i2, j2)
    //
    //
    // Adding two triangles: 
    //
    //  "-1"       -----o "-2"		Note that we want the triangles
    //      o-----/     |			oriented ccw., i.e.,
    //      |           |			(-2, -1, 0), and (-2, 0, 1).
    //      |           |
    //      o-----------o
    //   (i0, j0)    (i1, j1)
    //

    // This triangle can be used in all three modes.
    add_triangle(vert_p, mesh, (int)vert.size()-2, (int)vert.size()-1, i0*(dn+1)+j0);
    //col.push_back(Vector3D(1.0, 0.6, 0.6));

    if ((i1==-1) || (j1==-1))
      // Only one triangle to add.
      return;

    if ((i2==-1) || (j2==-1))
      {
	// Two triangles to add, i.e., one more.
	  add_triangle(vert_p, mesh, (int)vert.size()-2, i0*(dn+1)+j0, i1*(dn+1)+j1);
	//col.push_back(Vector3D(0.0, 0.0, 1.0));
      }
    else
      {
	// Three triangles to add, i.e., two more.
	  add_triangle(vert_p, mesh, i0*(dn+1)+j0, (int)vert.size()-1, i1*(dn+1)+j1);
	//col.push_back(Vector3D(0.6, 1.0, 0.6));
      
	  add_triangle(vert_p, mesh, (int)vert.size()-2, i0*(dn+1)+j0, i2*(dn+1)+j2);
	//col.push_back(Vector3D(0.6, 0.6, 1.0));
      }
  }






  //==============================================================================================================
  //
  // 081208: Splitting triangles. The original split_quad should probably have been implemented by two calls to
  //         such a routine as this, it would have made for smaller and more compact code.
  //
  //         This subroutine assumes that either one or two points are inside, and the rest outside, so that
  //         exactly two intersections (between triangle edges and trim curve) are to be found, and either one
  //         or two triangles are to be added.
  //
  //         (If there are 2n, for n>1, edge-curve intersections, the triangle should be split recursively, or
  //         the a priori refinement was too coarse...)
  //
  //         100218: The first corner, c1, is assumed to be inside, and we may then assume that one of three
  //                 possible situations has occured: c2 inside and c3 outside, or c2 outside and c3 inside,
  //                 or both c2 and c3 outside. (If all are inside, we shouldn't be here in the first place.)
  //
  //                 Since at least one corner is inside and one is outside, per assumption, we may also
  //                 assume now (this is made sure by the caller, who may re-order corners) that c1 is inside
  //                 and c2 is outside. This leaves just two cases:
  // 
  //                   1) c1 inside, c2 outside and c3 inside, or
  //                   2) c1 inside, c2 outside and c3 outside.
  // 
  //                 There should then be exactly two intersections, one on 1-2, and one on either 2-3 or on
  //                 1-3.
  //
  // 081208: Note that 'srf' is not evaluated to create new vertices, the old ones are simply interpolated! This
  //         could be improved!
  //
  // 081208: Hmm... Can this interpolation result in non-normalized normals?! Think so! But is it specified that
  //         they should be normalized?! Safest thing would probably be to normalize them... Not doing it for
  //         split_quad, but for split_tirangle...
  //
  // 090129: Fixed bug; the routine did not update 'vert_p', hence buggy later splitting of the newly produced
  //         triangles.
  //
  // 090203: Hmm... There is a problem here, when the trimming curve passes through one or more of the corners
  //         of the triangle. How can this be solved?
  //         First idea: Simply disallow the second intersection to be the same as the first.
  //                     One advantage of this, if possible, is that 'segment_contour_intersection_for_s2m' will
  //                     not have to be modified. Such a modification could affect everything else in a somewhat
  //                     unpredictable manner...
  //
  // 090217: Adding 'forced_skipping_of_second_edge' so that we get more control. This makes sure some triangles
  //         with more than two intersections are correctly (i.e., better) handled. Unfortunately, it doesn't
  //         solve all problems...
  //
  // 090218: Return value introduced. Will be set to true if further processing of produced triangles is needed.
  //
  // 100218: There is a problem with cases where the contour intersects (or touches) the triangles in more
  //         than two places. This may for instance happen when the contour has a "thin corner" like this:
  //
  //                                         \          |
  //                                   c2 +---\---------o c3  (c3=) i2 (2-3: s2=1)
  //                                       \   \       /|
  //                                         \  \     / |
  //                                           \ \   /  |
  //                                             \\ /   |
  //                                               X    |
  //                                             c1 \   |
  //                            (c1=) i1 (1-2: s1=0) \  |
  //                                                  \ |
  //                                                   \|
  //
  //         Here, c2 is outside, c1 and c3 inside (really, "on", but that info is not passed on to
  //         'split_triangle' as of today) and we get intersections i1 and i2 on the edges 1-2 and 2-3,
  //         respectively. This will result in the triangles (c1, i1, i2) and (c1, i2, c3), but 
  //
  // 100222: Another problem: Previously to getting here, a point may have been classified as "inside" due to
  //         in fact being "on" the contour, and when we later get here for splitting, the point was suddenly
  //         not being classified as "on" anymore because the contour segment which it is "on" is parallel to
  //         the ray being used in the "crossings loop" in the inside-test. To fix this, we must (can) snap
  //         points to ends of (such) segments. This is now done by sending the appropriate flag ('snap_ends')
  //         to 'segment_contour_intersection_for_s2m'. Currently, this routine is the only one setting that
  //         flag to true.
  //
  // 100223: Adding a flag for "innerness" of contours, because this is needed in order to take the
  //         appropriate action when the first intersection ('first_inters') is not found, see below.
  //
  // 100223: New situation not anticipated arises for 'bootm_face.g2'... c1 inside, c2 and c3 outside, none
  //         "on", but there are four intersections on the edges, and nothing really goes wrong, we just get
  //         an incorrect result due to this unforeseen and complex situation:
  //                                                                                                   
  //                                                                                                   
  //                                                 |          1    /contour                                 
  //                                                 |         /|   /                                  
  //                                                 |       /  |  /                                   
  //                                                 |     /    | /                                    
  //                                                 |   /  T   |/                                     
  //                                              i4 | /....... / i1                                     
  //                                                 |         /|                                      
  //                                               / |   Q    / |                                      
  //                                             /   |       /  |                                       
  //                                           3 ----|------/-- 2
  //                                              i3 |     / i2                                         
  //                                                 |    /                                           
  //                                                 |   /                                            
  //                                                 |  /                                             
  //                                                 | /                                              
  //                                                 |/                                               
  //
  //         The result will be the triangle T, and we miss the quad Q, which could have been composed of two
  //         triangles. This is just one possible hard configuration. An implicit assumption of the current
  //         approach is that the contours cannot intersect triangles in more than at most two places. Of
  //         course, a real life contour might intersect a triangle any number of times...
  //
  //         I have a feeling that 'forced_skipping_of_second_edge' possibly should be removed, and maybe
  //         replaced with some recursive refinement variation for cases where more than the expected maximal
  //         two intersections can be detected... Adding an experimental switch for this.
  //
  // 100223: In order to try not to break too much of the old code, and to improve readability, we'll try to
  //         resolve the situation with another processing stage: Before entering the main processing of the
  //         "old" 'split_triangle', we will check if we can find three distinct intersections of the edges of
  //         the triangle. If so, we split the triangle first, then call the "old" routine for each of these
  //         triangles. The simplest way to do this is to just split and return 'true' (meaning that further
  //         processing must be done) and then let 'trim_a_triangle' take care of the rest!
  //
  // 100224: It is really, really time to get rid of 'forced_skipping_of_second_edge'. In a very few cases, it
  //         may still produce better results, even though most of the time it is now (with the new
  //         three-intersecting-edges-preprocessing stage) obsolete... If those few cases could be solved by
  //         other means, it would make for a nice cleanup in the removal of
  //         'forced_skipping_of_second_edge'...
  //
  // 100224: As an attempt to fix the above mentioned situation, and possibly many similar ones, a new
  //         fallback for triangles of case 's==0' in 'trim_a_triangle' will be added. It is also an advantage
  //         in not cluttering up the current 'split_triangle', which is complex enough as it is.
  //
  //==============================================================================================================




  //==============================================================================================================
  //
  // 100223: For increased readability...
  //         Not using references, in fear of what push_backing on the vertex lists will do... 
  //         (Hmm, should be safe, I think...)
  //
  // 100224: Hmm... It is probably not necessary to do this linear interpolation with 's', as
  //         'segment_contour_intersection_for_s2m' seems to return the same result in its 'x' and 'y'
  //         parameters... Remember to fix this after the stuff works as it should...
  //
  // 100225: Now evaluating the surface instead of interpolating triangle corners in 3D.
  //
  //==============================================================================================================

  inline void push_an_intersection(const shared_ptr<ParamSurface> srf,
				   const Vector3D &c1, const Vector3D &c2,
				   const Vector2D &c1_p, const Vector2D &c2_p,
				   const Vector3D &n1, const Vector3D &n2,
				   const double &s,
				   vector<Vector3D> &vert, vector<Vector2D> &vert_p,
				   vector<Vector3D> &norm,
				   vector<int> &bd)
  {
    int dim = srf->dimension();
    vert_p.push_back(  (1.0-s)*c1_p + s*c2_p );

#ifndef EVAL_SRF_NOT_INTERP
    vert  .push_back(  (1.0-s)*c1   + s*c2   );
    Vector3D tmp =     (1.0-s)*n1   + s*n2;
    tmp.normalize();
    norm.push_back(tmp);
    bd.push_back(1);
#else
    vector<Point> res(3);
    srf->point(res, vert_p.back()[0], vert_p.back()[1], 1); // Could we have used 0 here? 

    // The point on the surface:
    if (dim == 3)
      vert.push_back( Vector3D(res[0].begin()) );
    else
      vert.push_back(Vector3D(res[0][0], res[0][1], 0.0));

    // The normal in the same point:
    Point nrm;
    if (dim == 3)
      nrm = res[1].cross(res[2]);
    else 
      nrm = Point(0.0, 0.0, 1.0);
    Vector3D tmp = Vector3D(nrm);
    tmp.normalize();
    norm.push_back(tmp);
    
    // This is Vibeke's, don't know what it's for...
    bd.push_back(1);
#endif

  }

    


  bool split_triangle(const int c1_indx, const int c2_indx, const int c3_indx,
		      shared_ptr<ParamSurface> srf,
		      vector< Vector3D > &vert, vector<Vector2D> &vert_p,
		      vector< int > &bd, vector< Vector3D > &norm,
		      vector<int> &newmesh, 
		      const vector< Vector3D > &trim_curve_p,
		      const vector<int> &contour,
		      const bool dbg = false,
		      /* const */ bool forced_skipping_of_second_edge = false,
		      const bool inner_trimming_curve = false) // 100223: See above.
  {
    const double abs_eps = 1e-12;	// 100218: Used for distance between two intersections in the param domain.

    // 100223: Testing a new approach. Want to get rid of 'forced_skipping_of_second_edge', but trying to do
    //         it one step at a time, through these temporary switches...

    const bool try_without_forced_stuff = true;
    
    if (try_without_forced_stuff)
      forced_skipping_of_second_edge = false;

    // 090129: Not using references in fear of what 'push_back' will do...
    const Vector2D c1_p=vert_p[c1_indx], c2_p=vert_p[c2_indx], c3_p=vert_p[c3_indx];
    const Vector3D c1=vert[c1_indx], c2=vert[c2_indx], c3=vert[c3_indx];
    const Vector3D n1=norm[c1_indx], n2=norm[c2_indx], n3=norm[c3_indx];

    if (dbg)
      {
	puts("\n\n  ----- split_triangle starting -----------------------------------------------------------------");
	printf("\n    The corners:\n\n");
	printf("hold on; plot(%f, %f, 'rd', 'markersize', %d, 'linewidth', 2); hold off\n", c1_p[0], c1_p[1], 14);
	printf("hold on; plot(%f, %f, 'gd', 'markersize', %d, 'linewidth', 2); hold off\n", c2_p[0], c2_p[1], 14);
	printf("hold on; plot(%f, %f, 'bd', 'markersize', %d, 'linewidth', 2); hold off\n\n", c3_p[0], c3_p[1], 14);
      }

    double x1, y1, s1, x2, y2, s2;
    
    if (dbg) printf("    Looking for intersection on edge 1-2...\n");
    const bool first_inters =
      segment_contour_intersection_for_s2m(c1_p[0], c1_p[1], c2_p[0], c2_p[1], 
					   &trim_curve_p[0][0], // Pointer to vertices is enough, size not needed
					   contour, x1, y1, s1, 
					   true 		// 100222: New 'snap_ends' flag.
					   DBG_FLAG);
    if (!first_inters)
      {
	if (!inner_trimming_curve)
	  {
	    // 100223: Disabling this warning, for with the new branch depending on "innerness" of contour,
	    //         the situation will probably (hopefully) be correctly treated in all cases now.
	    if ( (dbg) && (0) )
	      puts("    Warning:\n  It is likely that the 'inside'-test gave a false positive for a corner of a\n"
		   "    triangle that is really fully outside a trimming curve. We continue with this\n"
		   "    assumption, and discard the triangle in question, without attempting to clip it.\n");
	  }
	else
	  // 100223: The new and alternate branch, which should be more likely to be correct for inner curves.
	  add_triangle(vert_p, newmesh, c1_indx, c2_indx, c3_indx);
	if (dbg) puts("  ----- split_triangle done, no new triangles produced ---------------------------------\n\n");
	return false;
      }
  
    //
    // The first intersection is now between corner c1 and c2. The second may be between c2 and c3, or between
    // c3 and c1. Per assumption, c1 is *inside* the contour. (And c2 outside, c3 unknown.)
    //
    // 100223: But remember, now, that there may in fact be a lot more than the anticipated/assumed two
    //         intersections... trying to cope better with such situations from now on.
    //
    
    // Pushing inters1.
    push_an_intersection(srf, c1, c2, c1_p, c2_p, n1, n2, s1, vert, vert_p, norm, bd);
    if (dbg)
      printf("    First intersection found: s1=%f\n\n      hold on; plot(%f, %f, 'm*', 'markersize', "
	     "7, 'linewidth', 2); hold off\n\n", s1, vert_p.back()[0], vert_p.back()[1]);
    
    // 100218: Note that we only have to *test* for a 2-3 intersection, since if there is none, we should have a 3-1.
    if (dbg) printf("    Looking for intersection on edge 2-3...\n");
    const bool second_inters =
      segment_contour_intersection_for_s2m(c2_p[0], c2_p[1], c3_p[0], c3_p[1], &trim_curve_p[0][0], contour, x2, y2, s2,
					   true // 100222: New 'snap_ends'-flag
					   DBG_FLAG);

    const double i1_i2_dist_squared = second_inters ? (x2-x1)*(x2-x1) + (y2-y1)*(y2-y1) : 0.0;
    if (dbg) printf("    second_inters (2-3) =%d, s2=%f, dist=%g\n", second_inters, s2, sqrt(i1_i2_dist_squared));
    const bool inters1_and_inters2_distinct = i1_i2_dist_squared > abs_eps*abs_eps;


    //--------------------------------------------------------------------------------------------------------------
    //
    // 100223: When three intersections between triangle edges and the contour are detected, the triangle is
    //         split in four, and flagged for further trimming in a new iteration. This should take care of a
    //         lot of previously unsolved cases with minimal new code!
    //
    //--------------------------------------------------------------------------------------------------------------

    bool third_inters;
    const bool look_for_a_third_inters = first_inters && second_inters && inters1_and_inters2_distinct;
    if (dbg) printf("    look_for_a_third_inters = %d\n", look_for_a_third_inters);
    double x3, y3, s3;
    if ( look_for_a_third_inters )
      {
	if (dbg) printf("    PRE_SPLITTING_3_INTERS: Looking for intersection on edge 3-1...\n");
	third_inters = segment_contour_intersection_for_s2m(c3_p[0], c3_p[1], c1_p[0], c1_p[1],
							    &trim_curve_p[0][0], contour, x3, y3, s3,
							    true // 100222: New 'snap_ends'-flag
							    DBG_FLAG);
	if (third_inters)
	  {
	    const double i1_i3_dist_squared = (x3-x1)*(x3-x1) + (y3-y1)*(y3-y1);
	    const bool inters1_and_inters3_distinct = i1_i3_dist_squared > abs_eps*abs_eps;
	    const double i2_i3_dist_squared = (x3-x2)*(x3-x2) + (y3-y2)*(y3-y2);
	    const bool inters2_and_inters3_distinct = i2_i3_dist_squared > abs_eps*abs_eps;
	    
	    if (inters1_and_inters3_distinct && inters2_and_inters3_distinct)
	      {
		// puts("JEPP, FANT ET KJIPT LITE TRIANGEL...");
		
		// Pushing inters2.
		push_an_intersection(srf, c2, c3, c2_p, c3_p, n2, n3, s2, vert, vert_p, norm, bd);
		// Pushing inters3.
		push_an_intersection(srf, c3, c1, c3_p, c1_p, n3, n1, s3, vert, vert_p, norm, bd);
		
		// Triangles...
		add_triangle(vert_p, newmesh, c1_indx, (int)vert.size()-3, (int)vert.size()-1 DBG_FLAG);
		add_triangle(vert_p, newmesh, c2_indx, (int)vert.size()-2, (int)vert.size()-3 DBG_FLAG);
		add_triangle(vert_p, newmesh, c3_indx, (int)vert.size()-1, (int)vert.size()-2 DBG_FLAG);
		add_triangle(vert_p, newmesh, (int)vert.size()-3, (int)vert.size()-2, (int)vert.size()-1 DBG_FLAG);
		
		return true;
	      }
	  }
      }

    if ( second_inters && inters1_and_inters2_distinct && (!forced_skipping_of_second_edge) )
      {
	//----------------------------------------------------------------------------------------------------
	//
	// The second intersection is between corners c2 and c3, and it is distinct from the first
	// intersection. We can try to make the triangles as round as possible...
	//
	// 100223: What happens when 'forced_skipping_of_second_edge=true', is that this branch is not taken,
	//         although a perfectly fine intersection distinct from inters1 was found on 2-3. That will
	//         cause an intersection on 3-1 to be searched for, and used instead. (And who knows what will
	//         happen if that 3-1 intersection does not exist...) If we were to enter this branch, it
	//         seems a lot of triangles will be wrongly discarded, for some reasone. Cannot actually see
	//         how that will happen in this branch, right now, though... Could be that one of the four
	//         'add_triangle' calls below try to add a degenerate triangle. Hmm... seems not so unlikely.
	//
	//----------------------------------------------------------------------------------------------------

	// Pushing inters2.
	push_an_intersection(srf, c2, c3, c2_p, c3_p, n2, n3, s2, vert, vert_p, norm, bd);
	if (dbg)
	  printf("    branch 1:\n\n    hold on; plot(%f, %f, 'g*', 'markersize', 7, 'linewidth', 2); hold off\n\n", 
		 vert_p.back()[0], vert_p.back()[1]);
      
	if ( (x2-c1_p[0])*(x2-c1_p[0])+(y2-c1_p[1])*(y2-c1_p[1]) < (x1-c3_p[0])*(x1-c3_p[0])+(y1-c3_p[1])*(y1-c3_p[1]) )
	  // Splitting between c1 and inters2.
	  {
	    if (dbg) printf("    branch 2: Adding triangles (c1, i1, i2) and (c1, i2, c3).\n");
	    if (dbg) printf("      adding the first:\n");
	    add_triangle(vert_p, newmesh, c1_indx, (int)vert.size()-2, (int)vert.size()-1 DBG_FLAG);
	    if (dbg) printf("      adding the second:\n");
	    add_triangle(vert_p, newmesh, c1_indx, (int)vert.size()-1, c3_indx DBG_FLAG);
	  }
	else
	  // Splitting between c3 and inters1.
	  {
	    if (dbg) printf("    branch 3: Adding triangles (c1, i1, c3) and (i1, i2, c3).\n");
	    add_triangle(vert_p, newmesh, c1_indx, (int)vert.size()-2, c3_indx DBG_FLAG);
	    add_triangle(vert_p, newmesh, (int)vert.size()-2, (int)vert.size()-1, c3_indx DBG_FLAG);
	  }
      }
    else
      {
	//----------------------------------------------------------------------------------------------------
	// 
	// Ok, no second intersection 2-3, then it should be on 3-1. If not, that's an error... (There should
	// always be exactly two intersections!) (100223: Yeah... well... not entirely true, really...)
	//
	//----------------------------------------------------------------------------------------------------

	if (dbg) printf("    branch 4:\n    Looking for intersection on edge 3-1...\n");
	bool tmp;
	if (look_for_a_third_inters) // If so, we have already called 'segment_contour_intersection_for_s2m'...
	  tmp = third_inters, x2 = x3, y2 = y3, s2 = s3;
	else
	  tmp = segment_contour_intersection_for_s2m(c3_p[0], c3_p[1], c1_p[0], c1_p[1], 
						     &trim_curve_p[0][0], contour, x2, y2, s2,
						     true); // 100222: New 'snap_ends'-flag

	// 100224: Hmm... tmp only used for debugging output? (Hence the name 'tmp'?!)

	if (dbg)
	  printf("  2nd try: 3-1-inters=%d, dist=%g\n  forced_skipping_of_second_edge=%d\n\n"
		 "    hold on; plot(%f, %f, 'g*', 'markersize', 8, 'linewidth', 2); hold off\n\n",
		 tmp?1:0, sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1)), forced_skipping_of_second_edge, x2, y2);

	//
	// 100224: Ok, before, we had 'if (!forced...) { ... }' here, because it would be detected as a
	//         problem if we got here, i.e., one 1-2-intersection and no 2-3-intersection *without* the
	//         caller expecting this and therefore using 'forced...'.
	//
	//         If we now want to get rid of the 'forced...' flag, we cannot trap this as a potential
	//         problem, since we have no way of knowing if the caller expected it. If we get here, we must
	//         simply assume that it is correct, and continue with the splitting... Therefore: We wrap the
	//         test in the new 'do_ignore_forced_skipping'-test:
	//

	if (!try_without_forced_stuff)
	  if (!forced_skipping_of_second_edge)
	    {
	      MESSAGE("  Warning:\nShould not happen 2 ! Where did the second intersection go?!"
		      //"\nJust keeping the full triangle...\n");
		      "\n  Discarding the full triangle...\n");
	      fflush(stdout);
	      fflush(stderr);
	      //add_triangle(vert_p, newmesh, c1_indx, c2_indx, c3_indx);
	      
	      if (dbg)
		{
		  if (degenerate_triangle(vert_p[c1_indx], vert_p[c2_indx], vert_p[c3_indx] DBG_FLAG))
		    puts("Hmm... triangle is degenerate, so maybe appropriate to discard it then...");
		  int m=5*(c1_indx%4+1);
		  printf("\n\nwwwwwwwwwwww\n\n\n"
			 "hold on; plot(%f, %f, 'rs', 'markersize', %d, 'linewidth', 2); hold off\n"
			 "hold on; plot(%f, %f, 'gs', 'markersize', %d, 'linewidth', 2); hold off\n"
			 "hold on; plot(%f, %f, 'bs', 'markersize', %d, 'linewidth', 2); hold off\n\n",
			 c1_p[0], c1_p[1], m, c2_p[0], c2_p[1], m, c3_p[0], c3_p[1], m);
		  puts("  ----- split_triangle done ------------------------------------------------------------\n\n");
		}
	      return false;
	    }
	
	//==============================================================================================================
	// 090217: Now, if the first intersection is through the first corner, and there is another "non-corner"
	//         intersection in the "middle" edge (and 'forced_skipping_of_second_edge' is true) we do some
	//         special treatment: We split the triangle in two new ones covering the old. Note that these
	//         must then be recursively re-processed by the caller! To make this simpler, we introduce a
	//         boolean return value, which will be true when further processing may be required.
      
	if ( ((fabs(s1)<1e-12) || (fabs(s2-1.0)<1e-12)) && forced_skipping_of_second_edge )
	  {
	    double x3, y3, s3;
	    if (dbg) printf("    Looking for intersection on edge 2-3 again...\n");
	    const bool tmp2=
	      segment_contour_intersection_for_s2m(c2_p[0], c2_p[1], c3_p[0], c3_p[1], 
						   &trim_curve_p[0][0], contour, x3, y3, s3,
						   true); // 100222: New 'snap_ends' flag.
	    if ( tmp2 && (fabs(s3)>1e-12) && (fabs(s3-1.0)<1-1e-12) )
	      {
		// Pushing inters3.
		push_an_intersection(srf, c2, c3, c2_p, c3_p, n2, n3, s3, vert, vert_p, norm, bd);
		// Hmm... now we have added a totally unneeded point, for intersection 1...
		add_triangle(vert_p, newmesh, c1_indx, (int)vert.size()-1, c3_indx);
		add_triangle(vert_p, newmesh, c1_indx, c2_indx, (int)vert.size()-1);
	      
		if (dbg) puts("  ----- split_triangle done -----------------------------------------------------\n\n");
		return true;
	      }
	  }
      
	// Pushing inters2.
	push_an_intersection(srf, c3, c1, c3_p, c1_p, n3, n1, s2, vert, vert_p, norm, bd);
	if (dbg)
	  printf("  s2=%g\n  hold on; plot(%f, %f, 'b*', 'markersize', 7, 'linewidth', 2); hold off\n",
		 s2, vert_p[vert_p.size()-1][0], vert_p[vert_p.size()-1][1]);
      
	// There is only one triangle to push on the list in this case.
	if (dbg) printf("    adding the only triangle...\n");
	add_triangle(vert_p, newmesh, c1_indx, (int)vert.size()-2, (int)vert.size()-1);

      } // end of the block for the second intersection being of 3-1 kind instead of 2-3...

    if (dbg) puts("  ----- split_triangle done ----------- xxx ------------------------------------------------\n\n");
    return false;
  }






#ifdef SPLIT_LARGE_TRIANGLES

  //==============================================================================================================
  //
  // 100224: Making a version for splitting a triangle which would otherwise be discarded, due to all corners
  //         being classified as outside the contour. Such a triangle might still contain intersections.
  //
  //         Hmm... there is a danger of inifinite recursion here, how do we handle this?
  //
  // 100225: Testing for the need of splitting in this function too. (Note that in the
  //         'trim_a_triangle'/'split_triangle'-case this is divided between the two functions.)
  //
  //         New name: 'split_triangle_with_all_corners_outside', better reflecting the functionality.
  //
  //         Assuming all corners are outside. Not assuming anything about their "on-ness". (100225: Returning
  //         of none are "on".)
  // 
  //         Adding functionality for new cases as they appear. First case:
  //
  //           A corner is "on", and there is an *inner* intersection on the opposite edge. Then we split the
  //           triangle in two here. Hopefully, this approach is so conservative that nothing will be broken
  //           of those things already working...
  //
  //
  //        Return value: False if no further processing of produced triangles are needed. Else true.
  //
  //==============================================================================================================

  bool split_triangle_with_all_corners_outside(const int c1_indx, const int c2_indx, const int c3_indx,
					       const bool c1_inside, const bool c2_inside, const bool c3_inside,
					       const bool c1_on, const bool c2_on, const bool c3_on,
					       shared_ptr<ParamSurface> srf,
					       vector< Vector3D > &vert, vector<Vector2D> &vert_p,
					       vector< int > &bd, vector< Vector3D > &norm,
					       vector<int> &newmesh, 
					       const vector< Vector3D > &trim_curve_p,
					       const vector<int> &contour,
					       const bool dbg = false,
					       const bool inner_trimming_curve = false)
  {
    const double abs_eps = 1e-12;	// 100218: Used for distance between two intersections in the param domain.
    const double snap_eps = 1e-5;

    // 090129: Not using references in fear of what 'push_back' will do...
    const Vector2D c1_p=vert_p[c1_indx], c2_p=vert_p[c2_indx], c3_p=vert_p[c3_indx];
    const Vector3D c1=vert[c1_indx], c2=vert[c2_indx], c3=vert[c3_indx];
    const Vector3D n1=norm[c1_indx], n2=norm[c2_indx], n3=norm[c3_indx];

    if (dbg)
      {
	puts("\n\n  ----- split_triangle_with_all_corners_outside ------------------------------------------------");
	printf("\n    The corners:\n\n");
	printf("hold on; plot(%f, %f, 'rd', 'markersize', %d, 'linewidth', 2); hold off\n", c1_p[0], c1_p[1], 14);
	printf("hold on; plot(%f, %f, 'gd', 'markersize', %d, 'linewidth', 2); hold off\n", c2_p[0], c2_p[1], 14);
	printf("hold on; plot(%f, %f, 'bd', 'markersize', %d, 'linewidth', 2); hold off\n\n", c3_p[0], c3_p[1], 14);
      }
    
    if (!(c1_on || c2_on || c3_on))
      // 100225: At this time, we don't do anything with such a triangle.
      return false;
    
    // 100225: Ok, at least one corner is "on", so let's make sure c1 is "on"
    if (!c1_on)
      {
	if (c2_on)
	  return split_triangle_with_all_corners_outside(c2_indx, c3_indx, c1_indx,
							 c2_inside, c3_inside, c1_inside,
							 c2_on, c3_on, c1_on,
							 srf, vert, vert_p, bd, norm, newmesh, trim_curve_p,
							 contour, dbg, inner_trimming_curve);
	else
	  return split_triangle_with_all_corners_outside(c3_indx, c1_indx, c2_indx,
							 c3_inside, c1_inside, c2_inside,
							 c3_on, c1_on, c2_on,
							 srf, vert, vert_p, bd, norm, newmesh, trim_curve_p,
							 contour, dbg, inner_trimming_curve);
      }
    
    // 100225: Now c1 is "on". Checking for an intersection (inner) on edge 2-3.

    double x1, y1, s1;
    if (dbg) printf("    Looking for intersection on edge 2-3...\n");
    const bool inters =
      segment_contour_intersection_for_s2m(c2_p[0], c2_p[1], c3_p[0], c3_p[1], 
					   &trim_curve_p[0][0], // Pointer to vertices is enough, size not needed
					   contour, x1, y1, s1, 
					   true 		// 100222: New 'snap_ends' flag.
					   DBG_FLAG);
    if (!inters)
      return false;

    const double dist2_squared = (x1-c2_p[0])*(x1-c2_p[0]) + (y1-c2_p[1])*(y1-c2_p[1]);
    const double dist3_squared = (x1-c3_p[0])*(x1-c3_p[0]) + (y1-c3_p[1])*(y1-c3_p[1]);
    const bool interior = (dist2_squared>snap_eps*snap_eps) && (dist3_squared>snap_eps*snap_eps);
    if (dbg) printf("    dist2=%f, dist3=%f, interior=%d\n", sqrt(dist2_squared), sqrt(dist3_squared), interior);
    if (!interior)
      return false;
    
    // 100225: Ok, all conditions are met, we split.

    push_an_intersection(srf, c2, c3, c2_p, c3_p, n2, n3, s1, vert, vert_p, norm, bd);
    add_triangle(vert_p, newmesh, c1_indx, c2_indx, vert.size()-1 DBG_FLAG);
    add_triangle(vert_p, newmesh, c1_indx, vert.size()-1, c2_indx DBG_FLAG);
    if (dbg) puts("  ----- split_triangle_with_all_corners_outside, produced two new triangles ----------------\n\n");
    return true;
  }

#endif






  //==============================================================================================================
  //
  // 090202: Hmm... Should maybe introduce some tolerances here too...
  //
  // 090202: Note that this function is supposed to handle 3D points, even though we in this case have 2D points
  //         in the parameter domain. (In the 3D case, the point tested will be the projection onto the plane
  //         containing the triangle, if I remember correctly...) I think just replacing Vector3D with Vector2D
  //         actually will do the trick...
  //
  //==============================================================================================================

  inline bool pt_inside_tri(const Vector3D &p, const Vector3D &c1, const Vector3D &c2, const Vector3D &c3)
  {
    const Vector3D e=c2-c1, f=c3-c1, pc=p-c1;
    const double ee=e*e, ff=f*f, ef=e*f;
#if 0
    // This is faster on x86 than allocating new variables for pc*e and pc*f, it seems...
    const double u = ff*(pc*e) - ef*(pc*f);
    const double v = ee*(pc*f) - ef*(pc*e);
#else
    // Not on core2 is seems...
    const double pce=pc*e, pcf=pc*f;
    const double u = ff*pce - ef*pcf;
    const double v = ee*pcf - ef*pce;
#endif

    return ((u>=0.0) && (v>=0.0) && (u+v<=ee*ff-ef*ef));
  }

  inline bool pt_inside_tri(const Vector2D &p, const Vector2D &c1, const Vector2D &c2, const Vector2D &c3)
  {
    const Vector2D e=c2-c1, f=c3-c1, pc=p-c1;
    const double ee=e*e, ff=f*f, ef=e*f;
#if 0
    // This is faster on x86 than allocating new variables for pc*e and pc*f, it seems...
    const double u = ff*(pc*e) - ef*(pc*f);
    const double v = ee*(pc*f) - ef*(pc*e);
#else
    // Not on core2 is seems...
    const double pce=pc*e, pcf=pc*f;
    const double u = ff*pce - ef*pcf;
    const double v = ee*pcf - ef*pce;
#endif

    return ((u>=0.0) && (v>=0.0) && (u+v<=ee*ff-ef*ef));
  }

  //==============================================================================================================
  //
  // 090203: To avoid splitting triangles into almost degenerate ones.
  //
  //         Assuming 'pt_inside_tri' is true, so that we can measure distance to infinite extensions of edges
  //         for simplicity.
  //
  //==============================================================================================================

  inline double pt_dist_from_line(const Vector2D &p, const Vector2D &A, const Vector2D &B)
  {
    const double t = ((p-A)*(B-A))/((B-A)*(B-A));
    const Vector2D proj = A+t*(B-A);
    return sqrt((proj-p)*(proj-p));
  }

  inline double pt_dist_from_tri_edge(const Vector2D &p, const Vector2D &c1, const Vector2D &c2, const Vector2D &c3)
  {
    return std::min(pt_dist_from_line(p, c1, c2),
		    std::min(pt_dist_from_line(p, c2, c3), pt_dist_from_line(p, c3, c1)));
  }






  //==============================================================================================================
  //
  // 090217: Encapsulating this for clarity.
  //
  //==============================================================================================================

  bool trim_a_triangle(shared_ptr<ParamSurface> srf,
		       vector<Vector3D> &vert,
		       vector<Vector2D> &vert_p, 
		       vector< int > &bd, vector<Vector3D> &norm,
		       const vector<Vector3D> &trim_curve_p,
		       const vector<int> &contour,
		       const int c1_indx, const int c2_indx, const int c3_indx, 
		       vector<int> &mesh,
		       const bool inner_trimming_curve = false
#ifdef DBG
		       , const bool dbg = false
#endif
    )
  {
#ifndef DBG
    const bool dbg=false;
#endif
    if (dbg)
      printf("\n\n= trim_a_triangle starting =====================================================================\n");

    vector<int> new_mesh;

    const Vector3D c1=vert[c1_indx], c2=vert[c2_indx], c3=vert[c3_indx];
    const Vector3D n1=norm[c1_indx], n2=norm[c2_indx], n3=norm[c3_indx];
    const Vector2D c1_p=vert_p[c1_indx], c2_p=vert_p[c2_indx], c3_p=vert_p[c3_indx];
    int s0, s, s3;

    if (dbg)
      {
	printf("\n  The corners:\n\n");
	printf("hold on; plot(%f, %f, 'rx', 'markersize', %d, 'linewidth', 2); hold off\n", c1_p[0], c1_p[1], 14);
	printf("hold on; plot(%f, %f, 'gx', 'markersize', %d, 'linewidth', 2); hold off\n", c2_p[0], c2_p[1], 14);
	printf("hold on; plot(%f, %f, 'bx', 'markersize', %d, 'linewidth', 2); hold off\n\n", c3_p[0], c3_p[1], 14);
      }

    if (dbg) printf("  = Calling is_inside (point_inside_contour) for corner c1... ============================\n");
    const int c1_inside = is_inside(trim_curve_p, contour, c1_p[0], c1_p[1] DBG_FLAG);
    if (dbg) printf("  = Calling is_inside (point_inside_contour) for corner c2... ============================\n");
    const int c2_inside = is_inside(trim_curve_p, contour, c2_p[0], c2_p[1] DBG_FLAG);
    if (dbg) printf("  = Calling is_inside (point_inside_contour) for corner c3... ============================\n");
    const int c3_inside = is_inside(trim_curve_p, contour, c3_p[0], c3_p[1] DBG_FLAG);

    if (dbg) printf("  = Calling is_on_contour (point_on_contour) for corner c1... ============================\n");
    const int c1_on     = is_on_contour(trim_curve_p, contour, c1_p[0], c1_p[1] DBG_FLAG);
    if (dbg) printf("  = Calling is_on_contour (point_on_contour) for corner c2... ============================\n");
    const int c2_on     = is_on_contour(trim_curve_p, contour, c2_p[0], c2_p[1] DBG_FLAG);
    if (dbg) printf("  = Calling is_on_contour (point_on_contour) for corner c3... ============================\n");
    const int c3_on     = is_on_contour(trim_curve_p, contour, c3_p[0], c3_p[1] DBG_FLAG);

    if (dbg)
      printf("----------------------------------------------------------------------------------------------------\n"
	     "IN/ON-tests before considering whether curve is outer or inner:\n"
	     "c1_inside=%d, c2_inside=%d, c3_inside=%d, c1_on=%d, c2_on=%d, c3_on=%d\n"
	     "----------------------------------------------------------------------------------------------------\n",
	     c1_inside, c2_inside, c3_inside, c1_on, c2_on, c3_on);

    // 100211: s2==0 <=> no corner on curve itself.
    const int s2 = (c1_on << 0)  +  (c2_on << 1)  +  (c3_on << 2);

    const Vector2D centroid = (1.0/3.0)*(c1_p+c2_p+c3_p);
    const int c_inside  = is_inside(trim_curve_p, contour, centroid[0], centroid[1]);
    
    if (inner_trimming_curve)
      {
	// 100211: s0==0 <=> all corners are not inside. (Is "on edge" included in "inside"?)
	s0 = ((!c1_inside) << 0) + ((!c2_inside) << 1) + ((!c3_inside) << 2);

	// 100211: 
	s = (((!c1_inside) | c1_on) << 0)  +  (((!c2_inside) | c2_on) << 1)  +  (((!c3_inside) | c3_on) << 2);
	
	// 100211: s3==true <=> "centroid is not inside".
	s3 = !c_inside;

	if (dbg)
// 	  if ((s0!=7) || (s2!=0) || (s!=7))
	  printf("  s0=%d, s2=%d, s=%d, s3=%d\n", s0, s2, s, s3);
	
	if ( (s0==0) && ((s2==1) || (s2==2) || (s2==4)) )
	  {
	    if (dbg)
	      printf("No points on the inside, exactly one *on* the contour. Setting s=0 and skipping triangle. 1\n");
	    s=0;
	  }
	if ( ((s0==1) && (s2==1)) || ((s0==2) && (s2==2)) || ((s0==4) && (s2==4)) )
	  {
	    if (dbg)
	      printf("Two points not on the inside, the remaining both inside and *on*. Skipping triangle. 2\n");
	    s=0;
	  }

	if (dbg)
	  if ( (s==7) && (s3==0) )
	    {
	      printf("All corners on the inside, mean not, skipping triangle.\n");
	      // 090220: This is not so straightforward... Not all these triangles are either fully inside or
	      //         outside. They should really be split, but then we should split neighbours also, if we
	      //         ever are to get a valid triangulation. (Which we are currently not getting, since we
	      //         produce lots of duplicate nodes and so on...)
	      //s=0;
	    }

	if (dbg)
 	  if ( ((s2==3) && (s0==4)) ||
 	       ((s2==5) && (s0==2)) ||
 	       ((s2==6) && (s0==1))    )
 	    {
 	      printf("Two points on the curve, one inside. Skipping triangle.\n");
// 	      s=0;
 	    }
	
	if (dbg)
 	  if (s2==7)
 	    {
 	      printf("All points on contour, guessing that the triangle should be excluded. \n");
 	      printf("  s0=%d, s2=%d, s=%d, s3=%d\n", s0, s2, s, s3);
// 	      s=0;
 	    }

      }
    else
      {
	s0 = (c1_inside << 0) + (c2_inside << 1) + (c3_inside << 2);
	s = ((c1_inside | c1_on) << 0)  +  ((c2_inside | c2_on) << 1)  +  ((c3_inside | c3_on) << 2);
	s3 = c_inside;

	if (dbg)
	  printf("---> s0=%d s2=%d s=%d s3=%d\n", s0, s2, s, s3), fflush(stdout);

	if ( (s0==0) && ((s2==1) || (s2==2) || (s2==4)) )
	  {
	    if (dbg)
	      printf("No points on the inside, exactly one *on* the contour. Setting s=0 and skipping triangle. 4\n");
	    // s=0;
	  }

	if ( ((s0==1) && (s2==1)) || ((s0==2) && (s2==2)) || ((s0==4) && (s2==4)) )
	  {
	    if (dbg)
	      printf("Two points not inside, the remaining both inside and *on*. Skipping triangle. s=%d 5\n", s);
	    //s=0;
	    //s=7;
	  }

	if (dbg)
	  if ((s==7) && (!s3))
	    {
	      printf("====================> HUH?! s0=%d s2=%d s=%d s3=%d\n", s0, s2, s, s3);
// 	  s=0;
	    }
	
	if (dbg)
	  printf("---> s0=%d s2=%d s=%d s3=%d\n", s0, s2, s, s3), fflush(stdout);
	
      }
    
    if (degenerate_triangle(c1_p, c2_p, c3_p))
      {
	printf("Degenerate triangle! setting s=-1 and skipping triangle. 6\n");
	s=-1;
      }

	  
    //
    // 100223: Hmm... What is really the pattern with respect to the usage of
    //         'forced_skipping_of_second_edge'-values used below? 
    //
    //         Note that if s==0, the triangle will silently be forgotten...
    //
    // 100225: No, adding a stage for this situation too.
    //
    //         Hmm... Maybe it's time to stop letting "on" imply "inside"...?!
    //

#ifdef SPLIT_LARGE_TRIANGLES
    if (  (!(c1_inside || c2_inside || c3_inside))  &&  (!inner_trimming_curve)  &&  (c1_on || c2_on || c3_on)  )
      s=0;
#endif
    



    bool redo=false;
    if (dbg)
      printf("  Going into the switch, s=%d\n", s);
    switch (s)
      {
#ifdef SPLIT_LARGE_TRIANGLES
      case 0:
	if (dbg) printf("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
	if (!inner_trimming_curve)
	  redo = split_triangle_with_all_corners_outside(c1_indx, c2_indx, c3_indx,
							 c1_inside, c2_inside, c3_inside,
							 c1_on, c2_on, c3_on,
							 srf, vert, vert_p, bd, norm, new_mesh, trim_curve_p, contour,
							 dbg, inner_trimming_curve);
	if (dbg) printf("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
	break;
#endif
      case 5: // Corners 1 and 3 are inside, intersections must be on edges 1-2 and 2-3. Compatible with 1-case.
	// 090217: Not now, don't want the 'forced skip'.
 	redo=split_triangle(c1_indx, c2_indx, c3_indx, srf, vert, vert_p, bd, norm, new_mesh,
			    trim_curve_p, contour, dbg, false, inner_trimming_curve);
 	break;
      case 1: // Only corner 1 is inside, so intersections must be on edges 1-2 and 3-1.
	redo=split_triangle(c1_indx, c2_indx, c3_indx, srf, vert, vert_p, bd, norm, new_mesh,
			    trim_curve_p, contour, dbg, true, inner_trimming_curve);
	break;
      case 2: // Corner 2 is inside, intersections must be on edges 2-3 and 3-1. Can go straight into the 3-case.

	// 090217: Huh?! Should be 2-3 and 1-2, shouldn't it?!  Think maybe the comment was wrong, but
	//         conclusion and code correct.  Sadly, this will only work when the assumption about no
	//         "spurious" intersections is met. I.e., that there are no "superfluous"
	//         intersections. In many cases this does not hold. Trying to fix this by not combining
	//         cases here, and at the same time, not letting split_triangle choose to continue on the
	//         third edge if it does not find any intersection on the second edge. The problem seems
	//         to be that such a spurious intersection can be on edge 2, causing split_triangle not to
	//         look for the (correct and sought) intersection on edge 3...

	redo=split_triangle(c2_indx, c3_indx, c1_indx, srf, vert, vert_p, bd, norm, new_mesh, 
			    trim_curve_p, contour, dbg, true, inner_trimming_curve);
	break;
      case 3: // Corners 1 and 2 are inside, so intersections must be on edges 2-3 and 3-1.
	redo=split_triangle(c2_indx, c3_indx, c1_indx, srf, vert, vert_p, bd, norm, new_mesh,
			    trim_curve_p, contour, dbg, false, inner_trimming_curve);
	break;
      case 6: // Corners 2 and 3 are inside, intersections must be on edges 3-1 and 1-2. Compatible with 4-case.
	// 090217: Not now, don't want the 'forced skip'.
	redo=split_triangle(c3_indx, c1_indx, c2_indx, srf, vert, vert_p, bd, norm, new_mesh,
			    trim_curve_p, contour, dbg, false, inner_trimming_curve);
	break;
      case 4: // Corner 3 is inside, so intersections must be on edges 3-1 and 2-3.
	redo=split_triangle(c3_indx, c1_indx, c2_indx, srf, vert, vert_p, bd, norm, new_mesh,
			    trim_curve_p, contour, dbg, true, inner_trimming_curve);
	break;
      case 7: // All corners are inside, not splitting, just keeping it.
	add_triangle(vert_p, new_mesh, c1_indx, c2_indx, c3_indx);
	redo=false;
	break;
      }
    mesh=new_mesh;

    if (dbg)
      printf("  redo=%d\n= trim_a_triangle done ======================================================="
	     "==================\n\n\n", redo);

    return redo;
  }






  //==============================================================================================================
  //
  // More or less just a loop over 'triam_a_triangle' for a list of triangles.
  //
  // 100211: For debugging, a matlab script is written for each triangle in 'mesh', rendering the triangle
  //         plus those produced. Note that (currently) this is only done during the first iteration.
  //
  //         Magenta: This is being processed. One for each tats-file.
  //         Red:     These were produced from the red one.
  //
  //==============================================================================================================

  void trim_a_triangle_soup(shared_ptr<ParamSurface> srf, 
			    vector<Vector3D> &vert,
			    vector<Vector2D> &vert_p, 
			    vector< int > &bd, vector<Vector3D> &norm,
			    const vector<Vector3D> &trim_curve_p,
			    const vector<int> &contour,
			    vector<int> &mesh,
#ifdef DBG
			    vector<char> &col,
#endif
			    const bool inner_trimming_curve = false)
  {
    vector<int> new_mesh, mesh_to_reiterate, mesh_ok;
    const int maxiter=10;
    int iter=0;

#ifdef DBG
    system("rm tats*.m");
#endif

    do
      {
	
	  for (int ilim=(int)mesh.size()/3, i=0; i<ilim; i++)
	//for (int ilim=mesh.size()/3, i=0; i<1; i++)
	  {
	    const int c1_indx=mesh[3*i], c2_indx=mesh[3*i+1], c3_indx=mesh[3*i+2];
	    bool redo=trim_a_triangle(srf, vert, vert_p, bd, norm, trim_curve_p, contour,
				      c1_indx, c2_indx, c3_indx, new_mesh, inner_trimming_curve
#ifdef DBG
				      // 100223: Remember that when having more trimming curves, this may be
				      //         triggered once for each of them...
				      ,false // , i==423 // , false // , i==1069 // , false // , i==94
#endif
	      );
	    if (redo)
	      mesh_to_reiterate.insert(mesh_to_reiterate.end(), new_mesh.begin(), new_mesh.end());
	    else
	      mesh_ok.insert(mesh_ok.end(), new_mesh.begin(), new_mesh.end());
#ifdef DBG
	    // 100210: Showing both the triangle being processed (magenta) and the ones being produced (red.)
	    if ( (iter==0) || (iter==1) )
#if 0
	      if ( 
		
		( ((fabs(vert_p[mesh[3*i  ]][0]-(-22.56))<1e-5) && (fabs(vert_p[mesh[3*i  ]][1]-(30.888))<1e-5)) ||
		  ((fabs(vert_p[mesh[3*i+1]][0]-(-22.56))<1e-5) && (fabs(vert_p[mesh[3*i+1]][1]-(30.888))<1e-5)) ||
		  ((fabs(vert_p[mesh[3*i+2]][0]-(-22.56))<1e-5) && (fabs(vert_p[mesh[3*i+2]][1]-(30.888))<1e-5))    ) &&
		
		(vert_p[mesh[3*i]][0] > -26) &&
		(vert_p[mesh[3*i]][0] < -17) &&
		(vert_p[mesh[3*i]][1] >  28) &&
		(vert_p[mesh[3*i]][1] <  36)    )
#endif
	      {
		char fname[1000];
		if (iter==0)
		  sprintf(fname, "tats%04d.m", i);
		else
		  sprintf(fname, "xtats%04d.m", i);
		FILE *f = fopen(fname, "w");
		fprintf(f, "hold on\ngrid on\n");
		const int u=0, v=1;
		const char col1 = iter==0 ? 'm' : 'b';
		const char col2 = iter==0 ? 'r' : 'g';
		fprintf(f, "patch([%f; %f; %f], [%f; %f; %f], '%c');\n",
			vert_p[mesh[3*i]][u], vert_p[mesh[3*i+1]][u], vert_p[mesh[3*i+2]][u], 
			vert_p[mesh[3*i]][v], vert_p[mesh[3*i+1]][v], vert_p[mesh[3*i+2]][v], col1);
		fprintf(f, "line([%f; %f; %f; %f], [%f; %f; %f; %f], 'color', 'k');\n",
			vert_p[mesh[3*i]][u], vert_p[mesh[3*i+1]][u], vert_p[mesh[3*i+2]][u], vert_p[mesh[3*i]][u], 
			vert_p[mesh[3*i]][v], vert_p[mesh[3*i+1]][v], vert_p[mesh[3*i+2]][v], vert_p[mesh[3*i]][v]);
		for (int j=0; j<int(new_mesh.size())/3; j++)
		  {
		    fprintf(f, "patch([%f; %f; %f], [%f; %f; %f], '%c');\n",
			    vert_p[new_mesh[3*j]][u], vert_p[new_mesh[3*j+1]][u], vert_p[new_mesh[3*j+2]][u], 
			    vert_p[new_mesh[3*j]][v], vert_p[new_mesh[3*j+1]][v], vert_p[new_mesh[3*j+2]][v], col2);
		    fprintf(f, "line([%f; %f; %f; %f], [%f; %f; %f; %f], 'color', 'k');\n",
			    vert_p[new_mesh[3*j  ]][u], vert_p[new_mesh[3*j+1]][u],
			    vert_p[new_mesh[3*j+2]][u], vert_p[new_mesh[3*j  ]][u], 
			    vert_p[new_mesh[3*j  ]][v], vert_p[new_mesh[3*j+1]][v], 
			    vert_p[new_mesh[3*j+2]][v], vert_p[new_mesh[3*j  ]][v]);
		  }
		fprintf(f, "xlabel('magenta/blue=old triangle, red/green=new ones');\n");
		fprintf(f, "hold off\n");
		fclose(f);
	      }
#endif
	  }
	
#ifdef DBG
	if (iter==0)
	  if (inner_trimming_curve)
	    col=vector<char>(mesh_ok.size()/3, 'm');
	  else
	    col=vector<char>(mesh_ok.size()/3, 'c');
#endif
      
	mesh=mesh_to_reiterate;
	iter++;

      }
    while ( (mesh.size()>0) && (iter<maxiter) );

    mesh=mesh_ok;
#ifdef DBG
    while (col.size()<mesh.size()/3)
      if (inner_trimming_curve)
	col.push_back('y');
      else
	col.push_back('b');
#endif
  }






  //==============================================================================================================
  //
  // 090218: Pulling this out in a routine of it's own.
  // 090219: Making sure duplicate points do not survive.
  //
  //==============================================================================================================

  void construct_corner_lists(shared_ptr<ParamSurface> srf, 
			      vector<shared_ptr<ParamCurve> >& crv_set,
			      const vector< Vector3D > &vert,
			      const vector< Vector2D > &vert_p,
			      const vector< Vector3D > &norm,
			      const vector<int> &mesh,
			      const vector< vector< Vector3D > > &trim_curve_all,
			      const vector< vector< Vector3D > > &trim_curve_p_all,
			      const int dn, const int dm, // Will generate an (dn+1)x(dm+1)-mesh...
			      const double bd_res_ratio,
			      const double eps,
			      const double cosangle_limit,
			      const double pt_edge_degen_limit,
			      const double pt_mult_def,
			      const double min_corner_dist, const double max_corner_dist,

			      vector<int> &skip_quad, 
			      vector< vector<Vector3D> > &quad_corner_trim_curve,
			      vector< vector<Vector2D> > &quad_corner_trim_curve_p,
			      int &quad_corners)
  {
    // 090130: "Corner-handling", 1st step: Identify quads needing special treatment, i.e. extra refinement.
    // 090203: Note that we really want to avoid making almost-degenerate triangles...
    //         Hmm... This should really be something relative to the extent of the quad...
    //const double delta=1e-8; // See above.
    skip_quad=vector<int>(dn*dm, 0);
    quad_corner_trim_curve.resize(dn*dm);
    quad_corner_trim_curve_p.resize(dn*dm);

    const RectDomain dom = srf->containingDomain();
    const double u0 = dom.umin();
    const double u1 = dom.umax();
    const double v0 = dom.vmin();
    const double v1 = dom.vmax();
    const double dv=(v1-v0)/dm, du=(u1-u0)/dn;

#ifdef DBG
    // 100222: Producing matlab-scripts for debugging.
    {
      const double u02=vert_p[0][0], v02=vert_p[0][1];
      const double u12=vert_p[(dn+1)*(dm+1)-1][0], v12=vert_p[(dn+1)*(dn+1)-1][1];
      if (fabs(u02-u0)>eps) THROW("Huh?!");
      if (fabs(u12-u1)>eps) THROW("Huh?!");
      if (fabs(v02-v0)>eps) THROW("Huh?!");
      if (fabs(v12-v1)>eps) THROW("Huh?!");
      printf("Parameter domain: [%g, %g] x [%g, %g]\n", u02, u12, v02, v12);
    }
    FILE *f1=fopen("c1.m", "w"); // corners in quads
    FILE *f2=fopen("c2.m", "w"); // corner skipped due to small angle
    FILE *f3=fopen("c3.m", "w"); // corner skipped, is a duplicate
    FILE *f4=fopen("c4.m", "w"); // all points
    fprintf(f1, "hold on\ngrid on\n");
    fprintf(f2, "hold on\ngrid on\n");
    fprintf(f3, "hold on\ngrid on\n");
    fprintf(f4, "hold on\ngrid on\n");
#endif

    for (int curve=0; curve<int(crv_set.size()); curve++)
      {
#ifdef CLC_DBG
	printf("curve %d\n", curve);
#endif
	const vector<Vector3D> &trim_curve_p=trim_curve_p_all[curve];
	const vector<Vector3D> &trim_curve=trim_curve_all[curve];
	Vector2D last_dir(1.0, 0.0); // Will be kept normalized. The initial value will not be used,
	// ref. 'a_first_point_added'
	bool a_first_point_added=false;
	Vector2D last_used_p(0.0, 0.0);
	const int contour_points=(int)trim_curve_p.size();

	for (int i=0; i<contour_points; i++)
	  {
	    double u=trim_curve_p[i][0], v=trim_curve_p[i][1];
//	  const int /* ip1=(i+1)%contour_points, */ im1=(i+contour_points-1)%contour_points;
#ifdef CLC_DBG
	    printf("  i=%4d u=%7.3f, v=%7.3f", i, u, v);
#endif
	    if ((u>=u0) && (u<=u1) && (v>=v0) && (v<=v1))
	      {
		const int p=std::min(int(floor((v-v0)/dv)), dm-1), q=std::min(int(floor((u-u0)/du)), dn-1);
		const int vert_indx=p*(dn+1)+q, quad_indx=p*dn+q;
		const double &u_left=vert_p[vert_indx][0], &v_left=vert_p[vert_indx][1];
		const double &u_right=vert_p[vert_indx+dn+2][0], &v_right=vert_p[vert_indx+dn+2][1];
// 	      if (u<u_left)
// 		{
// 		  printf("\nu=%g, u_left=%g\n", u, u_left);
// 		  THROW("Huh?!");
// 		}
// 	      if (v<v_left) THROW("Huh?!");
// 	      if ((u>u_right) && (q<n-1)) THROW("Huh?!");
// 	      if ((v>v_right) && (p<n-1)) THROW("Huh?!");
		// Note that due to rounding errors, discrete arithmetic, or floor's behaviour, u or v might
		// actually be outside the quads parameter (sub)domain. This should only happen on the right and
		// upper rim...
		// 090219: Doing a small adjustment for this...
		u=std::min(std::max(u_left, u), u_right);
		v=std::min(std::max(v_left, v), v_right);
	      
		// Need two things: 1) angle!=pi, 2) not duplicate point. And 3) not on the edges.
		// 090218: New criterium: Accumulating angle changes.
		// 090220: Taking care of weeding out duplicates in the calling function. (?)
	      
		if ( ((u-u_left<pt_edge_degen_limit) || (u_right-u<pt_edge_degen_limit) || 
		      (v-v_left<pt_edge_degen_limit) || (v_right-v<pt_edge_degen_limit)    ) && (0) )
		  {
		    //printf("  quad-corner too close to edge to split.\n");
		    // 090202: Hmmm... On second thought... We should still split here, it's just that we
		    //         shouldn't split into four new triangles.
		    // 090203: Hmm... On third thought... It's ok *not* to "corner split" here, because
		    //         everything will be ok when we do the "ordinary" curve-splitting later. The
		    //         point about the corners is to force splitting in the corner, and when the
		    //         corner happens to be on a mesh-edge, this will happen automagically!
		    // 090203: Note that even if we end up here, the angle, distance to previous point
		    //         etc. may still not warrant a splitting. This is likely often the case when the
		    //         trimming curve follows an edge of the surface!
		    // 090203: No! If we don't add a corner even though on an edge, splitting may not occur
		    //         later either.  (But why not? Was this not solved with the "boundary-fix" some
		    //         time earlier?)  Another solution is to split, and then remove (or not add)
		    //         degenerate triangles...
		  }
		else
		  {
		    const Vector2D current_p(trim_curve_p[i][0], trim_curve_p[i][1]);
		    const Vector2D previous_p = 
		      i>0 ? Vector2D(trim_curve_p[i-1][0], trim_curve_p[i-1][1]) : Vector2D(1e99, 1e99);
		    const double dist_since_last_used_p=sqrt((current_p-last_used_p)*(current_p-last_used_p));
		    const double dist_since_prev_p=sqrt((current_p-previous_p)*(current_p-previous_p));
		 
#ifdef DBG
		    fprintf(f4, "plot(%f, %f, 'kd', 'markersize', 5);\n", current_p[0], current_p[1]);
#endif
 
		    if ( (dist_since_prev_p>pt_mult_def) // a new distinct point
			 || (!a_first_point_added) )
		      {
			const double cosangle = (last_dir*(current_p-last_used_p))/dist_since_last_used_p;
			// (Note: We have made sure 'dist_since_last_used_p' is non-zero!)
#ifdef CLC_DBG
			printf(" [%7.3f %7.3f   %7.3f %7.3f]", last_dir[0], last_dir[1],
			       current_p[0]-last_used_p[0], current_p[1]-last_used_p[1]);
			printf(" angle=%7.3f", acos(std::min(cosangle, 1.0))/M_PI*180.0);
			printf(" cosangle=%6.3f", cosangle);
#endif
//		      if ( ((cosangle<cosangle_limit) || (b*b>max_corner_dist*max_corner_dist)) && 
//			   (b*b>min_corner_dist*max_corner_dist) )
			if ((cosangle<cosangle_limit) || (!a_first_point_added))
			  {
			    //const Vector2D previous_p(trim_curve_p[im1][0], trim_curve_p[im1][1]);
// 			  if ( (current_p-previous_p)*(current_p-previous_p) < 1e-15*1e-15 )
// 			    {
// #ifdef DBG
// 			      printf(" skipping point, same as previous. ");
// #endif
// 			    }
// 			  else
// 			    {
			    // Now we flag the quad as being split, and append the corner to the quad's list.
			    //printf("dist from previous: %g ", sqrt((current_p-previous_p)*(current_p-previous_p)));
			    skip_quad[quad_indx]=1;
			    quad_corner_trim_curve[quad_indx].push_back(trim_curve[i]);
			    quad_corner_trim_curve_p[quad_indx].push_back(current_p);
			    quad_corners++;
			    last_dir=current_p-last_used_p;
//printf("\n last_dir 1 =%g %g\n", last_dir[0], last_dir[1]);
			    last_dir.normalize();
//printf("\n last_dir 2 =%g %g\n", last_dir[0], last_dir[1]);
			    last_used_p=current_p;
			    a_first_point_added=true;
#ifdef DBG
			    fprintf(f1, "plot(%f, %f, 'm*');\n", current_p[0], current_p[1]);
			    char *tmpstr=new char[100];
			    sprintf(tmpstr, "%d", i);
			    fprintf(f1, "text(%f, %f, '%s');\n", current_p[0], current_p[1]+0.00005, tmpstr);
#  ifdef CLC_DBG
			    printf(" (c1) USED");
#  endif
#endif
// 			    }
			  }
#ifdef DBG
			else
			  {
			    fprintf(f2, "plot(%f, %f, 'b+');\n", current_p[0], current_p[1]);
			    char *tmpstr=new char[100];
			    sprintf(tmpstr, "%d", i);
			    fprintf(f2, "text(%f, %f, '%s');\n", current_p[0], current_p[1]+0.00010, tmpstr);
#  ifdef CLC_DBG
			    printf(" (c2) TOO SMALL ANGLE");
#  endif
			  }
#endif
		      }
		    else
		      {
			// Use it, but not if already used!!!!
			if (dist_since_last_used_p>1e-14)
			  {
			  
#ifdef DBG
			    fprintf(f3, "plot(%f, %f, 'go', 'markersize', 4);\n", current_p[0], current_p[1]);
			    char *tmpstr=new char[100];
			    sprintf(tmpstr, "%d", i);
			    fprintf(f3, "text(%f, %f, '%s');\n", current_p[0], current_p[1]+0.00015, tmpstr);
#  ifdef CLC_DBG
			    printf(" (c3) DUPLICATE: %f", dist_since_last_used_p);
#  endif
#endif
		      
			    // Multiple point, definitely use this one, since we assume it is due to a kink in
			    // the curve.
			    skip_quad[quad_indx]=1;
			    quad_corner_trim_curve[quad_indx].push_back(trim_curve[i]);
			    quad_corner_trim_curve_p[quad_indx].push_back(current_p);
			    quad_corners++;
			    last_dir=current_p-last_used_p;
//printf("\n last_dir 3 =%g %g\n", last_dir[0], last_dir[1]);
			    last_dir.normalize();
//printf("\n last_dir 4 =%g %g\n", last_dir[0], last_dir[1]);
			    last_used_p=current_p;
			    a_first_point_added=true;
			  }
// 		      else
// 			{
// 			  printf("\n@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ %f\n", dist_since_last_used_p);
// 			}
		      
		      }
		  }
	      } // end of if-test ensuring that [u, v] is in the parameter domain.
#ifdef CLC_DBG
	    printf("\n");
#endif
	  } // end of i-loop over points on curve 'curve'.
      } // end of 'curve'-loop
  
#ifdef DBG
    fprintf(f1, "hold off\n");
    fprintf(f2, "hold off\n");
    fprintf(f3, "hold off\n");
    fprintf(f4, "hold off\n");
    fclose(f1);
    fclose(f2);
    fclose(f3);
    fclose(f4);
#endif
  }






  //==============================================================================================================
  //
  // Main routine returning a trimmed mesh in form of a list of triangles.
  //
  // 081206: Adding support for more than one curve.
  //
  //==============================================================================================================

  void make_trimmed_mesh(shared_ptr<ParamSurface> srf,
			 vector<shared_ptr<ParamCurve> >& crv_set,
			 vector< Vector3D > &vert,
			 vector< Vector2D > &vert_p,
			 vector< int > &bd,
			 vector< Vector3D > &norm,
			 //vector< Vector3D > &col,		// colour // 090130: Not in use
			 vector<int> &mesh,
			 vector< Vector3D > &trim_curve_0,	// Output
			 vector< Vector3D > &trim_curve_p_0,	// Output
			 const int dn,				// Initially dn+1 nodes in the u-direction
			 const int dm,				// Initially dm+1 nodes in the v-direction
			 //vector< Vector3D > &extra_v,		// 090130: Not in use
			 double bd_res_ratio)
  {
    const double duplicate_tolerance = 1e-8; // 100210: Absolute number. Was 1e-12.

    vert.resize((dn+1)*(dm+1));
    vert_p.resize((dn+1)*(dm+1));
    bd.resize((dn+1)*(dm+1));
    norm.resize((dn+1)*(dm+1));
    mesh.resize(0);

    int dim = srf->dimension();
    const RectDomain dom = srf->containingDomain();
    const double u0 = dom.umin();
    const double u1 = dom.umax();
    const double v0 = dom.vmin();
    const double v1 = dom.vmax();
    int i, j;
    const double ustep = (u1 - u0)/(dn-1);
    const double vstep = (v1 - v0)/(dm-1);
  
    //--------------------------------------------------------------------------------------------------------------
    //
    // Discretizing the trim curve...
    // 040621: We use the 'trim_curve_p' to store parameter pairs with z=0.0 in the form of triples for the
    //         points in 'trim_curve'.
    // 081208: Now discretizing a set of curves.
    //
    //--------------------------------------------------------------------------------------------------------------

    vector< vector<Vector3D> > trim_curve_all(crv_set.size()), trim_curve_p_all(crv_set.size());
    for (int c=0; c<int(crv_set.size()); c++)
      {
	vector<Vector3D> &trim_curve   = trim_curve_all[c];
	vector<Vector3D> &trim_curve_p = trim_curve_p_all[c];

	shared_ptr<SplineCurve> crv(crv_set[c]->geometryCurve());

	vector<double> corner_pars;
	int kk = crv->order();
	int kn = crv->numCoefs();
	std::vector<double>::iterator st = crv->knotsBegin();
	int kstat = 0;
	corner_pars.push_back(st[kk-1]);
	int knot_ind = kk;
	while (knot_ind < kn)
	  {
	    int knot_mult = 1;
	    while (st[knot_ind] == st[knot_ind+knot_mult])
	      ++knot_mult;
	    if (kstat < 0)
	      {
		THROW("Unexpected incident.");
	      }
	    if (knot_mult > kk - 2) // A kink
	      corner_pars.push_back(st[knot_ind]);
	  
	    knot_ind += knot_mult;
	  }
      
	corner_pars.push_back(st[kn]);
      
	int n2 = 200/(int)(corner_pars.size() - 1); //200; //std::max(200, 4*n); @@sbr Should be const.
	n2 = std::max(n2, 2);

	// 100212: For debugging/testing:
	// int n2 = 30/(corner_pars.size() - 1); //200; //std::max(200, 4*n); @@sbr Should be const.

	for (i = 0; i < int(corner_pars.size()) - 1; ++i)
	  {
	    double t0=corner_pars[i];
	    double t1=corner_pars[i+1];
	    double tstep = (t1 - t0)/(n2-1);
	    //   for (i=0; i<n2; i++)
	    double t = t0;
	    double ref_step = tstep;
	  
	    while (t <= t1)
	      {
		//       double t=i/double(n2);	// Assuming the curve is closed, so we don't
		// need t==1.0...
		//       t = ( crv->et[crv->ik-1] * (1.0-t) + 
		// 	    crv->et[crv->in]   * t          );
	      
		Point result; 
		Point result2(3);
		crv->point(result, t);
	      
		if ((bd_res_ratio > 0.0) && (trim_curve_p.size() > 0) &&
		    ((fabs(trim_curve_p[trim_curve_p.size()-1][0] - result[0]) > bd_res_ratio*ustep) ||
		     (fabs(trim_curve_p[trim_curve_p.size()-1][1] - result[1]) > bd_res_ratio*vstep)))
		  {
		    t -= ref_step; // Return to t of last approved parameter pt.
		    ref_step *= 0.5;
		    t += ref_step;
		  }
		else
		  // 090219: We only allow a point to be added if it is distinctly different from the last one
		  //         added!
		  {
//   		  if ( (trim_curve_p.size()==0) ||
//   		       ( ((result[0]-trim_curve_p[trim_curve_p.size()-1][0])*
//   			  (result[1]-trim_curve_p[trim_curve_p.size()-1][1]) ) > 1e-14*1e-14 ) )
		    {
		      trim_curve_p.push_back(Vector3D(result[0], result[1], 0.0));
// 		      printf("added %g %g\n", result[0], result[1]);
		      srf->point(result2, result[0], result[1]);
		      if (dim == 2)
			result2[2] = 0.0;
		      trim_curve.push_back(Vector3D(result2.begin()));
		    }
		    if (t == t1)
		      break;
		    t += tstep;
		    if (t > t1) {
		      tstep = t1 - (t - tstep);
		      t = t1;
		    }
		    ref_step = tstep;
		  }
// 	      printf("==================== t=%g\n", t);
	      }
	  }
      } // closing of loop over the curves



#ifdef DBG
    // 100222: Just to separate different runs from one another...
    printf("\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n");
    printf("\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n");
    printf("\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n");
#endif

  


#ifdef DBG
    puts("====================================================================================================");
    printf("Weeding out duplicates from contours, tolerance=%e\n", duplicate_tolerance);
#endif
    
    //--------------------------------------------------------------------------------------------------------------
    //
    // 090220: Weeding out duplicates. Note that we want the duplicates taken out for the lists used for the
    //         insideness-testing routines, but we want them in for the corner construction! Hence... making two
    //         sets of lists...
    // 100212: There was a bug here, when more than two points were equal, duplicates would survive.
    //
    //--------------------------------------------------------------------------------------------------------------

    vector< vector<Vector3D> > trim_curve_all_orig=trim_curve_all, trim_curve_p_all_orig=trim_curve_p_all;

    for (int c=0; c<int(crv_set.size()); c++)
      {
 	// printf("Length before: %d %d\n", int(trim_curve_all[c].size()), int(trim_curve_p_all[c].size()));

	ASSERT2(trim_curve_all[c].size()>0, printf("Should not be here with degenerate curve.\n"));
	ASSERT2(trim_curve_p_all[c].size()>0, printf("Should not be here with degenerate curve.\n"));

	vector<Vector3D> trim_curve, trim_curve_p;
	trim_curve.push_back(trim_curve_all[c][0]);
	trim_curve_p.push_back(trim_curve_p_all[c][0]);
	
	for (int i=1; i<int(trim_curve_all[c].size()); i++)
	  {
	    const double last_u = trim_curve_p[trim_curve.size()-1][0];
	    const double last_v = trim_curve_p[trim_curve.size()-1][1];
	    const double dist_squared = ( (trim_curve_p_all[c][i][0]-last_u) * (trim_curve_p_all[c][i][0]-last_u) + 
					  (trim_curve_p_all[c][i][1]-last_v) * (trim_curve_p_all[c][i][1]-last_v)   );
	    if ( dist_squared > duplicate_tolerance*duplicate_tolerance )
	      {
		// printf("%3d:%5d: prev point was (%f, %f), new is (%f, %f), dist=%e, pushing.\n",
		// c, i, last_u, last_v, trim_curve_p_all[c][i][0], trim_curve_p_all[c][i][1], sqrt(dist_squared));
		trim_curve.push_back(trim_curve_all[c][i]);
		trim_curve_p.push_back(trim_curve_p_all[c][i]);
	      }
	    // else
	    // printf("%3d:%5d: dist = %e, point skipped\n", c, i, sqrt(dist_squared));
	  }

	// printf("Length after:  %d %d\n", int(trim_curve.size()), int(trim_curve_p.size()));

	// 100212: Must now handle the first point! It can still be coincident with the last one.
	const double last_u = trim_curve_p[trim_curve.size()-1][0];
	const double last_v = trim_curve_p[trim_curve.size()-1][1];
	const double dist_squared = ( (trim_curve_p[0][0]-last_u) * (trim_curve_p[0][0]-last_u) + 
				      (trim_curve_p[0][1]-last_v) * (trim_curve_p[0][1]-last_v)   );
	if ( dist_squared <= duplicate_tolerance*duplicate_tolerance )
	  {
	    trim_curve.pop_back();
	    trim_curve_p.pop_back();
	  }

	// Finally, we replace the original curves with the new ones.
	trim_curve_p_all[c] = trim_curve_p;
	trim_curve_all[c] = trim_curve;

	// printf("Length after:  %d %d\n", int(trim_curve.size()), int(trim_curve_p.size()));
      }



  
    //--------------------------------------------------------------------------------------------------------------
    //
    // 100212: A check to make sure we have no duplicates.
    //
    //--------------------------------------------------------------------------------------------------------------

    {
      bool duplicates_found = false;
      for (int c=0; c<int(crv_set.size()); c++)
	{
	  vector<Vector3D> &trim_curve   = trim_curve_all[c];
	  vector<Vector3D> &trim_curve_p = trim_curve_p_all[c];
	  
	  for (int i=0; i<int(trim_curve.size()); i++)
	    {
		const int im1 = (i-1+(int)trim_curve.size())%(int)trim_curve.size();
	      const double dist_to_previous = 
		sqrt( (trim_curve_p[i][0]-trim_curve_p[im1][0])*(trim_curve_p[i][0]-trim_curve_p[im1][0]) + 
		      (trim_curve_p[i][1]-trim_curve_p[im1][1])*(trim_curve_p[i][1]-trim_curve_p[im1][1])   );

	      // printf("%2d:%5d: dist = %e\n", c, i, dist_to_previous);
	      if ( dist_to_previous <= duplicate_tolerance )
		{
#ifdef DBG
		  printf("%2d:%5d: DUPLICATE POINT! dist_to_previous = %e\n", c, i, dist_to_previous);
#endif
		  duplicates_found = true;
		}
	    }
	}
      if (duplicates_found)
	MESSAGE("\n\n  WARNING: Something is fishy. At this point all duplicates should have been eliminated.\n"
		"  This may lead to a funny mesh...\n\n");
    }
    
    
    

    //--------------------------------------------------------------------------------------------------------------
    //
    // Setting up the structure for fast inside-testing, and the 'contour' cursor array.
    // 081208: Again, extending to a set of contours...
    // 100212: The "fast inside-testing" stuff is long gone.
    //
    //--------------------------------------------------------------------------------------------------------------

    vector< vector<int> > contour_all(crv_set.size());
    for (int c=0; c<int(crv_set.size()); c++)
      {
	vector<int> &contour = contour_all[c];
	const vector<Vector3D> &trim_curve_p = trim_curve_p_all[c];
      
	contour.resize(trim_curve_p.size());
	for (i=0; i<(int)trim_curve_p.size(); i++)
	  contour[i]=3*i;
      }
  
#ifdef DBG
    // Just printing out matlab-code for rendering of all contours.
    if (1)
      for (int c=0; c<int(crv_set.size()); c++)
	{
	  const vector<int> &contour           = contour_all[c];
	  const vector<Vector3D> &trim_curve_p = trim_curve_p_all[c];

	  char *tmp=new char[1000];
	  sprintf(tmp, "d%d.m", c);
	  FILE *f=fopen(tmp, "w");
	  for (int i=0; i<int(trim_curve_p.size())-1; i++)
	    {
	      const double &x0 = trim_curve_p[contour[i]/3  ][0], &y0=trim_curve_p[contour[i]/3  ][1];
	      const double &x1 = trim_curve_p[contour[i]/3+1][0], &y1=trim_curve_p[contour[i]/3+1][1];
	      fprintf(f, "line([%f; %f], [%f; %f], 'color', 'b', 'linewidth', 3);\n", x0, x1, y0, y1);
	    }
	  fclose(f);
	}
#endif


  
    //--------------------------------------------------------------------------------------------------------------
    //
    // Cache'ing all 'is_inside' results...
    // 081208: Again, extending to a set of contours...
    //
    // 100210: It seems like a (harmless) "bug" that the 'vert-arrays are filled once for every curve.
    //
    //--------------------------------------------------------------------------------------------------------------

    vector< vector<int> > inside_all(crv_set.size());
    for (int c=0; c<int(crv_set.size()); c++)
      {
	const vector<int> &contour = contour_all[c];
	const vector<Vector3D> &trim_curve_p = trim_curve_p_all[c];
	vector<int> &inside = inside_all[c];
	inside = vector<int>((dn+1)*(dm+1));
	double uv[2], s;
	vector<Point> res(3);
	Point nrm;

	ASSERT2(dim==2 || dim==3, printf("Huh?! dim=%d\n", dim));

	for (i=0; i<=dm; i++)
	  {
	    double t=i/double(dm);
	    uv[1]=v0*(1.0-t) + v1*t;
	    for (j=0; j<=dn; j++)
	      {
		s=j/double(dn);
		uv[0]=u0*(1.0-s) + u1*s;
	      
		inside[i*(dn+1)+j] = is_inside(trim_curve_p, contour, uv[0], uv[1]);
		
		// 090115:
		if (s2m_with_boundary)
		  inside[i*(dn+1)+j] |= is_on_contour(trim_curve_p, contour, uv[0], uv[1]);

		// 090203: Enable this to see the trimming curve inside the untrimmed surface, for debugging purposes.
		// 100210: This does not work, or this enabling is not enough.
		// inside[i*(dn+1)+j] = 1; // !!!
		
		srf->point(res, uv[0], uv[1], 1);
		if (dim == 2)
		  nrm = Point(0.0, 0.0, 1.0);
		else
		  nrm = res[1].cross(res[2]);

		// 100210: Why on earth is this done for every curve?!
		if (dim == 3)
		  vert[i*(dn+1)+j] = Vector3D(res[0].begin());
		else
		  vert[i*(dn+1)+j] = Vector3D(res[0][0], res[0][1], 0.0);
		vert_p[i*(dn+1)+j] = Vector2D(uv);

		bd[i*(dn+1)+j]= (i==0 || i==dm || j==0 || j==dn) ? 1 : 0;
		if (nrm.length() < 1.0e-12)
		  nrm.setValue(0.0, 0.0, 1.0);
		norm[i*(dn+1)+j] = Vector3D(nrm.begin());
		norm[i*(dn+1)+j].normalize();
	      }
	  }
      }


#ifdef DBG
    {
      FILE *f5=fopen("r1.m", "w"); // regular mesh
      fprintf(f5, "hold on\ngrid on\n");
      for (int i=0; i<dm; i++)
	for (int j=0; j<=dn; j++)
	  fprintf(f5, "line([%f; %f], [%f; %f], 'color', 'k');\n", 
		  vert_p[i*(dn+1)+j][0], vert_p[(i+1)*(dn+1)+j][0],
		  vert_p[i*(dn+1)+j][1], vert_p[(i+1)*(dn+1)+j][1]);
      for (int i=0; i<=dm; i++)
	for (int j=0; j<dn; j++)
	  fprintf(f5, "line([%f; %f], [%f; %f], 'color', 'k');\n", 
		  vert_p[i*(dn+1)+j][0], vert_p[i*(dn+1)+j+1][0],
		  vert_p[i*(dn+1)+j][1], vert_p[i*(dn+1)+j+1][1]);
      fprintf(f5, "hold off\n");
      fclose(f5);
    }
#endif


  
    //--------------------------------------------------------------------------------------------------------------
    //
    // 090130: We want to preserve corner points on the contour. The idea: Identify all quads containing such
    //         corners, then treat these quads separately: Split in the first corner, producing four new
    //         triangles. These must then be processed further, until they contain no corner points. Then these
    //         triangles must still be processed the "old way".
    // 
    // 090203: Note that by checking that we do not make degenerate primitives, it does not become as important
    //         to check for duplicate corner points any more.
    //
    //--------------------------------------------------------------------------------------------------------------


    const double pt_edge_degen_limit=1e-6;	// 090203: Should really be relative to the extent of the
  						//         geometry.  Used to avoid constructing
  						//         almost-degenerate polygons. Applied in the
  						//         parameter domain.

//     const double cosangle_limit=	 	// 090203: The angle must be greater than this for a "corner"
//      cos(M_PI/180.0);
     const double cosangle_limit=		// 090203: The angle must be greater than this for a "corner"
       cos(M_PI/180.0*0.1);			//         to be defined as a corner, i.e., for splitting to
						//         be done.

    const double pt_mult_def=1e-14;		// 090203: In the parameter domain, two points closer than
						//         this are identified.
  
    const double eps=1e-13;                       // 090203: Used for zero-tests for distances in the parameter domain.

    //const double tau=1e-12;                     // 090203: Used for zero-tests for distances in the geometric space.

    const double max_corner_dist=			// 090218: When corners are this far away, they are used no matter what
      100.0*std::max((u1-u0)/dn, (v1-v0)/dm);		//         the angle is. Must be larger than 'min_corner_dist'.

    const double min_corner_dist=			// 090219: When corners are this close, they are not used no
      0.5*std::min((u1-u0)/dn, (v1-v0)/dm);		//         matter what the angle is.





#ifdef DBG
    puts("====================================================================================================");
    puts("Generating list of corners to force splitting in.");
#endif



    // 090218: Each of these will have dim nxn. 
    vector<int> skip_quad;				// Contains a 1 if quad is to be skipped, i.e., it is
							// replaced by custom list of triangles, since it did
							// contain corners or such.

    vector< vector<Vector3D> > quad_corner_trim_curve;	// The 3D corner point, fetched from the proper
							// trimming curve. Is this really used afterward?
							// Yes. Need it for updating 'vert'.

    vector< vector<Vector2D> > quad_corner_trim_curve_p;	// The 2D parameter corner point, constructed from a
    // 3D point by ignoring the third coordinate. Fetched
    // from the proper trimming curve. Not smart to call
    // it '...trim_curve...', should have been just
    // 'quad_corner_p' or something...

    int quad_corners=0;
    construct_corner_lists(srf, crv_set, vert, vert_p, norm, mesh,
			   trim_curve_all_orig, trim_curve_p_all_orig,
			   dn, dm, bd_res_ratio, eps, cosangle_limit, pt_edge_degen_limit, pt_mult_def, min_corner_dist,
			   max_corner_dist, skip_quad, quad_corner_trim_curve, quad_corner_trim_curve_p, quad_corners);
  
  



    //--------------------------------------------------------------------------------------------------------------
    //
    // 100210: Ok, in this stage, quads containing corners are split. They are initially split into four
    //         triangles, which may again be further split in the innermost loop. This is to take care of cases
    //         where quads contain more than just one "corner points".
    //
    //--------------------------------------------------------------------------------------------------------------

#ifdef CORNER_SPLITTING

#ifdef DBG
    puts("====================================================================================================");
    printf("\"Pre-splitting\" mesh, %d corners was found.\n", quad_corners);
#endif
  
    //
    // 090206: We build structures like 'quad_corner' for needed data. The problem with the "merged list"
    //         approach is that information for each curve, like where it ends, must be preserved. This is
    //         because the curves are closed.
    //

    // Now splitting first quads, then triangles.

#ifdef DBG
    printf("dn=%d, dm=%d\n", dn, dm);
#endif
    for (int i=0; i<dn*dm; i++)
      {
	//       printf("i=%3d: corners: %d\n", i, int(quad_corner[i].size()));
	if (quad_corner_trim_curve_p[i].size()>0)
	  {
	    // First, split the quad.
	  
	    const int p=i/dn, q=i%dn;
	  
	    vert.push_back(quad_corner_trim_curve[i][0]);
	    vert_p.push_back(quad_corner_trim_curve_p[i][0]);
	    bd.push_back(1);
	    // Unfortunately, we don't have any normals readily available. What we can do, is to just take the
	    // average of the quad. What we *should* do, is to evaluate the surface...
	    {
	      Vector3D new_norm = norm[p*(dn+1)+q] + norm[p*(dn+1)+q+1] + norm[(p+1)*(dn+1)+q] + norm[(p+1)*(dn+1)+q+1];
	      new_norm.normalize();
	      norm.push_back(new_norm);
	    }
	  
	    int first_new_triangle=(int)mesh.size()/3;
	  
	    // printf("i=%d, p=%d, q=%d, ndx=%d: splitting this quad in %g, %g.\n", i, p, q, p*(n+1)+q,
	    //	  vert_p[vert_p.size()-1][0], vert_p[vert_p.size()-1][1]);
	  
	    add_triangle(vert_p, mesh, (int)vert.size()-1, (p+1)*(dn+1)+q  ,  p   *(dn+1)+q  );
	    add_triangle(vert_p, mesh, (int)vert.size()-1, (p+1)*(dn+1)+q+1, (p+1)*(dn+1)+q  );
	    add_triangle(vert_p, mesh, (int)vert.size()-1,  p   *(dn+1)+q+1, (p+1)*(dn+1)+q+1);
	    add_triangle(vert_p, mesh, (int)vert.size()-1,  p   *(dn+1)+q  ,  p   *(dn+1)+q+1);
	  
	    // Now, if there are more corners, we have to go through the triangles. And we have to do this again
	    // for each new corner, and apply the split-testing and eventual splitting to all newly generated
	    // triangles also. We assume that only one triangle can contain a given corner, so when a corner is
	    // used, it does not need to be tested for again.
	  
	    if (int(quad_corner_trim_curve_p[i].size())>1)
	      {
	      
		//printf("  Corners in total: %d, continuing with triangle-splitting.\n", int(quad_corner[i].size()));
	      
		for (int j=1; j<int(quad_corner_trim_curve_p[i].size()); j++)
		  {
		    bool done_splitting=false;
		    // 100210: Seems to be enough to enable this to skip the triangle-splitting stage.
		    // done_splitting=true;

		    Vector2D &corner = quad_corner_trim_curve_p[i][j];
		    //printf("    j=%d, checking corner %f, %f.\n", j, corner[0], corner[1]);
		  
		    // Ok, must now check, and possibly split triangles from 'first_new_triangle'.
		    // Note that we quit after we split, since we want to start over again asap.
		    for (int k=first_new_triangle; (k<int(mesh.size())/3) && (!done_splitting); k++)
		      {
			//printf("      k=%d\n", k);
			if (
			  pt_inside_tri(corner, vert_p[mesh[3*k]], vert_p[mesh[3*k+1]], vert_p[mesh[3*k+2]])
			  // && (pt_dist_from_tri_edge(corner, vert_p[mesh[3*k]], 
			  // vert_p[mesh[3*k+1]], vert_p[mesh[3*k+2]])
			  // > pt_edge_degen_limit)
			  // 090204: This test caused corners on the outer boundary not to trigger
			  //         splitting. Enabling this test is probably not a good idea, unless these
			  //         splitting causes other problems...
			  )
			  // Splitting the triangle.
			  {
			    //printf("      splitting, i=%d, j=%d, k=%d.\n", i, j, k);
			    vert.push_back(quad_corner_trim_curve[i][j]);
			    vert_p.push_back(corner);
			    bd.push_back(1);
			    const int a=mesh[3*k], b=mesh[3*k+1], c=mesh[3*k+2];
			    // Unfortunately, we don't have any normals readily available. What we can do, is
			    // to just take the average of the quad. What we should do, is to evaluate the
			    // surface...
			    {
			      Vector3D tmp = norm[a] + norm[b] + norm[c];
			      tmp.normalize();
			      norm.push_back(tmp);
			    }
			    mesh.erase(mesh.begin()+3*k, mesh.begin()+3*k+3);
			    add_triangle(vert_p, mesh, (int)vert.size()-1, a, b);
			    add_triangle(vert_p, mesh, (int)vert.size()-1, b, c);
			    add_triangle(vert_p, mesh, (int)vert.size()-1, c, a);
			    done_splitting=true; // No need to check more triangles for a corner already used.
			    // 090204: Hmm... maybe not entirely true. If splitting was done on an edge, there may
			    //         be a triangle on the other side of that edge that also needs splitting...?!
			  }
		      }
		  }
	      }
	  }
      }
  
#endif // end of #ifdef CORNER_SPLITTING

#ifdef DBG
    vector<char> col(mesh.size()/3, 'g');
#endif
  



#ifdef SPLIT_NON_CORNER_QUADS

#ifdef DBG
    puts("====================================================================================================");
    puts("Splitting remaining quads according to the first trimming curve.");
#endif

    for (i=0; i<dm; i++)
      for (j=0; j<dn; j++)
	if (!skip_quad[i*dn+j])
	  {
	    // if ( (inside1[i*(n+1)+j]!=0) && (inside1[i*(n+1)+j]!=1) )
	    // MESSAGE("noe er muffens, inside1 er hverken 0 eller 1?!");
	    // if ( (inside[i*(n+1)+j]!=0) && (inside[i*(n+1)+j]!=1) )
	    // MESSAGE("noe er muffens, inside er hverken 0 eller 1?!");
	  
	    int s=((inside_all[0][ i   *(dn+1)+j  ] << 3) +	// 8
		   (inside_all[0][ i   *(dn+1)+j+1] << 2) +	// 4
		   (inside_all[0][(i+1)*(dn+1)+j  ] << 1) +	// 2
		   (inside_all[0][(i+1)*(dn+1)+j+1] << 0));	// 1
// 	  if (s!=0) // ((i==6) && (j==17))
// 	    printf("s=%d, i=%d, j=%d\n", s, i, j);
	  
// 090216: keeping all to test if the missing triangle was made at all. Ok, it was made, but then wrongly split here.
//s=15;
	  
// 090216: checking the lower right corner for debugging of problematic quad...
//if ((fabs(vert_p[i*(n+1)+j+1][0]-.01)<1e-5) && (fabs(vert_p[i*(n+1)+j+1][1]-0.65)<1e-5))
//  if ((i==12) && (j==3))
//    printf("oooooooooooooooooooooooooooooooooooooooooooooooooo> s=%d\n", s), s=15;
	  
	    // 100223: Shit! Didn't realize this override was here... This means that all the quad-splitting
	    //         code below is actually not used... Wonder when this was done...
	    s=15;
	    
	    switch (s)
	      {
	      case 15:
		// The whole rectangle is inside.
	      {
		add_triangle(vert_p, mesh, (i+1)*(dn+1)+j, i*(dn+1)+j  ,  i   *(dn+1)+j+1);
		add_triangle(vert_p, mesh, (i+1)*(dn+1)+j, i*(dn+1)+j+1, (i+1)*(dn+1)+j+1);
	      }
	      break;
	      //-------------------------------------------------------------------
	      case 10:
		// The two left (smaller u, smaller j) points are inside.
	      {
		//               first segm this col
		vector<int> si0; si0.push_back(i  ); si0.push_back(i+1);
		vector<int> sj0; sj0.push_back(j  ); sj0.push_back(j  );
		vector<int> si1; si1.push_back(i  ); si1.push_back(i+1);
		vector<int> sj1; sj1.push_back(j+1); sj1.push_back(j+1);
		//
		//   s0[1]   B    s1[1]	(Husk: i+1 er her!)
		//      o----x------o
		//      |           |
		//      |           |
		//      |           |
		//      o----x------o
		//   s0[0]   A    s1[0]	(Husk: i er her!)
		//
		split_quad(i+1, j,		// s0[1]
			   i, j,	 	// s0[0]
			   -1, -1,
			   // Then the two intersection points are
			   si0, sj0, si1,	// somewhere on s0[0]-s1[0] (A) and
			   sj1,		// somewhere on s0[1]-s1[1] (B).
			   srf, vert, vert_p, bd, norm, mesh, dn, dm,
			   trim_curve_p_all[0], contour_all[0], s2m_with_boundary);
	      }
	      break;
	      //-------------------------------------------------------------------
	      case 5:
		// The two right (larger u, larger j) points are inside.
	      {
		//               first segm this col
		vector<int> si0; si0.push_back(i+1); si0.push_back(i  );
		vector<int> sj0; sj0.push_back(j  ); sj0.push_back(j  );
		vector<int> si1; si1.push_back(i+1); si1.push_back(i  );
		vector<int> sj1; sj1.push_back(j+1); sj1.push_back(j+1);
		//
		//   s0[0]   A    s1[0]	(Husk: i+1 er her!)
		//      o----x------o
		//      |           |
		//      |           |
		//      |           |
		//      o----x------o
		//   s0[1]   B    s1[1]	(Husk: i er her!)
		//
		split_quad(i, j+1,		// s1[1]
			   i+1, j+1,	// s1[0]
			   -1, -1,
			   // Then the two intersection points are
			   si0, sj0, si1,	// somewhere on s0[0]-s1[0] (A) and
			   sj1,		// somewhere on s0[1]-s1[1] (B).
			   srf, vert, vert_p, bd, norm, mesh, dn, dm,
			   trim_curve_p_all[0], contour_all[0], s2m_with_boundary);
	      }
	      break;
	      //-------------------------------------------------------------------
	      case 12:
		// The two lower (smaller v, smaller i) points are inside.
	      {
		// (si0[0], sj0[0]) - (si1[0], sj1[0]) is the first segment, etc.
		//               first segm this col
		vector<int> si0; si0.push_back(i  ); si0.push_back(i+1);
		vector<int> sj0; sj0.push_back(j+1); sj0.push_back(j  );
		vector<int> si1; si1.push_back(i+1); si1.push_back(i  );
		vector<int> sj1; sj1.push_back(j+1); sj1.push_back(j  );
		//
		//   s0[1]        s1[0]	(Husk: i+1 er her!)
		//      o-----------o
		//      |           |
		//    B x           x A
		//      |           |
		//      o-----------o
		//   s1[1]        s0[0]	(Husk: i er her!)
		//
		split_quad(i, j,		// s1[1]
			   i, j+1,	 	// s0[0]
			   -1, -1,
			   // Then the two intersection points are
			   si0, sj0, si1,	// somewhere on s0[0]-s1[0] (A) and
			   sj1,		// somewhere on s0[1]-s1[1] (B).
			   srf, vert, vert_p, bd, norm, mesh, dn, dm,
			   trim_curve_p_all[0], contour_all[0], s2m_with_boundary);
	      }
	      break;
	      //-------------------------------------------------------------------
	      case 1:
		// The upper, right (larger i, larger j) point is inside.
	      {
		//               first segm this col
		vector<int> si0; si0.push_back(i+1); si0.push_back(i+1);
		vector<int> sj0; sj0.push_back(j+1); sj0.push_back(j+1);
		vector<int> si1; si1.push_back(i+1); si1.push_back(i  );
		vector<int> sj1; sj1.push_back(j  ); sj1.push_back(j+1);
		//
		//   s1[0]     A  s0[0]=s0[1]
		//      o-------x---o
		//      |        \  |
		//      |         \ x B
		//      |           |
		//      o-----------o
		//                s1[1]
		//
		split_quad(i+1, j+1,	// s0[0]
			   -1, -1, -1, -1,
			   // Then the two intersection points are
			   si0, sj0, si1,	// somewhere on s0[0]-s1[0] (A) and
			   sj1,		// somewhere on s0[1]-s1[1] (B).
			   srf, vert, vert_p, bd, norm, mesh, dn, dm,
			   trim_curve_p_all[0], contour_all[0], s2m_with_boundary);
	      }
	      break;
	      //-------------------------------------------------------------------
	      case 4:
		// The lower, right (smaller i, larger j) point is inside.
	      {
		//               first segm this col
		vector<int> si0; si0.push_back(i+1); si0.push_back(i  );
		vector<int> sj0; sj0.push_back(j+1); sj0.push_back(j  );
		vector<int> si1; si1.push_back(i  ); si1.push_back(i  );
		vector<int> sj1; sj1.push_back(j+1); sj1.push_back(j+1);
		split_quad(i, j+1,
			   -1, -1, -1, -1,
			   // Then the two intersection points are
			   si0, sj0, si1,	// somewhere on s0[0]-s1[0] (A) and
			   sj1,		// somewhere on s0[1]-s1[1] (B).
			   srf, vert, vert_p, bd, norm, mesh, dn, dm,
			   trim_curve_p_all[0], contour_all[0], s2m_with_boundary);
	      }
	      break;
	      //-------------------------------------------------------------------
	      case 8:
		// The lower, left (smaller i, smaller j) point is inside.
	      {
		//               first segm this col
		vector<int> si0; si0.push_back(i  ); si0.push_back(i  );
		vector<int> sj0; sj0.push_back(j  ); sj0.push_back(j  );
		vector<int> si1; si1.push_back(i  ); si1.push_back(i+1);
		vector<int> sj1; sj1.push_back(j+1); sj1.push_back(j  );
		split_quad(i, j,
			   -1, -1, -1, -1,
			   // Then the two intersection points are
			   si0, sj0, si1,	// somewhere on s0[0]-s1[0] (A) and
			   sj1,		// somewhere on s0[1]-s1[1] (B).
			   srf, vert, vert_p, bd, norm, mesh, dn, dm,
			   trim_curve_p_all[0], contour_all[0], s2m_with_boundary);
	      }
	      break;
	      //-------------------------------------------------------------------
	      case 2:
		// The upper, left (larger i, smaller j) point is inside.
	      {
		//               first segm this col
		vector<int> si0; si0.push_back(i+1); si0.push_back(i+1);
		vector<int> sj0; sj0.push_back(j  ); sj0.push_back(j  );
		vector<int> si1; si1.push_back(i  ); si1.push_back(i+1);
		vector<int> sj1; sj1.push_back(j  ); sj1.push_back(j+1);
		split_quad(i+1, j,
			   -1, -1, -1, -1,
			   // Then the two intersection points are
			   si0, sj0, si1,	// somewhere on s0[0]-s1[0] (A) and
			   sj1,		// somewhere on s0[1]-s1[1] (B).
			   srf, vert, vert_p, bd, norm, mesh, dn, dm,
			   trim_curve_p_all[0], contour_all[0], s2m_with_boundary);
	      }
	      break;
	      //-------------------------------------------------------------------
	      case 3:
		// The two upper points are inside.
	      {
		// Same as for 12, but swapped the two segments.
		// Also using (i+1, j+1) and (i+1, j) as the inside points.
		vector<int> si0; si0.push_back(i+1); si0.push_back(i  );
		vector<int> sj0; sj0.push_back(j  ); sj0.push_back(j+1);
		vector<int> si1; si1.push_back(i  ); si1.push_back(i+1);
		vector<int> sj1; sj1.push_back(j  ); sj1.push_back(j+1);
		split_quad(i+1, j+1,	// s1[1]
			   i+1, j, 		// s0[0]
			   -1, -1,
			   // Then the two intersection points are
			   si0, sj0, si1,	// somewhere on s0[0]-s1[0] (A) and
			   sj1,		// somewhere on s0[1]-s1[1] (B).
			   srf, vert, vert_p, bd, norm, mesh, dn, dm,
			   trim_curve_p_all[0], contour_all[0], s2m_with_boundary);
	      }
	      break;
	      //-------------------------------------------------------------------
	      case 14:
		// All except the upper right point are inside.
	      {
		//               first segm this col
		vector<int> si0; si0.push_back(i+1); si0.push_back(i+1);
		vector<int> sj0; sj0.push_back(j+1); sj0.push_back(j+1);
		vector<int> si1; si1.push_back(i  ); si1.push_back(i+1);
		vector<int> sj1; sj1.push_back(j+1); sj1.push_back(j  );
		split_quad(i, j,		//
			   i+1, j, 	//
			   i, j+1,	//
			   // Then the two intersection points are
			   si0, sj0,	// somewhere on s0[0]-s1[0] (A) and
			   si1, sj1,	// somewhere on s0[1]-s1[1] (B).
			   srf, vert, vert_p, bd, norm, mesh, dn, dm,
			   trim_curve_p_all[0], contour_all[0], s2m_with_boundary);
	      }
	      break;
	      //-------------------------------------------------------------------
	      case 11:
		// All except the lower right point are inside.
	      {
		//               first segm this col
		vector<int> si0; si0.push_back(i  ); si0.push_back(i  );
		vector<int> sj0; sj0.push_back(j+1); sj0.push_back(j+1);
		vector<int> si1; si1.push_back(i  ); si1.push_back(i+1);
		vector<int> sj1; sj1.push_back(j  ); sj1.push_back(j+1);
		split_quad(i+1, j,	//
			   i+1, j+1, 	//
			   i, j,		//
			   // Then the two intersection points are
			   si0, sj0,	// somewhere on s0[0]-s1[0] (A) and
			   si1, sj1,	// somewhere on s0[1]-s1[1] (B).
			   srf, vert, vert_p, bd, norm, mesh, dn, dm,
			   trim_curve_p_all[0], contour_all[0], s2m_with_boundary);
	      }
	      break;
	      //-------------------------------------------------------------------
	      case 7:
		// All except the lower left point are inside.
	      {
		//               first segm this col
		vector<int> si0; si0.push_back(i  ); si0.push_back(i  );
		vector<int> sj0; sj0.push_back(j  ); sj0.push_back(j  );
		vector<int> si1; si1.push_back(i+1); si1.push_back(i  );
		vector<int> sj1; sj1.push_back(j  ); sj1.push_back(j+1);
		split_quad(i+1, j+1,	//
			   i, j+1, 	//
			   i+1, j,	//
			   // Then the two intersection points are
			   si0, sj0,	// somewhere on s0[0]-s1[0] (A) and
			   si1, sj1,	// somewhere on s0[1]-s1[1] (B).
			   srf, vert, vert_p, bd, norm, mesh, dn, dm,
			   trim_curve_p_all[0], contour_all[0], s2m_with_boundary);
	      }
	      break;
	      //-------------------------------------------------------------------
	      case 13:
		// All except the upper left point are inside.
	      {
		//               first segm this col
		vector<int> si0; si0.push_back(i+1); si0.push_back(i+1);
		vector<int> sj0; sj0.push_back(j  ); sj0.push_back(j  );
		vector<int> si1; si1.push_back(i+1); si1.push_back(i  );
		vector<int> sj1; sj1.push_back(j+1); sj1.push_back(j  );
		split_quad(i, j+1,	//
			   i, j, 	//
			   i+1, j+1,	//
					// Then the two intersection points are
			   si0, sj0,	// somewhere on s0[0]-s1[0] (A) and
			   si1, sj1,	// somewhere on s0[1]-s1[1] (B).
			   srf, vert, vert_p, bd, norm, mesh, dn, dm,
			   trim_curve_p_all[0], contour_all[0], s2m_with_boundary);
	      }
	      break;
	      //-------------------------------------------------------------------
	      }
	  }

#ifdef DBG
    while (col.size()<mesh.size()/3)
      col.push_back('r');
#endif
  
#endif // end of #ifdef SPLIT_NON_CORNER_QUADS





#ifdef SPLIT_TRIANGLES_AT_OUTER_CURVE

#  ifdef DBG
    puts("====================================================================================================");
    puts("Splitting all triangles just generated, according to the first trimming curve.");
#  endif

    //
    // 090203: Here we need a new pass to split triangles produced by the new "corner-stuff", and now we should
    //         split according to the first trimming curve, only. (Note, this all is based on the assumption
    //         that we should keep whatever is *inside* curve 0, and *outside* curves 1, 2, ...)
    // 
    // 090203: Hmm... maybe this code can somehow be combined with the stage below, for more compact code?!
    //
    // 090203: Have a suspicion that h.g2 produces a denegerate triangle in this part of the code.
    //

    for (int c=0; c<std::min(1, int(crv_set.size())); c++)
      trim_a_triangle_soup(srf, vert, vert_p, bd, norm, trim_curve_p_all[c], contour_all[c], mesh
#  ifdef DBG
			   , col
#  endif
	);
  
#endif // end of #ifdef SPLIT_TRIANGLES_AT_OUTER_CURVE





#ifdef SPLIT_TRIANGLES_AT_INNER_CURVES

#ifdef DBG
    puts("====================================================================================================");
    puts("Splitting triangles according to the second and further trimming curves.");
#endif

    //printf("Triangles added in 'make_trimmed_mesh': %d\n", int(mesh.size())/3);

    //
    // 081208: Adding a post-processing step in which we take into consideration further trimming curves. Now we
    //         have a mesh consisting of a number of triangles. We loop through all of them, and split them
    //         similarly to what was done with the first trimming curve above.
    //
    //         (For max code reuse, clarity and orthogonality, we should really replace the above, older code,
    //         with a pass similar to the one following below.)
    //
    //         First a version without the cached inside-results, since these are kind of messed up now, when we
    //         must free us from the regular grid... can be fixed later...
    //          
    //           081210: What the h... did I mean by this? Shouldn't using the cached values be
    //                   straightforward?!
    //
    // 081210: Note that the number of vertices are not reduced, even when it could be. Should be relatively
    //         easy to add a pass for doing such reduction later, if needed...
    //

#ifdef DBG
    printf("\n\n================================================== crv_set.size=%d\n\n\n", int(crv_set.size()));
#endif

    for (int c=1; c<int(crv_set.size()); c++)
      //for (int c=3; c<5; c++)
	{
#ifdef DBG
	  printf("----------==========########## trimming against inner curve %d (range=[%d, %d]) ##########==========----------\n", c, 1, int(crv_set.size())-1);
#endif
	  trim_a_triangle_soup(srf, vert, vert_p, bd, norm, trim_curve_p_all[c], contour_all[c], mesh,
#ifdef DBG
			       col,
#endif
			       true);
	}
      
#endif // end of #ifdef SPLIT_TRIANGLES_AT_INNER_CURVES




#if 0

#  if 1
  
    // 100210: What was this?!?! probably some debugging... Seems to have been used for some hardcoded inner
    //         trimming curves... (number 2 and 3..)
  
    trim_a_triangle_soup(srf, vert, vert_p, bd, norm, trim_curve_p_all[2], contour_all[2], mesh, col, true);
    {
      FILE *f=fopen("e1.m", "w");
      for (int i=0; i<int(trim_curve_p_all[2].size())-1; i++)
	{
	  const double &x0=trim_curve_p_all[2][contour_all[2][i]/3  ][0];
	  const double &y0=trim_curve_p_all[2][contour_all[2][i]/3  ][1];
	  const double &x1=trim_curve_p_all[2][contour_all[2][i]/3+1][0];
	  const double &y1=trim_curve_p_all[2][contour_all[2][i]/3+1][1];
	  fprintf(f, "line([%f; %f], [%f; %f], 'color', 'r', 'linewidth', 3);\n", x0, x1, y0, y1);
	}
      fclose(f);
    }
#  endif

    puts("####################################################################################################");
    puts("####################################################################################################");
    puts("####################################################################################################");
  
  
#  if 0
    trim_a_triangle_soup(srf, vert, vert_p, bd, norm, trim_curve_p_all[3], contour_all[3], mesh, col, true);
    {
      FILE *f=fopen("e2.m", "w");
      for (int i=0; i<int(trim_curve_p_all[3].size())-1; i++)
	{
	  const double &x0=trim_curve_p_all[3][contour_all[3][i]/3  ][0];
	  const double &y0=trim_curve_p_all[3][contour_all[3][i]/3  ][1];
	  const double &x1=trim_curve_p_all[3][contour_all[3][i]/3+1][0];
	  const double &y1=trim_curve_p_all[3][contour_all[3][i]/3+1][1];
	  fprintf(f, "line([%f; %f], [%f; %f], 'color', 'r', 'linewidth', 3);\n", x0, x1, y0, y1);
	}
      fclose(f);
    }
#  endif
  
#endif






#ifdef DBG

    //--------------------------------------------------------------------------------------------------------------
    //
    // 100210: Showing all remaining triangles. Colour coding:
    //          
    //           yellow: 
    //           green:
    //           red:
    //           cyan:
    //           magenta:
    //
    //--------------------------------------------------------------------------------------------------------------

    FILE *f=fopen("vis.m", "w");
    fprintf(f, "clf\n");
    for (int i=0; i<int(mesh.size())/3; i++)
      {
	const int u=0, v=1;
#if 0
	if ( (vert_p[mesh[3*i]][u] > -40) &&
	     (vert_p[mesh[3*i]][u] < -10) &&
	     (vert_p[mesh[3*i]][v] >  25) &&
	     (vert_p[mesh[3*i]][v] <  60)    )
#endif
	  {
	    fprintf(f, "patch([%f; %f; %f], [%f; %f; %f], '%c');\n",
		    vert_p[mesh[3*i]][u], vert_p[mesh[3*i+1]][u], vert_p[mesh[3*i+2]][u], 
		    vert_p[mesh[3*i]][v], vert_p[mesh[3*i+1]][v], vert_p[mesh[3*i+2]][v], col[i]);
	    fprintf(f, "line([%f; %f; %f; %f], [%f; %f; %f; %f], 'color', 'k');\n",
		    vert_p[mesh[3*i]][u], vert_p[mesh[3*i+1]][u], vert_p[mesh[3*i+2]][u], vert_p[mesh[3*i]][u], 
		    vert_p[mesh[3*i]][v], vert_p[mesh[3*i+1]][v], vert_p[mesh[3*i+2]][v], vert_p[mesh[3*i]][v]);
	  }
      }
    fprintf(f, "zoom on\n");
    fprintf(f, "grid on\n");
    // fprintf(f, "axis('equal');\n");
    fprintf(f, "xlabel('u');\n");
    fprintf(f, "ylabel('v');\n");
    fprintf(f, "title('"
	    "green   = post-corner-splitting\\newline"
	    "red=post-first-curve-quad-splitting\\newline"
	    "cyan    = post-outer-curve-triangle-splitting, first iteration\\newline"
	    "magenta = post-inner-curve-triangle-splitting, first iteration\\newline"
	    "yellow  = post-inner-curve-triangle-splitting\\newline"
	    "blue    = post-outer-curve-triangle-splitting');\n");
    // fprintf(f, "set(gca, 'ydir', 'reverse');\n"); // 100210: To get the same view as default 'goview'.
    fclose(f);
  
#endif
      



    //
    // 081210: Note that it does not make as much sense as before to return these now, if we have more than one
    //         trimming curve. But the caller is probably not using them... As it is, we return only for the
    //         outer trimming curve.
    //
    trim_curve_0   = trim_curve_all[0];
    trim_curve_p_0 = trim_curve_p_all[0];

    return;
  }






} // namespace Go
