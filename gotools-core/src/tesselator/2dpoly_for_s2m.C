#include "GoTools/tesselator/2dpoly_for_s2m.h"
using std::vector;
#include <functional>
#include <algorithm>

#include <stdio.h>






#ifdef DBG			// 100222: A construction for combining debugging facilities, readability and
				//         optimization in the release-version.
#  define DBG_FLAG ,dbg
#else
#  define DBG_FLAG
#endif






namespace Go
{





/** struct for sorting a list of 2D-values.
 */
  struct less_largest_x:public std::binary_function<short_list_short_list, short_list_short_list, bool>
  {
  private:
    const vector<int> *contour_ptr;
    const double *vertices_ptr;
    int segments;
    double z;
    const int t;
  
  public:
    less_largest_x(const vector<int> &c, const double * const v,
		   const double &aux_value=0.0,
		   const bool transposed=false):
      contour_ptr(&c), vertices_ptr(v), z(aux_value), t(transposed)
      {
	  segments=(int)contour_ptr->size();
      }
  
    bool operator()(const short_list_short_list &a, const short_list_short_list &b)
      {
	const int a_segment=a.first; // Index into the contour-array.
	double ax0, ax1;
	if (a_segment==-1)
	  ax0=ax1=z;
	else
	  {
	    const int &a_start=(*contour_ptr)[a_segment];
	    // Note that *(contour_ptr+a_segment) is wrong, it's a pointer to vector<int>, not a pointer to int!
	    //sg const int &a_end=(*contour_ptr)[(a_segment+1)%segments];
	    const int &a_end = (*contour_ptr)[a_segment+1==segments ? 0 : a_segment+1];
	    ax0=vertices_ptr[a_start+t];
	    ax1=vertices_ptr[a_end  +t];
	  }
	const int b_segment=b.first;
	double bx0, bx1;
	if (b_segment==-1)
	  bx0=bx1=z;
	else
	  {
	    const int &b_start=(*contour_ptr)[b_segment];
	    //sg	  const int &b_end=(*contour_ptr)[(b_segment+1)%segments];
	    const int &b_end = (*contour_ptr)[b_segment+1==segments ? 0 : b_segment+1];
	    bx0=vertices_ptr[b_start+t];
	    bx1=vertices_ptr[b_end  +t];
	  }

	//printf("%10.3f %10.3f %10.3f %10.3f %d\n", ax0, ax1, bx0, bx1,
	//(std::max(ax0, ax1)<std::max(bx0, bx1)));
	return (std::max(ax0, ax1)<std::max(bx0, bx1));
      }
  };






/** struct for sorting a list of 2D-values.
 */

  struct less_largest_y:public std::binary_function<short_list, short_list, bool>
  {
  private:
    const vector<int> *contour_ptr;
    const double *vertices_ptr;
    int segments;
    double z;
    const int t;
  
  public:
    less_largest_y(const vector<int> &c, const double * const v,
		   const double &aux_value=0.0,
		   const bool transposed=false):
      contour_ptr(&c), vertices_ptr(v), z(aux_value), t(!transposed)
      {
	  segments=(int)contour_ptr->size();
      }
  
    bool operator()(const short_list &a, const short_list &b)
      {
	const int a_segment=a.first;
	double ay0, ay1;
	if (a_segment==-1)
	  ay0=ay1=z;
	else
	  {
	    const int &a_start=(*contour_ptr)[a_segment];
	    //	  const int &a_end=(*contour_ptr)[(a_segment+1)%segments];
	    const int &a_end = (*contour_ptr)[a_segment+1==segments ? 0 : a_segment+1];
	    ay0=vertices_ptr[a_start+t];
	    ay1=vertices_ptr[a_end  +t];
	  }
	const int b_segment=b.first;
	double by0, by1;
	if (b_segment==-1)
	  by0=by1=z;
	else
	  {
	    const int &b_start=(*contour_ptr)[b_segment];
	    //sg const int &b_end=(*contour_ptr)[(b_segment+1)%segments];
	    const int &b_end = (*contour_ptr)[b_segment+1==segments ? 0 : b_segment+1];
	    by0=vertices_ptr[b_start+t];
	    by1=vertices_ptr[b_end  +t];
	  }

	return (std::max(ay0, ay1)<std::max(by0, by1));
      }
  };






//==============================================================================================================
//
// 090129: The 'transposed' parameter seems not to be in used in spline2mesh. Removing it for efficiency.
//
//==============================================================================================================

  struct less_smallest_y:public std::binary_function<short, short, bool>
  {
  private:
    const vector<int> *contour_ptr;
    const double *vertices_ptr;
    int segments;
    double z;
    //const int t;

  public:
    less_smallest_y(const vector<int> &c, const double * const v, const double &aux_value=0.0
		    //, const bool transposed=false
      ):
      contour_ptr(&c), vertices_ptr(v), z(aux_value)//, t(!transposed)
      {
	  segments=(int)contour_ptr->size();
      }
  
    inline bool operator()(const short &a_segment, const short &b_segment)
      {
#if 0
	double ay0, ay1;
	if (a_segment==-1)
	  ay0=ay1=z;
	else
	  {
	    const int &a_start=(*contour_ptr)[a_segment];
	    //sg const int &a_end=(*contour_ptr)[(a_segment+1)%segments];
	    const int &a_end = (*contour_ptr)[a_segment+1==segments ? 0 : a_segment+1];
	    ay0=vertices_ptr[a_start+t];
	    ay1=vertices_ptr[a_end  +t];
	  }
	double by0, by1;
	if (b_segment==-1)
	  by0=by1=z;
	else
	  {
	    const int &b_start=(*contour_ptr)[b_segment];
	    //sg const int &b_end=(*contour_ptr)[(b_segment+1)%segments];
	    const int &b_end = (*contour_ptr)[b_segment+1==segments ? 0 : b_segment+1];
	    by0=vertices_ptr[b_start+t];
	    by1=vertices_ptr[b_end  +t];
	  }

	return (std::min(ay0, ay1)<std::min(by0, by1));
#else
	double ay0, ay1;
	if (a_segment==-1)
	  ay0=ay1=z;
	else
	  {
	    const int &a_start=(*contour_ptr)[a_segment]+1;
	    const int &a_end = (*contour_ptr)[a_segment+1==segments ? 0 : a_segment+1]+1;
	    ay0=vertices_ptr[a_start];
	    ay1=vertices_ptr[a_end];
	  }
	double by0, by1;
	if (b_segment==-1)
	  by0=by1=z;
	else
	  {
	    const int &b_start=(*contour_ptr)[b_segment]+1;
	    const int &b_end = (*contour_ptr)[b_segment+1==segments ? 0 : b_segment+1]+1;
	    by0=vertices_ptr[b_start];
	    by1=vertices_ptr[b_end];
	  }

	return (std::min(ay0, ay1)<std::min(by0, by1));
#endif
      }
  };






//==============================================================================================================
//
// 081206: Is this used at all? Or replaced by 'point_inside_contour_for_s2m'?!
//
// 090115: Such a routine will always be a problem, because different applications will need different biases
//         for falling down with points either on the inside or outside for ambiguous cases...
//
//         Modifying the current functions for different such cases, or increasing "tolerance" for the
//         appropriate tests, is not easy to do without a full review of all the "sorted_segments"-stuff.
//
//         Trying another approach, with a separate "point_on_contour_corners"-routine. 
//
//         If needed, this should be easy to extend to "point_on_contour", but note that such a routine,
//         implemented in the obvious, naive way, would negate the speed advantages of the
//         "sorted_segments"-stuff.
//
// 090129: Not so sure any more about the usefulness of 'sorted_segments'... Note that it was designed for
//         something totally different. I'm not sure the assumptions upon which its usefulness was based,
//         apply any more...
//
// 100211: There is some confusion with regard to tolerances here (and probably elsewhere in related
//         functions.) The problem, at least here and now, seems to be a usage of a tolerance involving the
//         adding of this tolerance to a number before a comparison, rather than the subtraction of it
//         followed by a comparison to zero... I.e.,
//
//           a < b+eps
//
//         does not work for b large (-5.4) and eps small (1e-12) because b+eps==b in this case. What should
//         we do instead?
//
//           a-b < eps ?
//
//         Can we be sure that the compiler does not produce the same code for these two cases? Or should we
//         simply make sure eps is large enough for b+eps != b for all possible (!) cases of b?!
//
//         Going for a solution with relative tolerances instead, and also increasing the tolerance used where
//         the problem appeared.
//
// 100212: Contrary to what I believed, this routine is passed contours with duplicated corner points.
//         Adding a test and skipping zero-length segments.
//
//==============================================================================================================

  inline bool point_inside_contour(const double x0, const double y0,
				   const double * const vertices,
				   const vector<int> &contour
#ifdef DBG
				   , bool dbg /* =false */
#endif
    )
  {
#ifndef DBG
    const bool dbg=false; // This way, all tests on 'dbg' is removed compile-time when DBG is not defined.
#endif

    //
    // At last we enter the old loop for counting intersections.
    //
    // 090117: Note that the rest of the code is very similar to 'segment_contour_intersection_for_s2m'. In
    //         fact, we could have called that one. (Maybe we should have?) The difference is that now we know
    //         that the segment is horizontal, (the ray,) and hence the code here can be made (slightly) simpler
    //         and faster. One disadvantage is that the two functions (this and
    //         'segment_contour_intersection_for_s2m') must agree on how to handle degenerate cases and so
    //         on. Otherwise strange problems will occur.

    int crossings=0;

    // const double eps=1e-8; // 100210: Increasing from 1e-12 to 1e-10.
    
    const double abs_eps = 1e-8;	// 100223: Trying to convert to these
    // const double snap_eps = 1e-5;

    const int n=(int)contour.size();
  
    if (dbg)
      printf("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ n=%d\n", n);

    for (int i=0; i<n; i++)
      {
	const int j=(i+1)%n, pre_i=(i+n-1)%n;
	const double &preb=vertices[contour[pre_i]+1];
	const double &prea=vertices[contour[pre_i]];
	const double &a=vertices[contour[i]], &b=vertices[contour[i]+1];
	const double &c=vertices[contour[j]], &d=vertices[contour[j]+1];
      
	const double len_squared = (a-c)*(a-c) + (b-d)*(b-d);
	if ( (len_squared<abs_eps*abs_eps) && (dbg) )
	  printf("\n\n  HUH?! Segment with length %e, duplicate corner!\n\n\n", sqrt(len_squared));

// 	if (dbg)
// 	  printf("  i=%3d: a=%g b=%g c=%g d=%g y0=%g abs_eps=%g !horiz=%d y0-inrange=%d(%d)\n", 
// 		 i, a, b, c, d, y0, abs_eps, 
// 		 (fabs(b-d)>abs_eps),
// 		 ((b-abs_eps<=y0) && (y0<d-abs_eps)),
// 		 ((d+abs_eps<y0) && (y0<b+abs_eps)));
// 	if ( dbg && ((i==77) || (i==79)) )
// 	  {
// 	    printf("\n\n\nhold on; plot(%f, %f, 'md', 'markersize', 10, 'linewidth', 3); hold off;\n", x0, y0);
// 	    printf("line([%f; %f], [%f; %f], 'color', 'm', 'linewidth', 5);\n\n\n", a, c, b, d);
// 	  }
	
	//
	// 100223: Remember, crossings are counted when occuring at the *start* of segments, not the ends,
	//         *except* when the curve has a vertical turn in this "start", i.e., there is a corner
	//         "pointing" upward or downward.
	//

	const bool not_horizontal_segment = fabs(b-d) > abs_eps;
	const bool y0_in_range_1 = (b-y0<=abs_eps) && (abs_eps<d-y0);	// y0 in [b-eps, d-eps)
	const bool y0_in_range_2 = (abs_eps<y0-d) && (y0-b<=abs_eps);	// y0 in (d+eps, b+eps]
	const bool y0_in_range = y0_in_range_1 || y0_in_range_2;

	if ( not_horizontal_segment && y0_in_range )
	  {
	    const double t=(y0-b)/(d-b);	// Where on the segment does the ray cross? t=0 for b, and t=1 for d.
	  
	    if (dbg)
	      printf("  i=%3d: a=%g b=%g c=%g d=%g y0=%g abs_eps=%g !horiz=%d y0-inrange=%d %d\n", 
		     i, a, b, c, d, y0, abs_eps, 
		     not_horizontal_segment, y0_in_range_1, y0_in_range_2);
	    
	    const bool interior = ((t>=abs_eps) && (t<=(1.0-abs_eps)));
	    const bool start = (fabs(t)<=abs_eps);
	    const bool going_up = ((preb<b) && (b<d));
	    const bool going_down = ((preb>b) && (b>d));
	    const bool not_turning_vertically = going_up || going_down;
	    const bool actually_crossing = (interior || (start && not_turning_vertically));
	    const double intersection_x = a+t*(c-a);
	    if (dbg)
	      printf("\n\n          t=%g, interior: %d, start: %d, going_up: %d, going_down: %d, "
		     "actually_crossing: %d, inters_x: %g",
		     t, interior, start, going_up, going_down, actually_crossing, intersection_x);
	    if (dbg)
	      printf("\n          preb=%g, b=%g, d=%g    x-values: %g %g %g", preb, b, d, prea, a, c);
	    if ( actually_crossing && (intersection_x > x0+abs_eps) )
	      {
		crossings++;
		if (dbg) printf(", CROSSES ");
		if (dbg)
		  {
		    printf("\n\n\nhold on; plot(%f, %f, 'md', 'markersize', 10, 'linewidth', 3); hold off;\n", x0, y0);
		    printf("line([%f; %f], [%f; %f], 'color', 'm', 'linewidth', 5);\n\n\n", a, c, b, d);
		  }
	      }
	    else
	      if (dbg)
		{
		  if (actually_crossing)
		    printf("\n          intersection_x too small: %e vs. x0=%e, inters_x-x0=%e", 
			   intersection_x, x0, intersection_x-x0);
		  printf("\n\n\nhold on; plot(%f, %f, 'bd', 'markersize', 10, 'linewidth', 3); hold off;\n", x0, y0);
		  printf("line([%f; %f], [%f; %f], 'color', 'b', 'linewidth', 5);\n\n\n", a, c, b, d);
		}
	    if (dbg) printf("\n");
	  }
      
      }
    if (dbg) printf("  Crossings in total: %d, i.e., %s\n", crossings, ((crossings & 1)==1) ? "INSIDE" : "OUTSIDE");
  
    return ((crossings & 1)==1);
  }






//==============================================================================================================
//
// 090115: See also comment for 'point_inside_contour' above.
// 090117: This must be (re)checked before being used...
//
//==============================================================================================================

  bool point_on_contour_corner(const double x0, const double y0,
			       const double * const vertices, const vector<int> &contour)
  {
    const double eps=1e-13; // 090115: Used for zero-tests for distances in the parameter domain.
    //         Hmm... these should *really*, *really* be taken from some global variable
    //         or something
    for (int ilim=(int)contour.size(), i=0; i<ilim; i++)
      {
	const double &corner_u=vertices[contour[i]], &corner_v=vertices[contour[i]+1];
	const double dist_squared = (corner_u-x0)*(corner_u-x0) + (corner_v-y0)*(corner_v-y0);

	if (dist_squared<eps*eps)
	  return true;
      }
    return false;
  }






//==============================================================================================================
//
// 090115: See also comment for 'point_inside_contour' above.
//         For a first test, just doing brute force search. Note that this will negate the fast algorithm used
//         for the inside-tests. The current function should do something similar, if this turns out to be a
//         good idea.
// 090129: This can probably be done in point_inside_contour at very little extra cost...
// 100213: Seems to be a problem with this routine, for some special cases... investigating...
//
//==============================================================================================================

  inline bool point_on_contour(const double x0, const double y0,
			       const double * const vertices, const vector<int> &contour
#ifdef DBG
			       , const bool dbg=false
#endif
    )
  {
#ifndef DBG
    const bool dbg=false; // This way, all tests on 'dbg' is removed compile-time when DBG is not defined.
#endif

    //const double eps=1e-12;	// 090115: Used for zero-tests for degeneracy and parallellity tests.
    // 090203: Here it is also used for testing if a point is on a curve segment.
    //const double eps=1e-8; 	// 090203: Need 1e-8 for proper meshing of 'bin_p1_3.g2'.
    //const double eps=1e-6; 	// 100213: See comments below.
    const double eps=1e-5; 	// 100214: Reverting to 1e-8, think the need to have 1e-6 is really another problem
				// 100218: Keeping 1e-5 after discussion with Vibeke.

    const double tau=1e-12; 	// 090115: Used for zero-tests for distances in the parameter domain.
  				// 100214: Changed from 1e-14 to 1e-12, not needed right now, but 1e-14 seems small.

    const int n=(int)contour.size();
    if (dbg) printf("  point_on_contour for (%f, %f), %d segments.\n", x0, y0, n);
    for (int i=0; i<n; i++)
      {
	const int j=(i+1)%n;
	const double &x2=vertices[contour[i]], &y2=vertices[contour[i]+1];
	const double &x3=vertices[contour[j]], &y3=vertices[contour[j]+1];
	const double contour_segment_length_squared = (x3-x2)*(x3-x2) + (y3-y2)*(y3-y2);

	if (contour_segment_length_squared>tau*tau)
	  {
	    // We project (x0, y0) onto the segment, and measure the distance.
	    // The projection will be (x2, y2) + s*[(x3, y3)-(x2, y2)].
	    const double s = ( (x0-x2)*(x3-x2) + (y0-y2)*(y3-y2) ) / contour_segment_length_squared;

	    // 100213: Hmm... Checking for s in (-tau, 1+tau) is no good, if 1+tau==1. Better to use (-eps,
	    //         1+eps) when eps is larger than tau, but this does not solve the problem! Which is that
	    //         's' is not a measure of distance. Should rather use
	    //         s*sqrt(contour_segment_length_squared)! 
	    //
	    // 100214: But, then, remember that the interval should be [0,
	    //         sqrt(contour_segment_length_squared)].
	    //
	    // 100213: Hmm again. Seems that 1e-8 (and 1e-7!) for some reason is too restrictive. (How can this
	    //         be?!)  Changing to 1e-6. 
	    //
	    // 100214: Suspect that the reason for 1e-8 being to large is that something else is really wrong
	    //         with the contour of the test case. (Folding back onto itself or something strange?!)

	    const double s2 = s * sqrt(contour_segment_length_squared);

	    if ( (dbg) && (0) )
	      {
		printf("    %3d: (x2, y2) = (%f, %f), (x3, y3) = (%f, %f)\n", i, x2, y2, x3, y3);
		printf("line([%f; %f], [%f; %f], 'color', 'y', 'linewidth', 5);\n", x2, x3, y2, y3);
		printf("         s=%f, s2=%f\n", s, s2);
	      }

	    // if ((s>-tau) && (s<1.0+tau))
	    // if ((s>=-eps) && (s<=1.0+eps))
	    // if ((s2>=-eps) && (s2<=1.0+eps))
	    if ( (s2>=-eps) && (s2<=sqrt(contour_segment_length_squared)+eps) )
	      {
		const double x = x2 + s*(x3-x2), y = y2 + s*(y3-y2);
		const double dist_squared = (x-x0)*(x-x0) + (y-y0)*(y-y0);
		if (dbg) printf("\t\t\t\t\t\t\t\t\t\tdist to segment = %e\n", sqrt(dist_squared));

		// 100214: That the tolerance required in the next test is as large as 1e-4 really shows that
		//         the trimming curve is crappy, so by using 1-e4 we are really just hiding a symptom...

		if (dist_squared<eps*eps)
		  //if (dist_squared<1e-6*1e-6)
		  //if (dist_squared<1e-4*1e-4)
		  {
		    if (dbg) printf("TRUE\n");
		    return true;
		  }
	      }
	  }
	else
	  {
	    //
	    // 090117: The contour segment is degenerate, and will be treated as a point.
	    //
	    // 100214: We should not get here, duplicate points should have been removed, but it's no harm in
	    //         keeping the branch I guess...
	    //
	    const double mx = 0.5*(x2+x3), my = 0.5*(y2+y3);
	    const double dist_squared = (x0-mx)*(x0-mx) + (y0-my)*(y0-my);
	    if (dist_squared<eps*eps)
	      {
		if (dbg) printf("TRUE\n");
		return true;
	      }
	  }
      }

    if (dbg) printf("FALSE\n");
    return false;
  }






  //==============================================================================================================
  //
  // 100219: See comments with corresponding dates in 'segment_contour_intersection_for_s2m()'.
  //
  //==============================================================================================================

  bool point_on_segment(const double x0, const double y0,
			const double x2, const double y2,
			const double x3, const double y3,
			const double zero_eps,
			const double snap_eps
#ifdef DBG
			, const bool dbg = false
#endif
    )
  {
#ifndef DBG
    const bool dbg = false; // This way, all tests on 'dbg' is removed compile-time when DBG is not defined.
#endif

    const double t = ( (x0-x2)*(x3-x2) + (y0-y2)*(y3-y2) ) / ( (x3-x2)*(x3-x2) + (y3-y2)*(y3-y2) );
    const double x = x2 + t*(x3-x2) - x0;
    const double y = y2 + t*(y3-y2) - y0;
    const double dist = sqrt( x*x+y*y );
    if ( (dist<snap_eps) && (t>=-zero_eps) && (t-zero_eps<=1.0) )
      {
	if (dbg)
	  {
	    printf("    END TESTS 2: t=%f, dist=%f\n", t, dist);
	    // printf("\nline([%f; %f], [%f; %f], 'color', 'k', 'linewidth', 6);\n", x0, x1, y0, y1);
	    printf("line([%f; %f], [%f; %f], 'color', 'y', 'linewidth', 6);\n", x2, x3, y2, y3);
	  }
	return true;
      }

    return false;
  }






  //==============================================================================================================
  //
  // 081208: Comment added today.
  //
  //         The contour consists of n segments, so we are finding an intersection between two straight line
  //         pieces. The corners of the contour polygon are stored in 'vertices'. 
  //         The 'x0', 'y0', 'x1' and 'y1' seem to be points in the parameter domain of the surface/curve.
  //
  // 081209: There seems to be a problem with segments of length 0 in the contour... why haven't this happened
  //         before? Or is it just a symptom of something else being the matter?
  //
  // 090115: The return value is true if an intersection between the segment and contour is found.
  //
  // 090117: Note that the contour may (of course) intersect the segment more than once. Just one intersection
  //         is returned, and it is pretty indetermined which one it will be.
  //
  // 090117: If we modify the behaviour of "insideness-tests" in spline2mesh, this must be mirrored
  //         here. Otherwise, 'split_quad' will not find the intersections it expect, and fail!
  //         Or - we could do it in the caller, e.g., in 'split_quad'.
  //
  // 090203: There is a fundamental problem with the whole approach. In some cases, a quad (or it might as well
  //         be any polygon, e.g., a triangle) has an edge that will intersect the trimming curve ('contour') in
  //         more than one place. For example, a contour which follows an iso-parameter line, say a border of
  //         the surface, and then turns sharpely back, making out a wedge. Somewhere near the tip of this edge,
  //         this situation will always occur, independent of tesselation level etc. If the segment of a quad
  //         (triangle) normal to the border is determined to intersect the contour in the end, i.e., at the
  //         border of the surface, the wedge tip will be cut away wrongly...
  //         How can this be solved??
  //         Not sure if this is a universal fix, but for now, we will favor inner intersections before
  //         "end-of-segment" intersections...
  //         This also means that we will not be able to return as soon as we find an intersection, slowing
  //         things down somewhat. Considerably?
  //         A compromise: When first intersection found is of type "end", continue until a "better" is found.
  //
  // 100218: This situation has again surfaced. Why does it not work, if it was fixed before?
  //         Making the FIX090203 permanent.
  //         Increasing tolerances.
  //
  // 100219: Hmm... I do not entirely like the way things have turned here. By favoring inner intersections,
  //         we are in a way handling the case where there are many intersections on an edge, and choosing one
  //         of them, as opposed to finding *all* and also treating *all*...
  //
  // 100219: There was a "bug" here: Originally the routine was not supposed to classify any intersection on
  //         parallel segments as an intersection. However, this is needed when called from 'split_triangle'
  //         (as opposed to from a ray-contour problem of a point-inside-contour context,) where we want
  //         intersections to be found when the ends of the triangle edge lie *on* (or, actually, 1e-5 or so
  //         from) such a parallel (to the edge) contour segment! Fixing this now, and adding a flag, so that
  //         it will not be the default behaviour. ('snap_ends', see also comments in the code below.)
  //
  //         100225: Should 'snap_ends' be turned to default=true, and then removed altogether when proving
  //                 itself to be a good idea?
  //
  // 100222: The question is, then, who (and when) should turn this on? Following the conservative approach,
  //         which is to make it default to false, and turn it (minimally) on in order to solve the current
  //         (trim3.g2) case...
  //
  // 100222: One thing, which is not entirely satisfactory with the whole implementation:
  //
  //           The "insidedness/onness" for points vs. contours used in routines for determination of which,
  //           and how, triangles are to be clipped, is not the same, nor equivalent, to those used during the
  //           actual clipping. The reason for this is, e.g., that these two operations due to their
  //           implementation/organization require slightly different definitions of insidedness. This could
  //           (if true, then it should) maybe be fixed.
  //
  //==============================================================================================================

  bool segment_contour_intersection_for_s2m(const double x0, const double y0,
					    const double x1, const double y1,
					    const double * const vertices,
					    const vector<int> &contour,
					    double &x, double &y, // Resulting intersection
					    double &s, // param for the segment containing the inters, [0, 1].
					    const bool snap_ends /* = false */
#ifdef DBG
					    , const bool dbg /* = false */
#endif
    )
  {
#ifndef DBG
    const bool dbg=false; // This way, all tests on 'dbg' is removed compile-time when DBG is not defined.
#endif

    if (dbg) puts("\n\n  # Entering segment_contour_intersection_for_s2m ########################################\n");

    //
    // Now we search for the (actually, just "a") piece of the contour which intersects the segment. Maybe it is
    // faster to just do this without the previous test, if we can use the 'sorted_segments' information.
    //

    int i;
    //const double tau=1e-14; // 090115: Used for zero-tests for degeneracy and parallellity tests.
    const double tau=1e-12; // 100218:

    const double parallellity_eps = 1e-10; // 100219

    const double eps=1e-13; // 090115: Used for zero-tests for distances in the parameter domain.
    //const double eps2=1e-8; // 090203: This is used to determine whether or not an intersection is at one of the ends
    const double eps2=1e-5; // 100218:
    //         of an edge. (In which case we continue to look for an "interior" intersection,
    //         which will be favoured. Is this always sensible? Or is it possible to make a
    //         case in which the opposite behaviour should be preferred?
    // 090204: We need this as large as 1e-8 for bin_p1_3.g2 to be meshed properly.
    bool intersection_found_at_the_end_of_the_segment = false;
    double x_saved=1e99, y_saved=1e99, s_saved=1e99; // Initial values to help spot unforeseen problems...
    const int n=(int)contour.size();
    for (i=0; i<n; i++)
      {
	const int ii=i, jj=(ii+1)%n;
	const double &x2=vertices[contour[ii]], &y2=vertices[contour[ii]+1];
	const double &x3=vertices[contour[jj]], &y3=vertices[contour[jj]+1];
	double t;

#if 0
	if (dbg)
	  {
	    const double d02 = sqrt( (x2-x0)*(x2-x0) + (y2-y0)*(y2-y0) );
	    const double d03 = sqrt( (x3-x0)*(x3-x0) + (y3-y0)*(y3-y0) );
	    const double d12 = sqrt( (x2-x1)*(x2-x1) + (y2-y1)*(y2-y1) );
	    const double d13 = sqrt( (x3-x1)*(x3-x1) + (y3-y1)*(y3-y1) );
	    
	    const double eps3 = 1e-3;

	    if ( (d02<eps3) || (d03<eps3) || (d12<eps3) || (d13<eps3) || (i==111169) )
	      {
		printf("    i=%4d: END TESTS: d02=%f, d03=%f, d12=%f, d13=%f\n", i, d02, d03, d12, d13);
		printf("\nline([%f; %f], [%f; %f], 'color', 'k', 'linewidth', 6);\n", x0, x1, y0, y1);
		printf("line([%f; %f], [%f; %f], 'color', 'y', 'linewidth', 6);\n", x2, x3, y2, y3);
	      }

	    // testing if ends of first segment is in the other.

	    const double t0 = ( (x0-x2)*(x3-x2) + (y0-y2)*(y3-y2) ) / ( (x3-x2)*(x3-x2) + (y3-y2)*(y3-y2) );
	    double p = x2+t0*(x3-x2) - x0;
	    double q = y2+t0*(y3-y2) - y0;
	    const double dist0 = sqrt(p*p+q*q);
	    //if ( (dist0<1e-5) && (t0>=-eps) && (t0<=1.0+eps) )
	    if (i==69)
	      {
		printf("    i=%4d: END TESTS 2: t0=%f, dist=%f\n", i, t0, dist0);
		printf("\nline([%f; %f], [%f; %f], 'color', 'k', 'linewidth', 6);\n", x0, x1, y0, y1);
		printf("line([%f; %f], [%f; %f], 'color', 'y', 'linewidth', 6);\n", x2, x3, y2, y3);
	      }

	    p = x2+t0*(x3-x2) - x1;
	    q = y2+t0*(y3-y2) - y1;
	    const double dist1 = sqrt(p*p+q*q);
	    if (i==69)
	    //if ( (dist1<1e-5) && (t0>=-eps) && (t0<=1.0+eps) )
	      {
		printf("    i=%4d: END TESTS 3: t0=%f, dist=%f\n", i, t0, dist0);
		printf("\nline([%f; %f], [%f; %f], 'color', 'k', 'linewidth', 6);\n", x0, x1, y0, y1);
		printf("line([%f; %f], [%f; %f], 'color', 'y', 'linewidth', 6);\n", x2, x3, y2, y3);
	      }

	  }
#endif
      
	// 081209: Adding this test...
	const double contour_segment_length_squared = (x2-x3)*(x2-x3)+(y2-y3)*(y2-y3);
	// if (dbg2)	printf("i=%3d cont segm len=%f ", i, sqrt(contour_segment_length_squared));

	if (contour_segment_length_squared>tau*tau)
	  {

	    if (snap_ends)
	      {
		//
		// 100219: We check if the ends of the segment (the one given as input to the function) are
		//         close to the contour piece in question. If 'snap_ends' is true, we could have that
		//         situation even if the pieces are (close to) parallel. We do this testing now,
		//         before the parallellity-test, since that test is too strict to determine if this
		//         'snap_ends'-test should be performed. (It is too strict because it is meant to trap
		//         division with zero before they happen, mainly.)
		//
		//         Note also that if one of the ends is on the contour, this is an "end intersection",
		//         and this 'intersection_found_at_the_end_of_the_segment' is set and the search
		//         continued, in case there is also some inner intersection.
		//

		if (point_on_segment(x0, y0, x2, y2, x3, y3, eps, eps2)) // , prøve med true her og nedenfor, se om det alltid blir bra da.))
		  {
		    intersection_found_at_the_end_of_the_segment = true;
		    x_saved = x0;
		    y_saved = y0;
		    s_saved = 0.0;
		  }
		else
		  if (point_on_segment(x1, y1, x2, y2, x3, y3, eps, eps2 DBG_FLAG))
		    {
		      intersection_found_at_the_end_of_the_segment = true;
		      x_saved = x1;
		      y_saved = y1;
		      s_saved = 1.0;
		    }
	      }
	    
	    // 090115: The next one should be zero for contour segment and segment parallel. If so, no intersection...
	    const double parallellity = (x3-x2)*(y1-y0)-(y3-y2)*(x1-x0);
	    if ( fabs(parallellity) > parallellity_eps ) // i.e., not parallel...
	      {
		//
		// 100219: Let the (triangle) segment has ends A=(x0, y0) and B=(x1, y1), and the contour
		//         piece has ends C=(x2, y2) and D=(x3, y3).  With intersection S = A+s(B-A) =
		//         C+t(D-C), we get the following:
		//
		t = ((x0-x2)*(y1-y0)-(y0-y2)*(x1-x0)) / parallellity;
		if ((t>=-eps) && (t<=1.0+eps))
		  {
		    if (fabs(x1-x0)>fabs(y1-y0))
		      s = ((x2-x0)+t*(x3-x2))/(x1-x0);
		    else
		      s = ((y2-y0)+t*(y3-y2))/(y1-y0);
		    // if (dbg2) printf("s=%f ", s);
		    if ((s>=-eps) && (s<=1.0+eps))
		      {
			if (dbg)
			  {
			    const double ss = ((x0-x2)*(y3-y2)-(y0-y2)*(x3-x2)) / parallellity;
			    printf("    i=%4d: s=%f (%f), t=%f, parallellity=%f\n", i, s, ss, t, parallellity);
			    printf("\nline([%f; %f], [%f; %f], 'color', 'k', 'linewidth', 6);\n", x0, x1, y0, y1);
			    printf("line([%f; %f], [%f; %f], 'color', 'y', 'linewidth', 6);\n", x2, x3, y2, y3);
			    const double px = x0 + s*(x1-x0), py = y0 + s*(y1-y0);
			    const double qx = x2 + t*(x3-x2), qy = y2 + t*(y3-y2);
			    printf("hold on; plot(%f, %f, 'gs', 'markersize', 15, 'linewidth', 2);\n", px, py);
			    printf("plot(%f, %f, 'ks', 'markersize', 15, 'linewidth', 2); hold off;\n\n", qx, qy);
			  }
			s = std::max(std::min(s, 1.0), 0.0);
			x = x2+t*(x3-x2);
			y = y2+t*(y3-y2);
			//printf("***** s=%f t=%f x=%f y=%f\n", s, t, x, y);
			//printf("returning true from %s:%d...\n", __FILE__, __LINE__);
			if ((fabs(s)<eps2) || (fabs(s-1.0)<eps2))
			  {
			    intersection_found_at_the_end_of_the_segment=true;
			    if (dbg) printf("Intersection is at and end.\n");
			    x_saved=x;
			    y_saved=y;
			    s_saved=s;
			  }
			else
			  {
			    if (dbg) puts("  # Returning true from segment_contour_intersection_for_s2m # (a) #####\n");
			    return true;
			  }
		      }
		    else
		      if (dbg) printf("    i=%4d: NO INTERS, s OUT OF RANGE s=%f, t=%f, p=%g\n", i, s, t, parallellity);
		  }
		else
		  if (dbg) printf("    i=%4d: no inters, t out of range s=%f, t=%f, p=%g\n", i, s, t, parallellity);
	      }
	    else
	      if (dbg)
		{
		  printf("    i=%4d: PARALLEL parallellity=%e\n", i, parallellity);
		  // The egde of the triangle ("segment") and the piece of the contour:
		  printf("\nline([%f; %f], [%f; %f], 'color', 'k', 'linewidth', 6);\n", x0, x1, y0, y1);
		  printf("line([%f; %f], [%f; %f], 'color', 'y', 'linewidth', 6);\n\n", x2, x3, y2, y3);
		}
	  }
	else
	  {
	    if (dbg) printf("    Hmm... degen contour segment... should that still be possible?!\n");
	    // 090117: The contour segment is degenerate, and will be treated as a point.
	    const double mx = 0.5*(x2+x3), my = 0.5*(y2+y3);

	    // The projection will be (x0, y0) + s*[(x1, y1)-(x0, y0)].
	    s = ( (mx-x0)*(x1-x0)+(my-y0)*(y1-y0) ) / ( (x1-x0)*(x1-x0)+(y1-y0)*(y1-y0) );
	    // if (dbg2) printf("s=%f ", s);

	    if ((s>-eps) && (s<1.0+eps))
	      {
		x = x0 + s*(x1-x0);
		y = y0 + s*(y1-y0);
	      
		const double dist_squared = (x-mx)*(x-mx) + (y-my)*(y-my);
		if (dist_squared<eps*eps)
		  {
		    intersection_found_at_the_end_of_the_segment=true;
		    x_saved=x;
		    y_saved=y;
		    s_saved=s;
		  }
	      }
	  }
	// if (dbg2) printf("\n");
      }

    if (intersection_found_at_the_end_of_the_segment)
      {
	x=x_saved;
	y=y_saved;
	s=s_saved;
	if (dbg) puts("  # Returning true from segment_contour_intersection_for_s2m # (b) #######################\n");
	return true;
      }
  
    // 090205: This message is not longer as relevant (or correct!) Now (also before?) it is perfectly
    //         acceptable that this returns 'false' for the "second intersection, first try". The problem is
    //         that it should not happen in the "second try" done by 'split_triangle'... But the proper place to
    //         catch that problem is most likely in 'split_triangle' itself, so I'll disable this old message...
//   MESSAGE("\nCouldn't find proper intersection between segment\n"
// 	  "and contour, while the segment's ends\n"
// 	  "are determined to be on opposite sides of the contour!!!\n"
// 	  "If you see this, the results are very likely going to look like crap, unfortunately.\n"
// 	  "It is currently unknown how this could happen.\n");
//   printf("Segment length in the parameter domain: %g\n\n\n", sqrt((x0-x1)*(x0-x1) + (y0-y1)*(y0-y1)));

    if (dbg) puts("  # Returning false from segment_contour_intersection_for_s2m # (c) ##########################\n");

    return false;
  }






//==============================================================================================================
//
//==============================================================================================================

  vector< short_list_short_list > sort_2dpoly_segments(const double * const vertices,
						       const vector<int> &contour,
						       const bool transposed /* =false */)
  {
    vector< short_list_short_list > sorted_segments;
    const int n=(int)contour.size();

    const int debug=0;
  
    // 020615: Introducing the 'transposed' for making a version in which x and y swap roles.

    //
    // We cannot choose, it must contain indices into the 'contour' vector, if we were to copy that vector's
    // indices directly into 'vertices', we would lose information about one of the ends of all segments, as
    // soon as we started reordering the segments.  Hence, initializing with 'contour[i]' was wrong. It must be
    // 'i'.
    sorted_segments.resize(n);
    for (int i=0; i<n; i++)
	sorted_segments[i].first=(short)i;

    // Now we sort this "outer" list wrt. largest x-value of the segments.
    {
      less_largest_x compare_functor(contour, vertices, 0.0, transposed);
      std::sort(sorted_segments.begin(), sorted_segments.end(), compare_functor);
    }
  
    if ((debug) && (!transposed))
      {
	//short_list_short_list *first_cut;
	vector<short_list_short_list>::iterator first_cut;
      
	printf("\n");
	for (first_cut=sorted_segments.begin(); first_cut<sorted_segments.end(); first_cut++)
	  {
	    const int a=contour[first_cut->first], b=contour[(first_cut->first+1)%n];
	    printf("%3d: %5d-%5d: (%10.3f, %10.3f, %10.3f)-"
		   "(%10.3f, %10.3f, %10.3f) maxx: %10.3f\n",
		   int(first_cut-sorted_segments.begin()),
		   a, b,
		   vertices[a+0], vertices[a+1], vertices[a+2],
		   vertices[b+0], vertices[b+1], vertices[b+2],
		   std::max(vertices[a+0], vertices[b+0]));
	  }
      }
  
    // Next, we fill inn all the lists of each pair of this "outer" list, i.e., the "middle" lists.
    if (transposed)
      {
	THROW("Huh?! This should not be in use, I think... J.O. 090129");
	exit(0);
      }
    less_smallest_y lsy_comp(contour, vertices, 0.0); // , transposed);
    less_largest_y lly_comp(contour, vertices, 0.0); // , transposed);
  
    for (int i=0; i<n; i++)
      {
	sorted_segments[i].second=new vector<short_list>;
	sorted_segments[i].second->resize(n-i);
	for (int j=0; j<n-i; j++)
	  (*sorted_segments[i].second)[j].first=sorted_segments[i+j].first;
      
	// Now we sort this "middle" list wrt. largest y-value of the segments.
	std::sort(sorted_segments[i].second->begin(), sorted_segments[i].second->end(), lly_comp);
      
	if ((debug) && (!transposed))
	  {
	    vector<short_list>::iterator second_cut;
	  
	    printf("\n");
	    for (second_cut=sorted_segments[i].second->begin(); second_cut<sorted_segments[i].second->end(); second_cut++)
	      {
		const int a=contour[second_cut->first], b=contour[(second_cut->first+1)%n];
		printf("     %3d: %5d-%5d: (%10.3f, %10.3f, %10.3f)-"
		       "(%10.3f, %10.3f, %10.3f) maxy: %10.3f\n",
		       int(second_cut-sorted_segments[i].second->begin()),
		       a, b,
		       vertices[a+0], vertices[a+1], vertices[a+2],
		       vertices[b+0], vertices[b+1], vertices[b+2],
		       std::max(vertices[a+1], vertices[b+1]));
	      }
	  }
      
	// Next, we fill inn all the lists of each pair of this "middle" list, i.e., the "inner" lists.
	for (int j=0; j<n-i; j++)
	  {
	    (*sorted_segments[i].second)[j].second=new vector<short>;
	    (*sorted_segments[i].second)[j].second->resize(n-i-j);
	    for (int k=0; k<n-i-j; k++)
	      (*(*sorted_segments[i].second)[j].second)[k] = (*sorted_segments[i].second)[j+k].first;
	  
	    // Now we sort this "inner" list wrt. smallest y-value of the segments. Default 'compare' should be
	    // 'less'.
	    //
	    // 020614: Don't know why such a default comparison functor is mentioned. How could it ever work
	    //         with these compounded lists with the "strange" comparisons that are to be done?!
	    // 020615: Probably I thought that as the inner list's entries consist of just an integer, there is
	    //         no need for the complex compare_functor, but then I forgot that it's not the integers (or
	    //         shorts) that are to be used as key in the sorting, but that which they point to!
	    std::sort((*sorted_segments[i].second)[j].second->begin(),
		      (*sorted_segments[i].second)[j].second->end(),
		      lsy_comp);
	  
	    if ((debug) && (!transposed))
	      {
		vector<short>::iterator third_cut;
	      
		printf("\n");
		for (third_cut=(*sorted_segments[i].second)[j].second->begin();
		     third_cut<(*sorted_segments[i].second)[j].second->end();
		     third_cut++)
		  {
		    const int a=contour[*third_cut], b=contour[(*third_cut+1)%n];
		    printf("          %3d: %5d-%5d: (%10.3f, %10.3f, %10.3f)-"
			   "(%10.3f, %10.3f, %10.3f) miny: %10.3f\n",
			   int(third_cut-
			       (*sorted_segments[i].second)[j].second->begin()),
			   a, b,
			   vertices[a+0], vertices[a+1], vertices[a+2],
			   vertices[b+0], vertices[b+1], vertices[b+2],
			   std::min(vertices[a+1], vertices[b+1]));
		  }
	      }
	  
	  }
      }
  
//    if (transposed)
//      {
//        printf("\n");
//        print_sorted_segments2(sorted_segments, contour, vertices);
//      }
  
    return sorted_segments;
  }






//==============================================================================================================
//
// Function for testing of insidedness of point based on old piecewise linear contours.
//
// Note: The z-coo. is ignored by that routine. The curve is specified in the parameter domain.
//
// 081202: So what is the difference between 'trim_curve_p' and 'contour'?!
//         Ok, 'contour' is a list of indices into 'trim_curve_p', it seems.
//         (Suspecting that for usage here, 'contour' is simply {0, ..., n} for the proper 'n'...)
//         (Actually {0, 3, ..., 3*n}.)
//
// 100211: Removing traces of "sorted_segments". Also, making the 'contour' a reference. Why was it not?!?!
//
//==============================================================================================================

  int is_inside(const vector< Vector3D > &trim_curve_p,
		const vector<int> &contour,
		const double u, const double v
#ifdef DBG
		, const bool dbg /* =false */
#endif
    )
  {
#ifndef DBG
    // const bool dbg=false; // This way, all tests on 'dbg' is removed compile-time when DBG is not defined.
#endif


    // test 100212

    if (contour.size()==0)
      return 0;
  
#if 0 // 021012

    const double * const vertices = &trim_curve_p[0][0];

  
    vector<int> cont;
    cont.push_back(contour[0]);
    for (int i=0; i<int(contour.size()); i++)
      {
	const int j = (i+1)%int(contour.size());
	const double &a=vertices[contour[i]], &b=vertices[contour[i]+1];
	const double &c=vertices[contour[j]], &d=vertices[contour[j]+1];
	const double len_squared = (a-c)*(a-c) + (b-d)*(b-d);
	const double eps = 1e-8;
	if ( len_squared < eps*eps )
	  {
	    if (dbg)
	      printf("is_inside preprocessing: HUH?! Segment with length %e, duplicate corner! REMOVING IT\n",
		     sqrt(len_squared));
	  }
	else
	  cont.push_back(contour[i]);
      }


    for (int s=2; s<=10; s++)
      {
	vector<int> cont2;
	cont2.push_back(cont[0]);
	for (int i=0; i<int(cont.size()); i++)
	  {
	    const int j = (i+1)%int(cont.size());
	    const double &a=vertices[cont[i]], &b=vertices[cont[i]+1];
	    const double &c=vertices[cont[j]], &d=vertices[cont[j]+1];
	    const double len_squared = (a-c)*(a-c) + (b-d)*(b-d);
	    const double eps = 1e-8;
	    if ( len_squared < eps*eps )
	      {
		if (dbg)
		  printf("is_inside prep STAGE %d : HUH?! Segment %d with length %e, duplicate corner! REMOVING IT\n",
			 s, i, sqrt(len_squared));
	      }
	    else
	      cont2.push_back(cont[i]);
	  }
	cont=cont2;
      }
  
    return point_inside_contour(u, v, &trim_curve_p[0][0], cont, dbg);
#endif
  



    // 100212: The const_cast not needed any more.
    //return point_inside_contour(u, v, &trim_curve_p[0][0], const_cast<vector<int> &>(contour), dbg);
    return point_inside_contour(u, v, &trim_curve_p[0][0], contour DBG_FLAG);
  }






//==============================================================================================================
//
// 090115:
// 090117: This should be (re)checked before being used... Not currently in use.
//
//==============================================================================================================

  bool is_on_corner(const vector< Vector3D > &trim_curve_p, const vector<int> &contour, const double u, const double v)
  {
    return point_on_contour_corner(u, v, &trim_curve_p[0][0], const_cast<vector<int> &>(contour));
  }






//==============================================================================================================
//
// 090117:
//
//==============================================================================================================

  int is_on_contour(const vector< Vector3D > &trim_curve_p, const vector<int> &contour, const double u, const double v
#ifdef DBG
		    , const bool dbg /* = false */
#endif
    )
  {
    return point_on_contour(u, v, &trim_curve_p[0][0], const_cast<vector<int> &>(contour) DBG_FLAG);
  }






  //==============================================================================================================
  //
  // 090204:
  // 100218: Eliminating unneeded test. Increasing the tolerance.
  //
  //==============================================================================================================
  
  bool degenerate_triangle(const Vector2D &c1, const Vector2D &c2, const Vector2D &c3
#ifdef DBG
			   , const bool dbg /* = false */
#endif
    )
  {
#ifndef DBG
    const bool dbg = false;
#endif

    // const double eps=1e-14;
    const double eps=1e-9;
    const double eps_squared=eps*eps;
  
    Vector2D e=c2-c1, f=c3-c1;
    const double e_len_squared=e*e, f_len_squared=f*f, g_len_squared=(c2-c3)*(c2-c3);

    if (dbg) printf("      degen_tri: edge lengths: %f %f %f\n", sqrt(e*e), sqrt(f*f), sqrt(g_len_squared));

    if ( (e_len_squared<eps_squared) || (f_len_squared<eps_squared) || (g_len_squared<eps_squared) )
      return true;

    // Ok, even if all lengths are non-zero, it may still be degenerate. Then the area will be zero.
    // 
    // 100219: Hmm... come to think of it, this test should be the only one needed...  But on the other hand,
    //         if we use the division below, we must still test on the size of the denominator, so we might as
    //         well keep it as it is. (Or maybe exchange the whole thing for an area-test?!)

    const double ef=(e*f)/sqrt(e_len_squared*f_len_squared);
    
    if (dbg) printf("      degen_tri: ef=%e, fabs(ef)-1.0=%e\n", ef, fabs(ef)-1.0);
    
    if (fabs(fabs(ef)-1.0) < eps)
      return true;
    
    return false;
  }






} // namespace Go
