//===========================================================================
//                                                                           
// File: GapRemovalVolume.C                                                  
//                                                                           
// Created: Wed Nov  2 16:12:41 2011                                         
//                                                                           
// Author: Sverre Briseid <Sverre.Briseid@sintef.no>
//                                                                           
// Revision: $Id$
//                                                                           
// Description:
//                                                                           
//===========================================================================

#include "GoTools/trivariate/GapRemovalVolume.h"
#include <assert.h>


namespace Go
{

//===========================================================================
void
GapRemoval::removeGapSpline(shared_ptr<SplineVolume>& vol1, 
			    shared_ptr<SurfaceOnVolume>& bd_sf1,
			    double sf1_start1, double sf1_end1,
			    double sf1_start2, double sf1_end2,
			    shared_ptr<SplineVolume>& vol2, 
			    shared_ptr<SurfaceOnVolume>& bd_sf2,
			    double sf2_start1, double sf2_end1,
			    double sf2_start2, double sf2_end2,
			    Point vertex_ll, Point vertex_ur,
			    double epsge, int orientation)
//===========================================================================
{
    ;//MESSAGE("removeGapSpline() under construction.");

#if 0

  // @@sbr201208 Not sure if this function is needed as things are
  // done somewhat differently compared to the 2D case.

  // T-connection configurations etc between adjacent volumes are currently
  // not handled in all cases. Check that the surface represents an entire boundary for
  // both adjacent volumes
  Point v1_p1 = bd_sf1->volumeParameter(sf1_start1, sf1_start2);
  Point v1_p2 = bd_sf1->volumeParameter(sf1_end1, sf1_start2);
  Point v1_p3 = bd_sf1->volumeParameter(sf1_start1, sf1_end2);
  Point v1_p4 = bd_sf1->volumeParameter(sf1_end1, sf1_end2);

  Point v2_p1 = bd_sf2->volumeParameter(sf2_start1, sf2_start2);
  Point v2_p2 = bd_sf2->volumeParameter(sf2_end1, sf2_start2);
  Point v2_p3 = bd_sf2->volumeParameter(sf2_start1, sf2_end2);
  Point v2_p4 = bd_sf2->volumeParameter(sf2_end1, sf2_end2);

  Array<double, 6> vol1_span = vol1->parameterSpan();
  Array<double, 6> vol2_span = vol2->parameterSpan();
  // RectDomain dom1 = srf1->parameterDomain();
  // RectDomain dom2 = srf2->parameterDomain();

  double ptol = 1.0e-10;
  int bd1, bd2;  // Specifies the volume boundaries corresponding to 
  // the current faces
  // 0 = umin, 1 = umax, 2 = vmin,  3 = vmax, 4 = wmin, 5 = wmax
  int orientation1, orientation2;
  bool swap1, swap2; // Denotes whether the u- and v-dir of the
		     // oppsite object is swapped (wrt to bd_sf1 & bd_sf2).
  bd1 = bd_sf1->whichBoundary(epsge, orientation1, swap1);
  bd2 = bd_sf2->whichBoundary(epsge, orientation2, swap2);

  if (bd1 < 0 || bd2 < 0)
      return;  // Unexpected situation

  int const_dir1 = bd1/2;
  int const_dir2 = bd2/2;

  // @@@ VSK, We should have a special treatment of the surface corners to
  // avoid creating new gaps towards other surfaces meeting in the corner
  // That is not implemented yet.
  shared_ptr<SplineVolume> v_1 = vol1;
  shared_ptr<SplineVolume> v_2 = vol2;

  SplineSurface* bd_sf1_spline = bd_sf1->asSplineSurface();
  SplineSurface* bd_sf2_spline = bd_sf2->asSplineSurface();
  assert(bd_sf1_spline != NULL && bd_sf2_spline != NULL);

  // We check if we are at a corner points along volume in all four sf points.
  bool atcorner1 = (fabs(bd_sf1_spline->startparam_u() - sf1_start1) < ptol &&
		    fabs(bd_sf1_spline->endparam_u() - sf1_end1) < ptol &&
		    fabs(bd_sf1_spline->startparam_v() - sf1_start2) < ptol &&
		    fabs(bd_sf1_spline->endparam_v() - sf1_end2) < ptol);
  bool atcorner2 = (fabs(bd_sf2_spline->startparam_u() - sf2_start1) < ptol &&
		    fabs(bd_sf2_spline->endparam_u() - sf2_end1) < ptol &&
		    fabs(bd_sf2_spline->startparam_v() - sf2_start2) < ptol &&
		    fabs(bd_sf2_spline->endparam_v() - sf2_end2) < ptol);
  bool keep_first = false, keep_second = false;
  if (!atcorner1 && !atcorner2)
      return;  // Specific T-situation. Currently not handled
  else if (!atcorner1)
  {   // The inside of vol1 matches the outer bound of vol2.
      // In const dir we pick the hole domain.
      // Pick the relevant part of volume one, modify volume two
      double u1 = (const_dir1 == 0) ? vol1_span[0] : sf1_start1;
      double v1 = (const_dir1 == 1) ? vol1_span[2] : ((const_dir1 == 0) ? sf1_start1 : sf1_start2);
      double w1 = (const_dir1 == 2) ? vol1_span[4] : sf1_start2;
      double u2 = (const_dir1 == 0) ? vol1_span[1] : sf1_end1;
      double v2 = (const_dir1 == 1) ? vol1_span[3] : ((const_dir1 == 0) ? sf1_end1 : sf1_end2);
      double w2 = (const_dir1 == 2) ? vol1_span[5] : sf1_end2;
      shared_ptr<SplineVolume> vol3(v_1->subVolume(u1, v1, w1, u2, v2, w2, ptol));
      v_1 = vol3;
      vol1_span = v_1->parameterSpan();
      keep_first = true;
  }
  else if (!atcorner2)
  {
      // Pick the relevant part of surface two, modify surface one
      double u1 = (const_dir2 == 0) ? vol2_span[0] : sf2_start1;
      double v1 = (const_dir2 == 1) ? vol2_span[2] : ((const_dir2 == 0) ? sf2_start1 : sf2_start2);
      double w1 = (const_dir2 == 2) ? vol2_span[4] : sf2_start2;
      double u2 = (const_dir2 == 0) ? vol2_span[1] : sf2_end1;
      double v2 = (const_dir2 == 1) ? vol2_span[3] : ((const_dir2 == 0) ? sf2_end1 : sf2_end2);
      double w2 = (const_dir2 == 2) ? vol2_span[5] : sf2_end2;
      shared_ptr<SplineVolume> vol3(v_2->subVolume(u1, v1, w1, u2, v2, w2, ptol));
      v_2 = vol3;
      vol2_span = v_2->parameterSpan();
      keep_first = true;
  }





  bool opposite = false;
  double t1 = (bd1 == 0 ||  bd1 == 1) ? f1_p1[1] : f1_p1[0];
  double t2 = (bd1 == 0 ||  bd1 == 1) ? f1_p2[1] : f1_p2[0];
  double t3 = (bd2 == 0 ||  bd2 == 1) ? f2_p1[1] : f2_p1[0];
  double t4 = (bd2 == 0 ||  bd2 == 1) ? f2_p2[1] : f2_p2[0];
  if ((t2 - t1)*(t4 -t3) < 0.0)
      opposite = true;
  if ((same_orient1 && !same_orient2) || (!same_orient1 && same_orient2))
      opposite = !opposite;
  if (same_orientation != NULL && !(*same_orientation))
      opposite = !opposite;


  int orientation = 1; // @@sbr This must be set!
  try {
      GeometryTools::averageBoundaryCoefs(s1, bd1, keep_first, s2, bd2, keep_second, true, 
					  vertex1, true, vertex2, orientation);
  }
  catch(...)
  {
      return;
  }

  // Update boundary curves
  bool updated1, updated2;
  updated1 = bd_cv1->updateIsoCurves();
  updated2 = bd_cv2->updateIsoCurves();

  double mdist1, mdist2;
  int nmb_sample = 200;
  checkBoundaryDist(bd_cv1, bd_cv2, start1, end1, start2, end2,
		    nmb_sample, mdist1, mdist2);
  if (getenv("DEBUG") && (*getenv("DEBUG")) == '1')
  {
      std::cout << "removeGapSpline, distances: " << mdist1 << ", ";
      std::cout << mdist2 << std::endl;
  }

  // @@sbr201111 I guess a crv for domain to be altered may be more
  // appropriate.  Otherwise we expect the domain to correspond to
  // iso-curves.
#endif


}

}
