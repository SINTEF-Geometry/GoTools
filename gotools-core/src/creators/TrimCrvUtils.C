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


#include "GoTools/creators/TrimCrvUtils.h"

#include "GoTools/geometry/GeometryTools.h"
#include "GoTools/geometry/PointCloud.h"
#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/geometry/Utils.h"
#include "GoTools/geometry/BoundedSurface.h"
//#include "GoTools/lrsplines2D/LRSplineSurface.h"
#include "GoTools/creators/ApproxCurve.h"

#include <fstream>
#include <math.h>

using namespace Go;
using std::vector;
using std::cout;
using std::endl;
using std::ifstream;

//#define DEBUG



//===========================================================================
vector<double>  
TrimCrvUtils::readTrimPoints(ifstream& filein, Point& translate_vec)
//===========================================================================
{
  vector<double> pts_2d;
  translate_vec = Point(0.0, 0.0, 0.0);
  const bool translate_model = true;//false;
  if (translate_model)
    {
      MESSAGE("Switched off translation of the model to the origin!");
    }

  ObjectHeader header;
  header.read(filein);
  if (header.classType() == Class_PointCloud)
    {
      PointCloud3D pt_cloud;
      pt_cloud.read(filein);
      if (translate_model)
	{
	  translateToOrigin(pt_cloud, translate_vec);
	}
      const int dim = pt_cloud.dimension();
      assert(dim == 3);
      const int num_pts = pt_cloud.numPoints();
      vector<double> pts_3d(pt_cloud.rawData(), pt_cloud.rawData() + dim*num_pts);
      pts_2d.resize(num_pts*2);
      for (int ki = 0; ki < num_pts; ++ki)
	{
	  pts_2d[ki*2] = pts_3d[ki*dim];
	  pts_2d[ki*2+1] = pts_3d[ki*dim+1];
	}
    }
  else if (header.classType() == Class_LineCloud)
    {
      LineCloud line_cloud;
      line_cloud.read(filein);
      if (translate_model)
	{
	  translateToOrigin(line_cloud, translate_vec);
	}
      // Assuming the dimension is 3. Which is a requirement by LineCloud, strangely enough.
      const int dim = 3;
      const int num_lines = line_cloud.numLines();
      const int num_pts = num_lines + 1;
      pts_2d.resize(num_pts*2);
      const double* raw_data = line_cloud.rawData();
      pts_2d[0] = raw_data[0];
      pts_2d[1] = raw_data[1];
      for (int ki = 1; ki < num_pts; ++ki)
	{
	  pts_2d[ki*2] = raw_data[ki*dim*2-dim];
	  pts_2d[ki*2+1] = raw_data[ki*dim*2-dim+1];
	}
    }
  else
    {
      MESSAGE("Object  type " << header.classType() << " is not supported.");
    }

  return pts_2d;
}


//==========================================================================
vector<vector<double> > 
TrimCrvUtils::splitTrimPoints(vector<double>& pts_2d, 
			      double epsgeo, double kink_tol)
//==========================================================================
{
  vector<vector<double> > split_pts_2d;
  const bool do_not_split = false;//true;

  if (do_not_split)
    {
      MESSAGE("Turned off splitting the boundary in kinks.");
      split_pts_2d.push_back(pts_2d);
      return split_pts_2d;
    }
    
  try
    {
      split_pts_2d = TrimCrvUtils::splitCurvePointsInKinks(pts_2d, kink_tol);
    }
  catch (...)
    {
      MESSAGE("Failed ...");
      return split_pts_2d;
    }

  // Since we are reading trim pts we expect them to form a loop.
  makeConnectedLoop(split_pts_2d, epsgeo);
  return split_pts_2d;
}
    
//==========================================================================
void TrimCrvUtils::makeConnectedLoop(vector<vector<double> >& trim_seqs_2d, 
				     double epsgeo)
//==========================================================================
{
  // Ensure that the trim sequences join up in a loop

  int last_ind = trim_seqs_2d.size() - 1;
  // We make sure that first and last points are within epsgeo.
  Point start_pt(trim_seqs_2d[0][0], trim_seqs_2d[0][1]);
  Point end_pt(trim_seqs_2d[last_ind][trim_seqs_2d[last_ind].size()-2],
	       trim_seqs_2d[last_ind][trim_seqs_2d[last_ind].size()-1]);
  double dist_end_pts = start_pt.dist(end_pt);
  //cout << "Distance end points: " << dist_end_pts << endl;
  if (dist_end_pts > epsgeo)
    {
      //cout << "Adding point to make sure we create a closed loop." << endl;
      if (trim_seqs_2d[last_ind].size() > 4)
	trim_seqs_2d[last_ind].insert(trim_seqs_2d[last_ind].end(), 
				      start_pt.begin(), start_pt.end());
      else
	trim_seqs_2d[0].insert(trim_seqs_2d[0].begin(), end_pt.begin(),
			       end_pt.end());
    }
}

//===========================================================================
vector<vector<double> > 
TrimCrvUtils::extractConstParSeqs(vector<double>& pts_2d, 
				  int parix, double par,
				  int nmb_match,
				  double tol1, double tol2)
//===========================================================================
{
  vector<vector<double> > trim_seqs;
  int ki, kj, kr, kh;
  int prev = 0;
  for (kj=0; kj<pts_2d.size(); )
    {
      // Traverse until first occurance of point on specified constant
      // parameter using the expected strictest tolerance 
      for (; kj<pts_2d.size(); kj+=2)
	if (fabs(pts_2d[kj+parix]-par) <= tol1)
	  break;
	
      // Traverse all points on constant parameter curve
      for (kr=kj+2; kr<pts_2d.size(); kr+=2)
	if (fabs(pts_2d[kr+parix]-par) > tol2)
	  break;

      // Traverse back to find a good match end point
      for (; kr > kj+2; kr-=2)
	if (fabs(pts_2d[kr-2+parix]-par) <= tol1)
	  break;
	
      if ((kr-kj)/2 >= nmb_match)
	{
	  // A sequence is found
	  if (kj > prev)
	    {
	      vector<double> curr_seq(pts_2d.begin()+prev, pts_2d.begin()+kj+2);
	      trim_seqs.push_back(curr_seq);
	    }

	  vector<double> const_seq(pts_2d.begin()+kj, pts_2d.begin()+kj+2);
	  const_seq.insert(const_seq.end(), pts_2d.begin()+kr-2, 
			   pts_2d.begin()+kr);
	  trim_seqs.push_back(const_seq);
	    
	  prev = kr-2;
	  kj = kr;
	}
      else
	kj = kr;
    }
  if (prev < pts_2d.size()-2)
    {
      vector<double> curr_seq(pts_2d.begin()+prev, pts_2d.end());
      trim_seqs.push_back(curr_seq);
    }
	

  return trim_seqs;
}

//===========================================================================
SplineCurve* 
TrimCrvUtils::createVariableOffsetDistFunction(const SplineCurve& par_cv, 
					       const vector<double>& pts_2d)
//===========================================================================
{
  MESSAGE("Not implemented yet!");

  // @@sbr201410 Remember to only use points which are inside the surface domain! Which means
  // that we should either include that variable in the argument list or handle this on the
  // outside. Or perhaps remove those points?


  return 0;
}


//===========================================================================
SplineCurve* 
TrimCrvUtils::createOffsetTrimCurve(const SplineCurve& par_cv, 
				    const SplineCurve& offset_dist)
//===========================================================================
{
  MESSAGE("Not implemented yet!");

  // We offset 10 % more than the exact value.
  const double dist_mult = 1.1;

  return 0;
}


//===========================================================================
shared_ptr<SplineCurve> 
TrimCrvUtils::approximateTrimPts(vector<double> trim_pts, int dim, 
				 double epsgeo, int max_iter)
//===========================================================================
{
  int num_pts = trim_pts.size()/dim;

  // Check for the trivial case
  if (num_pts == 2)
    {
      Point pnt1(trim_pts[0], trim_pts[1]);
      Point pnt2(trim_pts[dim], trim_pts[dim+1]);
      double len = std::max(0.01, pnt1.dist(pnt2));
      shared_ptr<SplineCurve> curve(new SplineCurve(pnt1, 0.0,
						    pnt2, len));
      return curve;
    }
							  
  // Otherwise, approximate point sequence
  vector<double> pts_xy(2*num_pts);
  vector<double> params_xy(num_pts);
  // vector<double> pts_z(num_pts);
  // vector<double> params(num_pts);
  // pts_z[0] = trim_pts[dim-1];
  pts_xy[0] = trim_pts[0];
  pts_xy[1] = trim_pts[1];
  // params[0] = 0.0;
  params_xy[0] = 0.0;
  for (int ki = 1; ki < num_pts; ++ki)
    {
      Point prev_pt(trim_pts.begin() + (ki-1)*dim, trim_pts.begin() + ki*dim);
      Point curr_pt(trim_pts.begin() + ki*dim, trim_pts.begin() + (ki+1)*dim);
      double dist = prev_pt.dist(curr_pt);
      //	    params[ki] = params[ki-1] + dist;
      //	    pts_z[ki] = trim_pts[ki*dim+2];
      pts_xy[2*ki] = trim_pts[ki*dim];
      pts_xy[2*ki+1] = trim_pts[ki*dim+1];
      Point prev_par_pt(pts_xy.begin() + (ki-1)*2, pts_xy.begin() + ki*2);
      Point curr_par_pt(pts_xy.begin() + ki*2, pts_xy.begin() + (ki+1)*2);
      double dist_par = prev_par_pt.dist(curr_par_pt);
      params_xy[ki] = params_xy[ki-1] + dist_par;
    }


  ApproxCurve approx_cv_2d(pts_xy, params_xy, 2, epsgeo);
  double maxdist_2d = -1.0, avgdist_2d = -1.0;
  //cout << "Approximating the 2D curve." << endl;
  shared_ptr<SplineCurve> spline_cv_appr_2d = approx_cv_2d.getApproxCurve(maxdist_2d, avgdist_2d, max_iter);
  //cout << "Done approximating the 2D curve." << endl;

  return spline_cv_appr_2d;
}


//===========================================================================
//===========================================================================
void TrimCrvUtils::moveCurveCoefsInsideSurfDomain(const ParamSurface* sf,
						  vector<shared_ptr<SplineCurve> >& par_cvs)
//===========================================================================
{
  RectDomain rect_dom = sf->containingDomain();
  for (size_t ki = 0; ki < par_cvs.size(); ++ki)
    {
      const int dim = par_cvs[ki]->dimension();
      assert(dim == 2);
	    
      for (auto iter = par_cvs[ki]->coefs_begin(); iter != par_cvs[ki]->coefs_end(); iter += 2)
	{
	  bool inside_dom = sf->inDomain(iter[0], iter[1]);
	  if (!inside_dom)
	    {
	      //MESSAGE("The point (" << iter[0] << ", " << iter[1] << ") was not inside the domain!");
	      if (iter[0] < rect_dom.umin())
		{
		  iter[0] = rect_dom.umin();
		  //MESSAGE("Moved inside!");
		}
	      else if (iter[0] > rect_dom.umax())
		{
		  iter[0] = rect_dom.umax();
		  //MESSAGE("Moved inside!");
		}
	      if (iter[1] < rect_dom.vmin())
		{
		  iter[1] = rect_dom.vmin();
		  //MESSAGE("Moved inside!");
		}
	      else if (iter[1] > rect_dom.vmax())
		{
		  iter[1] = rect_dom.vmax();
		  //MESSAGE("Moved inside!");
		}
	    }
	}
    }
}


//===========================================================================
vector<shared_ptr<Line> > 
TrimCrvUtils::approximateCurve(const SplineCurve& par_curve, double epspar)
//===========================================================================
{
  // Currently splitting when we encounter a coef more than epspar away from the line segment.
  // A better approach resulting in fewer segments is to always look for the largest distance.
  // But this takes more time to perform.
  vector<shared_ptr<SplineCurve> > spline_segments;
  spline_segments.push_back(shared_ptr<SplineCurve>(par_curve.clone()));
  int curr_ind = 0;
  // Wost case scenario is that we must keep on splitting until there are no more inner knots left ...
  // At least the approach should be a finite algorithm (given that epspar is not 0.0).
  while (curr_ind < spline_segments.size())
    {
      SplineCurve* curr_segment = spline_segments[curr_ind].get();
      int max_dist_ind;
      double max_coef_dist;
      double max_coef_par;
      // @@sbr201410 Perhaps replace max_coef_par with max_coef_par? I.e. the closest knot par?
      // In that way we are guaranteed to succeed.
      largestLineCoefDistance(*curr_segment,
			      max_dist_ind, max_coef_dist, max_coef_par);
      if (max_coef_dist > epspar)
	{
	  // We must split the curve.
	  shared_ptr<SplineCurve> left(curr_segment->subCurve(curr_segment->startparam(), max_coef_par));
	  shared_ptr<SplineCurve> right(curr_segment->subCurve(max_coef_par, curr_segment->endparam()));
	  spline_segments[curr_ind] = left;
	  spline_segments.insert(spline_segments.begin() + curr_ind + 1, right);
	}
      else
	{
	  ++curr_ind;
	}
    }

  // We run through all the spline segments, creating line_segments from end points.
  vector<shared_ptr<Line> > line_segments;
  for (size_t ki = 0; ki < spline_segments.size(); ++ki)
    {
      Point from = spline_segments[ki]->ParamCurve::point(spline_segments[ki]->startparam());
      Point to = spline_segments[ki]->ParamCurve::point(spline_segments[ki]->endparam());
      line_segments.push_back(shared_ptr<Line>(new Line(from, to,
							spline_segments[ki]->startparam(),
							spline_segments[ki]->endparam())));
    }

  return line_segments;
}


//===========================================================================
void TrimCrvUtils::largestLineCoefDistance(const SplineCurve& curr_segment,
					   int& max_dist_ind, 
					   double& max_coef_dist, 
					   double& max_coef_par)
//===========================================================================
{
  MESSAGE("Under construction!");
  const int dim = curr_segment.dimension();
  assert(dim == 2);
	
  Line line(curr_segment.ParamCurve::point(curr_segment.startparam()),
	    curr_segment.ParamCurve::point(curr_segment.endparam()),
	    curr_segment.startparam(),
	    curr_segment.endparam());
  max_coef_dist = -1.0;
  max_dist_ind = -1;
  int num_coefs = curr_segment.numCoefs();
  for (int ki = 0; ki < num_coefs; ++ki)
    {
      Point curr_coef(curr_segment.coefs_begin()[ki*dim],
		      curr_segment.coefs_begin()[ki*dim+1]);
      double clo_t, clo_dist;
      Point clo_pt;
      line.closestPoint(curr_coef, line.startparam(), line.endparam(),
			clo_t, clo_pt, clo_dist);
      if (clo_dist > max_coef_dist)
	{
	  max_coef_dist = clo_dist;
	  max_dist_ind = ki;
	}
    }

  max_coef_par = curr_segment.basis().grevilleParameter(max_dist_ind);
}


//===========================================================================
shared_ptr<SplineCurve> 
TrimCrvUtils::clipToDomain(const SplineCurve& par_cv, const Domain& domain)
//===========================================================================
{
  MESSAGE("Not implemented yet!");
  shared_ptr<SplineCurve> dummy;

  return dummy;
}


//===========================================================================
std::vector<std::vector<double> > 
TrimCrvUtils::splitCurvePointsInKinks(const std::vector<double>& trim_pts_2d,
				      double kink_tol)
//===========================================================================
{
  vector<vector<double> > segment_pts;

  // We must handle noise in the input data and can not use the direction between two consecutive points to define
  // the kinks in the boundary curve. We use the average of a certain number of points in both directions to
  // estimate the direction.
  const int num_avg_pts_both_dirs = 5; //10;//50;//10;//5;

  // @@sbr201410 Do not compare with previous direction, longer back!
  //MESSAGE("Check direction change vs pts[ki-num_avg_pts_both_dirs]");

  // const double kink_tol = 5e-01; // 0.1 deg => 5.7 degrees.
  // If the direction changes more than a certain tol after num_avg_pts, we define this as a kink and must
  // locate the point for which this change of direction is most prominent.
  double tol = 1.0e-4;
  const int num_pts = (int)trim_pts_2d.size()/2;
  double dd = Utils::distance_squared(trim_pts_2d.begin(),
				      trim_pts_2d.begin()+2,
				      trim_pts_2d.end()-2);
  bool closed = (dd < tol);
  if (num_pts <= 2*num_avg_pts_both_dirs)
    {
      segment_pts.push_back(trim_pts_2d);
      return segment_pts;
    }

  vector<int> segment_ind; // We store all kink indices.
  segment_ind.push_back(0);
  Point prev_dir = averageDirection(trim_pts_2d, 0, num_avg_pts_both_dirs);
  //MESSAGE("num_pts: " << num_pts);
  vector<Point> avg_dirs(num_pts, Point(0.0, 0.0));
  vector<double> angles(num_pts, 0.0);
  for (int ki = num_avg_pts_both_dirs; ki < num_pts - num_avg_pts_both_dirs; ++ki)
    {
      int from = std::max(0, ki - num_avg_pts_both_dirs);
      int to = std::min(ki + num_avg_pts_both_dirs, num_pts - 1);

      Point avg_dir = averageDirection(trim_pts_2d, from, to);
      avg_dir.normalize();
      avg_dirs[ki] = avg_dir;
      if (avg_dir.length() == 0)
	{ // Typically this means the we are turning around 180 degrees in the current point.
	  MESSAGE("Something strange with the input I suspect!");
	  continue;
	}
      //	    double angle = prev_dir.angle(avg_dir);
      if (avg_dirs[ki-num_avg_pts_both_dirs].length() > 0.0 && avg_dir.length() > 0.0)
	{
	  double angle = avg_dirs[ki-num_avg_pts_both_dirs].angle(avg_dir);
	  angles[ki] = angle;
	  if (angle > kink_tol)
	    {
	      //MESSAGE("ki = " << ki << ", angle = " << angle);
	      segment_ind.push_back(ceil(ki-0.5*num_avg_pts_both_dirs));
	    }
	}
      //	    prev_dir = avg_dir;
    }

  segment_ind.push_back(num_pts - 1);

#ifdef DEBUG
  { // We write to file the point and the segments surrounding it.
    vector<double> corner_pts;
    for (size_t ki = 0; ki < segment_ind.size(); ++ki)
      {
	corner_pts.push_back(trim_pts_2d[segment_ind[ki]*2]);
	corner_pts.push_back(trim_pts_2d[segment_ind[ki]*2+1]);
	corner_pts.push_back(0.0); // We place the pts in the z-plane.
      }

    LineCloud trim_pts;

    PointCloud3D end_pts(corner_pts.begin(), corner_pts.size()/3);
    std::ofstream fileout_debug("tmp/ptset_debug.g2");
    //MESSAGE("Writing to file the end pts of segments.");
    // end_pts.writeStandardHeader(fileout_debug);
    fileout_debug << "400 1 0 4 255 0 0 255" << endl;
    end_pts.write(fileout_debug);
    //MESSAGE("Done writing to file the end pts of segments.");
  }
#endif

  bool split_into_segments = true;
  if (split_into_segments)
    {
      // We run through and split based on segment_ind.
      // When there is a consecutive row of indices we choose the middle one.
      vector<int> thinned_segment_ind;
      int counter = 0;
      for (size_t ki = 0; ki < segment_ind.size(); ++ki)
	{
	  if ((ki < segment_ind.size() - 1) && (segment_ind[ki] + 1 == segment_ind[ki+1]))
	    {
	      ++counter;
	    }
	  else
	    {
	      if (counter > num_pts/4 && closed)
		{
		  // Split in both end points to avoid a too
		  // restricted loop
		  if (segment_ind[ki] > counter)
		    thinned_segment_ind.push_back(segment_ind[ki] - counter);
		  if (segment_ind[ki] < num_pts - 1)
		    thinned_segment_ind.push_back(segment_ind[ki]);
		  counter = 0;
		}
	      else if (counter > 1)
		{
		  vector<double> ang(counter);
		  for (size_t ka=0; ka<counter; ++ka)
		    {
		      int kb = segment_ind[ki] - counter + ka;
		      int from = std::max(0, kb - num_avg_pts_both_dirs);
		      int to = std::min(kb + num_avg_pts_both_dirs, num_pts - 1);

		      // Point avg_dir1 = averageDirection(trim_pts_2d, 
		      // 				    from, kb);
		      // Point avg_dir2 = averageDirection(trim_pts_2d, 
		      // 				    kb, to);
		      Point avg_dir1(trim_pts_2d[2*kb]-trim_pts_2d[2*from], 
				     trim_pts_2d[2*kb+1]-trim_pts_2d[2*from+1]);
		      Point avg_dir2(trim_pts_2d[2*to]-trim_pts_2d[2*kb], 
				     trim_pts_2d[2*to+1]-trim_pts_2d[2*kb+1]);
		      ang[ka] = avg_dir1.angle(avg_dir2);
		    }
		  int stop_break = 1;

		  // We are at the end of a consecutive stream of indices. 
		  // Picking the most significant one.
		  //int avg_ind = segment_ind[ki] - floor(0.5*counter);
		      
		  //cout << "avg_ind: " << avg_ind << endl;
		  double max_val = ang[0];
		  int max_ind = 0;
		  for (int kb=1; kb<(int)ang.size(); ++kb)
		    {
		      if (ang[kb] > max_val)
			{
			  max_val = ang[kb];
			  max_ind = kb;
			}
		    }
		  thinned_segment_ind.push_back(segment_ind[ki] - counter + max_ind);
		  counter = 0;
		}
	      else
		thinned_segment_ind.push_back(segment_ind[ki]);
	    }
	}

      // First and last point should be included, but we make sure anyway.
      if (thinned_segment_ind.front() != 0)
	{
	  thinned_segment_ind.insert(thinned_segment_ind.begin(), 0);
	}
      if (thinned_segment_ind.back() != num_pts - 1)
	{
	  thinned_segment_ind.push_back(num_pts - 1);
	}

#ifdef DEBUG
      { // We write to file the point and the segments surrounding it.
	vector<double> corner_pts;
	for (size_t ki = 0; ki < thinned_segment_ind.size(); ++ki)
	  {
	    corner_pts.push_back(trim_pts_2d[thinned_segment_ind[ki]*2]);
	    corner_pts.push_back(trim_pts_2d[thinned_segment_ind[ki]*2+1]);
	    corner_pts.push_back(0.0); // We place the pts in the z-plane.
	  }

	LineCloud trim_pts;

	PointCloud3D end_pts(corner_pts.begin(), corner_pts.size()/3);
	std::ofstream fileout_debug("tmp/ptset_debug2.g2");
	MESSAGE("Writing to file the end pts of segments.");
	// end_pts.writeStandardHeader(fileout_debug);
	fileout_debug << "400 1 0 4 0 255 0 255" << endl;
	end_pts.write(fileout_debug);
	MESSAGE("Done writing to file the end pts of segments.");
      }
#endif

      int num_segments = thinned_segment_ind.size() - 1;
      segment_pts.resize(num_segments);
      for (size_t ki = 0; ki < thinned_segment_ind.size() - 1; ++ki)
	{
	  int first = thinned_segment_ind[ki];
	  int last = thinned_segment_ind[ki+1];
	  segment_pts[ki].insert(segment_pts[ki].end(),
				 trim_pts_2d.begin() + 2*first, trim_pts_2d.begin() + 2*(last + 1));
	}

      //	MESSAGE("Returning the input for now ...");
      //	segment_pts.push_back(trim_pts_2d);
    }
  else
    {
      MESSAGE("Returning the input for now ...");
      segment_pts.push_back(trim_pts_2d);
    }

  return segment_pts;
}


//===========================================================================
void TrimCrvUtils::translateToOrigin(GeomObject& go_object, 
				     Point& translate_vec)
//===========================================================================
{
  const int dim = go_object.dimension();
  BoundingBox bd_box = go_object.boundingBox();
  Point box_center = 0.5*(bd_box.low() + bd_box.high());
  translate_vec = -box_center;

  translateObject(go_object, translate_vec);
}


//===========================================================================
void TrimCrvUtils::translateObject(GeomObject& go_object, 
				   const Point& translate_vec)
//===========================================================================
{
  // We can not handle all types of objects.
  if (go_object.instanceType() == Class_PointCloud)
    {
      PointCloud3D& pt_cloud = dynamic_cast<PointCloud3D&>(go_object);
      vector<double> transl_vec(translate_vec.begin(), translate_vec.end());
      Vector3D transl2_vec(transl_vec);
      pt_cloud.translate(transl2_vec);
    }
  else if (go_object.instanceType() == Class_LineCloud)
    {
      LineCloud& line_cloud = dynamic_cast<LineCloud&>(go_object);
      GeometryTools::translateLineCloud(translate_vec, line_cloud);
    }
  else if (go_object.instanceType() == Class_SplineCurve)
    {
      SplineCurve& spline_cv = dynamic_cast<SplineCurve&>(go_object);
      spline_cv.translateCurve(translate_vec);	
    }
  else if (go_object.instanceType() == Class_SplineSurface)
    {
      SplineSurface& spline_sf = dynamic_cast<SplineSurface&>(go_object);
      GeometryTools::translateSplineSurf(translate_vec, spline_sf);
    }
  // else if (go_object.instanceType() == Class_LRSplineSurface)
  //   {
  //     LRSplineSurface& lr_spline_sf = dynamic_cast<LRSplineSurface&>(go_object);
  //     lr_spline_sf.translate(translate_vec);
  //   }
  else
    {
      MESSAGE("Object type " << go_object.instanceType() << " is not supported yet.");
    }
}


//===========================================================================
void TrimCrvUtils::translateSurfaceDomain(ParamSurface* sf, 
					  const Point& translate_vec)
//===========================================================================
{
  RectDomain dom = sf->containingDomain();
  BoundedSurface* bdsf = dynamic_cast<BoundedSurface*>(sf);
  if (bdsf)
    dom = bdsf->underlyingSurface()->containingDomain();

  sf->setParameterDomain(dom.umin()+translate_vec[0],
			 dom.umax()+translate_vec[0],
			 dom.vmin()+translate_vec[1],
			 dom.vmax()+translate_vec[1]);
  // if (sf->instanceType() == Class_SplineSurface)
  //   {
  //     SplineSurface* spline_sf = dynamic_cast<SplineSurface*>(sf);
  //     spline_sf->setParameterDomain(spline_sf->startparam_u() + translate_vec[0],
  // 				    spline_sf->endparam_u() + translate_vec[0],
  // 				    spline_sf->startparam_v() + translate_vec[1],
  // 				    spline_sf->endparam_v() + translate_vec[1]);
  //   }
  // else if (sf->instanceType() == Class_LRSplineSurface)
  //   {
  //     LRSplineSurface* lr_spline_sf = dynamic_cast<LRSplineSurface*>(sf);
  //     lr_spline_sf->setParameterDomain(lr_spline_sf->startparam_u() + translate_vec[0],
  // 				       lr_spline_sf->endparam_u() + translate_vec[0],
  // 				       lr_spline_sf->startparam_v() + translate_vec[1],
  // 				       lr_spline_sf->endparam_v() + translate_vec[1]);
  //   }
  // else
  //   {
  //     MESSAGE("Object type " << sf->instanceType() << " is not supported yet.");
  //   }
}

#if 0
//===========================================================================
void TrimCrvUtils::scaleZ(ParamSurface& sf, double scale_factor)
//===========================================================================
{
  assert(sf.dimension() == 3);
  if (sf.instanceType() == Class_LRSplineSurface)
    {
      LRSplineSurface& lr_spline_sf = dynamic_cast<LRSplineSurface&>(sf);
      auto iter_begin = lr_spline_sf.basisFunctionsBegin();
      auto iter_end = lr_spline_sf. basisFunctionsEnd();
      auto iter = iter_begin;
      while (iter != iter_end)
	{
	  double gamma = iter->second->gamma();
	  Point coef = iter->second->Coef();
	  coef[2] *= scale_factor;
	  iter->second->setCoefAndGamma(coef, gamma);
	  ++iter;
	}
    }
  else if (sf.instanceType() == Class_LRSplineSurface)
    {
      SplineSurface& spline_sf = dynamic_cast<SplineSurface&>(sf);

    }
  else
    {
      MESSAGE("Z scaling of surface type " << sf.instanceType() << " not supported yet.");
    }
}
#endif

//===========================================================================
Point TrimCrvUtils::averageDirection(const vector<double>& pts_2d, 
				     int from, int to)
//===========================================================================
{
  const int num_pts = pts_2d.size()/2;
  assert((from <= to) && (from >= 0) && (to <= num_pts - 1));

  // For each point we define the direction in the point as the
  // average over the linear segments at both sides.
  Point avg_dir(0.0, 0.0);
  //for (int ki = from; ki < to; ++ki)
  to = std::min(to, num_pts);
  for (int ki = from+1; ki < to-1; ++ki)
    {
      Point left_leg(0.0, 0.0), right_leg(0.0, 0.0);
      if (ki > 0)
	{
	  left_leg[0] = pts_2d[ki*2] - pts_2d[(ki-1)*2];
	  left_leg[1] = pts_2d[ki*2+1] - pts_2d[(ki-1)*2+1];
	}
      if (ki < num_pts - 1)
	{
	  right_leg[0] = pts_2d[(ki+1)*2] - pts_2d[ki*2];
	  right_leg[1] = pts_2d[(ki+1)*2+1] - pts_2d[ki*2+1];
	}
      // if (left_leg.length() > 0.0)
      // 	left_leg.normalize();
      // if (right_leg.length() > 0.0)
      // 	right_leg.normalize();
      Point curr_dir = 0.5*(left_leg + right_leg);
      if (curr_dir.length() > 0.0)
	curr_dir.normalize();
      avg_dir += curr_dir;
    }

  if (avg_dir.length() == 0.0)
    {
      //MESSAGE("Oops, zero length ... from = " << from << ", to = " << to);
    }
  else
    {
      avg_dir.normalize();
    }

  return avg_dir;
}

//===========================================================================
double getEpsgeo(const vector<double>& pts_2d)
//===========================================================================
{
  int num_pts = pts_2d.size()/2;
  double sum_dist = 0.0;
  for (int ki = 0; ki < num_pts - 1; ++ki)
    {
      double dist = sqrt((pts_2d[ki*2] - pts_2d[(ki+1)*2])*(pts_2d[ki*2] - pts_2d[(ki+1)*2]) +
			 (pts_2d[ki*2+1] - pts_2d[(ki+1)*2+1])*(pts_2d[ki*2+1] - pts_2d[(ki+1)*2+1]));
      sum_dist += dist;
    }
  double avg_dist = sum_dist/(num_pts - 1);

  return avg_dist;
}
