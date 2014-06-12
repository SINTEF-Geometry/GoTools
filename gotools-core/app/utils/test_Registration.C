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



#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include "GoTools/utils/RegistrationUtils.h"
#include "GoTools/geometry/ObjectHeader.h"



using namespace std;
using namespace Go;


// Random double in range [0.0, 1.0]
double rnd_d()
{
  return (double)rand() / RAND_MAX;
}


// Random, uniformly distributed point in ball, or on ball sphere
Point random_sphere_point(double radius, double edge_only)
{
  double h = radius * (rnd_d() * 2.0 - 1.0);
  double angle = rnd_d() * 2.0 * M_PI;
  double sqr_rad = radius * radius - h * h;
  if (sqr_rad < 0.0)
    sqr_rad = 0.0;  // Just in case, should not happen
  double rad = sqrt(sqr_rad);
  Point result(h, rad * cos(angle), rad * sin(angle));

  if (!edge_only)
    {
      double factor = pow(rnd_d(), 1.0/3.0);
      result *= factor;
    }

  return result;
}


// Sum of square distances between pairs of points
double sq_dist(const vector<Point>& p1, const vector<Point>& p2)
{
  int n_pts = p1.size();
  double result = 0.0;
  for (int i = 0; i < n_pts; ++i)
    result += p1[i].dist2(p2[i]);
  return result;
}


void registration_test_sets(vector<Point> &fixed, vector<Point> &movable, bool allow_rescaling, int n_points,
			    double max_point_dist, double max_disturb_dist, double max_translate_dist,
			    double min_rescale, double max_rescale, vector<vector<double> >& rot_matrix, Point& translate, double& rescale)
{
  // Create first collection of points
  fixed.resize(0);
  for (int i = 0; i < n_points; ++i)
    fixed.push_back(random_sphere_point(max_point_dist, false));

  // Copy first collection to second, with some disturbance
  movable.resize(0);
  for (int i = 0; i < n_points; ++i)
    movable.push_back(fixed[i] + random_sphere_point(max_disturb_dist, false));

  // Create random rotation, rescaling and translation, and apply them on second point set
  Point rot_axis = random_sphere_point(1.0, true);
  double rot_angle = 2.0 * M_PI * rnd_d();
  double cos_angle = cos(rot_angle);
  double sin_angle = sin(rot_angle);

  rot_matrix.resize(3);
  for (int i = 0; i < 3; ++i)
    {
      rot_matrix[i].resize(3);
      for (int j = 0; j < 3; ++j)
	rot_matrix[i][j] = rot_axis[i] * rot_axis[j] * (1.0 - cos_angle);
    }

  rot_matrix[0][0] += cos_angle;
  rot_matrix[1][1] += cos_angle;
  rot_matrix[2][2] += cos_angle;

  rot_matrix[1][2] -= sin_angle * rot_axis[0];
  rot_matrix[2][1] += sin_angle * rot_axis[0];
  rot_matrix[2][0] -= sin_angle * rot_axis[1];
  rot_matrix[0][2] += sin_angle * rot_axis[1];
  rot_matrix[0][1] -= sin_angle * rot_axis[2];
  rot_matrix[1][0] += sin_angle * rot_axis[2];

  rescale = 1.0;
  if (allow_rescaling)
    rescale = rnd_d() * (max_rescale - min_rescale) + min_rescale;

  translate = random_sphere_point(max_translate_dist, false);

  for (int i = 0; i < n_points; ++i)
    {
      Point old_p = movable[i];
      Point new_p(0.0, 0.0, 0.0);
      for (int i = 0; i < 3; ++i)
	for (int j = 0; j < 3; ++j)
	  new_p[i] += old_p[j] * rot_matrix[i][j];
      new_p *= rescale;
      new_p += translate;
      movable[i] = new_p;
    }
}


void reportRegistrationError(RegistrationInput params, RegistrationResult result)
{
  cout << endl << "Dropping registration parameters and results:" << endl;
  cout << "Input parameters:" << endl;
  cout << "  Max Newton interations       = " << params.max_newton_iterations_ << endl;
  cout << "  Newton tolerance             = " << params.newton_tolerance_ << endl;
  cout << "  Calculate Newton tol weights ? " << (params.calculate_tolerance_weights_ ? "Yes" : "No") << endl;
  cout << "    Rotation weight            = " << params.tolerance_weight_rotation_ << endl;
  cout << "    Translation weight         = " << params.tolerance_weight_translation_ << endl;
  cout << "    Rescaling weight           = " << params.tolerance_weight_rescale_ << endl;
  cout << "  Max solve interations        = " << params.max_solve_iterations_ << endl;
  cout << "  Solve tolerance              = " << params.solve_tolerance_ << endl;
  cout << "Result:" << endl;
  cout << "  Result type            : ";
  switch (result.result_type_)
    {
    case RegistrationOK:
      cout << "Registration was OK" << endl;
      break;
    case TooFewPoints:
      cout << "Less than three input points" << endl;
      break;
    case PointSetSizeDiff:
      cout << "Different size of input point sets" << endl;
      break;
    case AreaTooSmall:
      cout << "Could not find point tripples far enough from being colinear" << endl;
      break;
    case SolveFailed:
      cout << "Solving linear system failed" << endl;
      break;
    default:
      cout << "Unknown return type (should not happen)" << endl;
    }
  cout << "  Last Newton iteration  = " << result.last_newton_iteration_ << endl;
  cout << "  Last weighted change   = " << result.last_change_ << endl;
  cout << "  Last solve return code = " << result.solve_result_ << endl;
  cout << endl;
}


// Run registration test
void test_registration(bool allow_rescaling, int n_points = 100,
		       double max_point_dist = 100.0, double max_disturb_dist = 1.0, double max_translate_dist = 35.0,
		       double min_rescale = 0.2, double max_rescale = 10.0)
{
  cout << endl << "Running registration test" << endl;
  cout << "Numer of points = " << n_points << endl;
  cout << "Rescaling option = " << (allow_rescaling ? "On" : "Off") << endl;

  // Create point sets
  vector<Point> fixed;
  vector<Point> movable;
  vector<vector<double> > test_rot_matrix;
  Point test_translate;
  double test_rescale;
  registration_test_sets(fixed, movable, allow_rescaling, n_points,
			 max_point_dist, max_disturb_dist, max_translate_dist,
			 min_rescale, max_rescale, test_rot_matrix, test_translate, test_rescale);
  cout << "Ratio max disturb/point distance = " << (max_disturb_dist / max_point_dist) << endl;
  if (allow_rescaling)
    cout << "Chosen rescaling = " << test_rescale << endl;

  double distance_before = sq_dist(fixed, movable);
  cout << "Square distance before registration = " << distance_before << endl;

  // Perform registration
  RegistrationInput params;
  params.max_newton_iterations_ = 20;
  RegistrationResult reg_result = registration(fixed, movable, allow_rescaling, params);
  bool run_OK = reg_result.result_type_ == RegistrationOK;

  // Check registration
  if (!run_OK)
    {
      cout << "**** Error: Registration failed!!!" << endl;
      reportRegistrationError(params, reg_result);
    }
  else
    {
      vector<vector<double> > rotation = reg_result.rotation_matrix_;
      Point translate = reg_result.translation_;
      double rescale = reg_result.rescaling_;
      vector<Point> after_registration;
      for (int i = 0; i < n_points; ++i)
	{
	  Point p_bef = movable[i];
	  Point p_aft(translate);
	  for (int i = 0; i < 3; ++i)
	    for (int j = 0; j < 3; ++j)
	      p_aft[i] += p_bef[j] * rotation[i][j] * rescale;
	  after_registration.push_back(p_aft);
	}

      double distance_after = sq_dist(fixed, after_registration);
      cout << "Square distance after registration = " << distance_after << endl;

      // Check if this is the best result in neighbourhood
      int bad_count = 0;
      for (int i = 0; i < (allow_rescaling ? 14 : 12); ++i)
	{
	  // i = 0-2, small rotation around unit vector in coordinate system
	  // i = 3-5, as 0-2, but opposite rotation direction
	  // i = 6-8, small translation along unit vector
	  // i = 9-11, as 6-8 but opposite direction
	  // i = 12, small up-scaling
	  // i = 13, small down-scaling
	  vector<Point> disturbed;
	  for (int i = 0; i < after_registration.size(); ++i)
	    disturbed.push_back(after_registration[i]);
	  if (i < 6)
	    {
	      int ax0 = i;
	      double angle = 0.0001;
	      if (i >= 3)
		{
		  angle = -angle;
		  ax0 -= 3;
		}
	      double cos_t = cos(angle);
	      double sin_t = sin(angle);
	      int ax1 = (ax0 + 1) % 3;
	      int ax2 = (ax0 + 2) % 3;
	      for (int j = 0; j < disturbed.size(); ++j)
		{
		  double p1 = disturbed[j][ax1];
		  double p2 = disturbed[j][ax2];
		  disturbed[j][ax1] = p1 * cos_t + p2 * sin_t;
		  disturbed[j][ax2] = p2 * cos_t - p1 * sin_t;
		}
	    }
	  else if (i < 12)
	    {
	      int ax0 = i - 6;
	      double delta = 0.001;
	      if (i >= 9)
		{
		  delta = -delta;
		  ax0 -= 3;
		}
	      Point delta_pt(0.0, 0.0, 0.0);
	      delta_pt[ax0] = delta;
	      for (int j = 0; j < disturbed.size(); ++j)
		disturbed[j] += delta_pt;
	    }
	  else
	    {
	      double delta = 0.001;
	      if (i == 13)
		delta = -delta;
	      delta += 1.0;
	      for (int j = 0; j < disturbed.size(); ++j)
		disturbed[j] *= delta;
	    }
	  double distance_disturb = sq_dist(fixed, disturbed);
	  double delta_disturb = distance_disturb - distance_after;
	  if (delta_disturb < 0.0)
	    {
	      cout << "**** Error when checking disturbance for i = " << i << "   Change should be non-negative, but is " << delta_disturb << endl;
	      ++bad_count;
	    }
	}
      if (bad_count == 0)
	cout << "Neighbourhood tests indicate this is the locally best solution" << endl;
      else
	cout << "**** Error: Neighbourhood tests show " << bad_count << " improvement direction(s)" << endl;

      // Test if we get same result without performing common rotation, rescaling and translation
      vector<Point> untransformed;
      for (int i = 0; i < n_points; ++i)
	{
	  Point p = (movable[i] - test_translate) / test_rescale;
	  Point q(0.0, 0.0, 0.0);
	  for (int j = 0; j < 3; ++j)
	    {
	      double sum = 0.0;
	      for (int k = 0; k < 3; ++k)
		sum += test_rot_matrix[k][j] * p[k];   // Apply the inverse (= the transpose) of the rotation matrix
	      q[j] = sum;
	    }
	  untransformed.push_back(q);
	}

      RegistrationInput params2;
      RegistrationResult reg_result2 = fineRegistration(fixed, untransformed, allow_rescaling, params);
      bool run_OK2 = reg_result2.result_type_ == RegistrationOK;

      if (!run_OK2)
	{
	  cout << "**** Error: Untransformed registration test failed!!!" << endl;
	  reportRegistrationError(params2, reg_result2);
	}
      else
	{
	  vector<vector<double> > rotation2 = reg_result2.rotation_matrix_;
	  Point translate2 = reg_result2.translation_;
	  double rescale2 = reg_result2.rescaling_;
	  double untransform_diff = 0.0;
	  for (int i = 0; i < n_points; ++i)
	    {
	      Point p_bef = untransformed[i];
	      Point p_aft(translate2);
	      for (int i = 0; i < 3; ++i)
		for (int j = 0; j < 3; ++j)
		  p_aft[i] += p_bef[j] * rotation2[i][j] * rescale2;
	      untransform_diff += p_aft.dist2(after_registration[i]);
	    }
	  cout << "Square distance to untransformed registration = " << untransform_diff << endl;
	}
    }
}


void writePoints(ostream& os, const vector<Point>& pts, int r, int g, int b)
{
  os << endl << "400 1 0 4 " << r << " " << g << " " << b << " 255" << endl;
  int n_pts = pts.size();
  os << n_pts << endl;
  for (int i = 0; i < n_pts; ++i)
    os << pts[i] << endl;
}


void gui_test_registration(string filename, bool allow_rescaling, int n_points = 100,
			   double max_point_dist = 100.0, double max_disturb_dist = 1.0, double max_translate_dist = 35.0,
			   double min_rescale = 0.2, double max_rescale = 10.0)
{
  // Create point sets
  vector<Point> fixed;
  vector<Point> movable;
  vector<vector<double> > test_rot_matrix;
  Point test_translate;
  double test_rescale;
  registration_test_sets(fixed, movable, allow_rescaling, n_points,
			 max_point_dist, max_disturb_dist, max_translate_dist,
			 min_rescale, max_rescale, test_rot_matrix, test_translate, test_rescale);

  // Find registration
  RegistrationInput params;
  RegistrationResult reg_result = registration(fixed, movable, allow_rescaling, params);
  if (reg_result.result_type_ == RegistrationOK)
    {

      // Apply registration
      vector<vector<double> > rotation = reg_result.rotation_matrix_;
      Point translate = reg_result.translation_;
      double rescale = reg_result.rescaling_;
      vector<Point> moved;
      for (int i = 0; i < n_points; ++i)
	{
	  Point p_bef = movable[i];
	  Point p_aft(translate);
	  for (int i = 0; i < 3; ++i)
	    for (int j = 0; j < 3; ++j)
	      p_aft[i] += p_bef[j] * rotation[i][j] * rescale;
	  moved.push_back(p_aft);
	}

      // Write fixed points set (red), the other point set in original position (green)
      // and the other point set after applying registration (blue)
      ofstream outfile(filename);
      writePoints(outfile, fixed, 255, 0, 0);
      writePoints(outfile, movable, 0, 255, 0);
      writePoints(outfile, moved, 0, 0, 255);
      outfile.close();
    }
}


int main()
{
  //  unsigned int seed = time(NULL);
  unsigned int seed = 1400228652;
  cout << "Random seed = " << seed << endl;
  srand(seed);

  test_registration(false);
  test_registration(true);
  test_registration(false, 100, 10.0, 10.0);
  test_registration(true, 100, 10.0, 10.0);
  test_registration(false, 10000);
  test_registration(true, 10000);

  gui_test_registration("reg_test_rescale_off.g2", false, 100, 100.0, 50.0);
  gui_test_registration("reg_test_rescale_on.g2", true, 100, 100.0, 50.0);
}
