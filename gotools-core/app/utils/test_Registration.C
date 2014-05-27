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



#include "GoTools/utils/RegistrationUtils.h"
#include <stdlib.h>
#include <time.h>
#include <math.h>



using namespace std;
using namespace Go;


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

double sq_dist(const vector<Point>& p1, const vector<Point>& p2)
{
  int n_pts = p1.size();
  double result = 0.0;
  for (int i = 0; i < n_pts; ++i)
    result += p1[i].dist2(p2[i]);
  return result;
}


void simple_reg_old1(vector<Point> &fixed, vector<Point> &movable)
{
  fixed.resize(0);
  fixed.push_back(Point(1.0,0.0,0.0));
  fixed.push_back(Point(0.0,1.0,0.0));
  fixed.push_back(Point(0.0,0.0,1.0));

  movable.resize(0);
  movable.push_back(Point(1.0,0.0,0.0));
  movable.push_back(Point(0.0,1.0,0.0));
  movable.push_back(Point(0.0,0.0,1.0));
}

void simple_reg_old2(vector<Point> &fixed, vector<Point> &movable)
{
  fixed.resize(0);
  fixed.push_back(Point(1.0,0.0,0.0));
  fixed.push_back(Point(0.0,1.0,0.0));
  fixed.push_back(Point(0.0,0.0,1.0));

  movable.resize(0);
  movable.push_back(Point(0.0,1.0,0.0));
  movable.push_back(Point(0.0,0.0,1.0));
  movable.push_back(Point(1.0,0.0,0.0));
}

void simple_reg(vector<Point> &fixed, vector<Point> &movable)
{
  fixed.resize(0);
  fixed.push_back(Point(1.0,2.0,2.0) / 3.0);
  fixed.push_back(Point(2.0,1.0,-2.0) / 3.0);
  fixed.push_back(Point(-2.0,2.0,-1.0) / 3.0);

  movable.resize(0);
  movable.push_back(Point(2.0,0.0,0.0));
  movable.push_back(Point(0.0,2.0,0.0));
  movable.push_back(Point(0.0,0.0,2.0));
}

void hard_reg(vector<Point> &fixed, vector<Point> &movable)
{
  // Collect some points
  int n_pts = 100;
  fixed.resize(0);
  for (int i = 0; i < n_pts; ++i)
    fixed.push_back(random_sphere_point(100.0, false));

  // Copy them, but with some disturbance
  movable.resize(0);
  for (int i = 0; i < n_pts; ++i)
    movable.push_back(fixed[i] + random_sphere_point(1.0, false));
    //movable.push_back(fixed[i]);
  // Rotate, rescale and translate copy
  Point rot_axis = random_sphere_point(1.0, true);
  double rot_angle = 2.0 * M_PI * rnd_d();
  double cos_angle = cos(rot_angle);
  double sin_angle = sin(rot_angle);

  vector<vector<double> > rot_matrix(3);
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

  double rescale = rnd_d() * 9.8 + 0.2;
  Point translate = random_sphere_point(35.0, true);

  for (int i = 0; i < n_pts; ++i)
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


int main()
{
  //  unsigned int seed = time(NULL);
  unsigned int seed = 1400228652;
  cout << "Random seed = " << seed << endl;
  srand(seed);

  // Set up point sets
  vector<Point> fixed_pts, copy_pts;
  hard_reg(fixed_pts, copy_pts);
  //simple_reg(fixed_pts, copy_pts);
  int n_pts = fixed_pts.size();

  double error_before = sq_dist(fixed_pts, copy_pts);
  cout << "Error before raw registration = " << error_before << endl;

  vector<vector<double> > raw_rotation;
  Point raw_translate;

  bool run_OK = registration(fixed_pts, copy_pts, raw_rotation, raw_translate);

  if (!run_OK)
    cout << "Raw registration failed!!!" << endl;
  else
    {
      vector<Point> after_registration;
      for (int i = 0; i < n_pts; ++i)
	{
	  Point p_bef = copy_pts[i];
	  Point p_aft(raw_translate);
	  for (int i = 0; i < 3; ++i)
	    for (int j = 0; j < 3; ++j)
	      p_aft[i] += p_bef[j] * raw_rotation[i][j];
	  after_registration.push_back(p_aft);
	}

      double error_after = sq_dist(fixed_pts, after_registration);
      cout << "Error after raw registration = " << error_after << endl;

      //      for (int i = 0; i < n_pts; ++i)
      //	cout << "Point " << i << " is " << after_registration[i] << endl;

    }

}
