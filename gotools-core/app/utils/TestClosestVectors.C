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



#include <vector>
#include <fstream>
#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/geometry/GoTools.h"
#include "GoTools/geometry/GeomObject.h"
#include "GoTools/geometry/Factory.h"
#include "GoTools/geometry/Utils.h"
#include "GoTools/utils/RegistrationUtils.h"
#include "GoTools/utils/ClosestPointUtils.h"



using namespace Go;
using namespace std;


#define COLOR_INTERVALS 15


vector<int> r_vals;
vector<int> g_vals;
vector<int> b_vals;

int last_r, last_g, last_b;

void addColor(int r, int g, int b)
{
  r_vals.push_back(last_r = r);
  g_vals.push_back(last_g = g);
  b_vals.push_back(last_b = b);
}

void makeColors(int set_r, int set_g, int set_b)
{
  double start_r = (double)last_r + 0.5;
  double start_g = (double)last_g + 0.5;
  double start_b = (double)last_b + 0.5;
  double r_step = (double)(set_r * 255 - last_r) / (double)(COLOR_INTERVALS);
  double g_step = (double)(set_g * 255 - last_g) / (double)(COLOR_INTERVALS);
  double b_step = (double)(set_b * 255 - last_b) / (double)(COLOR_INTERVALS);

  for (int i = 1; i <= COLOR_INTERVALS; ++i)
    addColor((int)(start_r + r_step * (double)i), (int)(start_g + g_step * (double)i), (int)(start_b + b_step * (double)i));
}

void makeColors()
{
  addColor(0, 0, 255);
  makeColors(0, 1, 1);
  makeColors(0, 1, 0);
  makeColors(1, 1, 0);
  makeColors(1, 0, 0);
}


int main( int argc, char* argv[] )
{
  GoTools::init();

  if (argc < 4 || argc > 6) {
    cout << "Usage:  " << argv[0] << " surfaceFile pointFile use_id_reg [search_extend, default = 3] [diff_color_file, absent = no file]" << endl;
    return 1;
  }

  ifstream in_surf(argv[1]);
  ObjectHeader header;
  vector<shared_ptr<GeomObject> > surfaces;

  while (!in_surf.eof())
    {
      header.read(in_surf);
      shared_ptr<GeomObject> obj(Factory::createObject(header.classType()));
      obj->read(in_surf);
      surfaces.push_back(obj);
      Utils::eatwhite(in_surf);
    }
  in_surf.close();

  ifstream in_pts(argv[2]);
  vector<float> pts;

  while (!in_pts.eof())
    {
      for (int j = 0; j < 3; ++j)
	{
	  float f;
	  in_pts >> f;
	  pts.push_back(f);
	}
      Utils::eatwhite(in_pts);
    }
  in_pts.close();

  bool use_id_reg = atoi(argv[3]) == 1;
  int search_extend = 3;
  if (argc >= 5)
    search_extend = atoi(argv[4]);

  vector<vector<double> > regRotation;
  Point regTranslation;

  if (use_id_reg)
    {
      regTranslation = Point(0.0, 0.0, 0.0);
      regRotation.resize(3);
      for(int i = 0; i < 3; i++)
	{
	  regRotation[i].resize(3);
	  for (int j = 0; j < 3; ++j)
	    regRotation[i][j] = i==j;
	}
    }
  else
    {
      vector<Point> reg_pts_pc;
      vector<Point> reg_pts_sf;

      reg_pts_pc.push_back(Point(535.215, -516.438, -561.635));
      reg_pts_sf.push_back(Point(526.994, -461.239, -540.883));

      reg_pts_pc.push_back(Point(1419.07, -1378.61, -513.867));
      reg_pts_sf.push_back(Point(1533.7, -1447.84, -523.269));

      reg_pts_pc.push_back(Point(1899.35, 762.057, 345.142));
      reg_pts_sf.push_back(Point(1961.96, 775.782, 365.377));

      RegistrationInput regParameters;
      RegistrationResult regResult = registration(reg_pts_sf, reg_pts_pc, false, regParameters);
      /*
      RegistrationResult rescaleRegResult = registration(reg_pts_sf, reg_pts_pc, true, regParameters);
      cout << "Rescaling set to " << rescaleRegResult.rescaling_ << endl;
      */
      regRotation = regResult.rotation_matrix_;
      regTranslation = regResult.translation_;
    }

  shared_ptr<boxStructuring::BoundingBoxStructure> structure = preProcessClosestVectors(surfaces, 200.0);

  // 0 = all, 1 = every 100 starting at 0, 2 = every 10 starting at 0, 3 = special,
  // 4 = compare different calls, 5 = drop signed distences
  // All + 8 = same for old
  // int round_type = 0;
  int round_type = 0;

  int my_round_type = round_type & 7;
  bool old_also = round_type >= 8;

  vector<float> distances;
  int beyond = 3000000;

  if (my_round_type == 0)
    {
      // All
      distances = closestPointCalculations(pts, structure, regRotation, regTranslation, 0, 0, 1, beyond, search_extend);
      if (old_also)
	distances = closestVectorsOld(pts, structure, regRotation, regTranslation, 0, 0, 1, beyond, search_extend);
    }
  else if (my_round_type == 1)
    {
      // Every 100
      distances = closestPointCalculations(pts, structure, regRotation, regTranslation, 0, 0, 100, beyond, search_extend);
      if (old_also)
	distances = closestVectorsOld(pts, structure, regRotation, regTranslation, 0, 0, 100, beyond, search_extend);
    }
  else if (my_round_type == 2)
    {
      // Every 10
      distances = closestPointCalculations(pts, structure, regRotation, regTranslation, 0, 0, 10, beyond, search_extend);
      if (old_also)
	distances = closestVectorsOld(pts, structure, regRotation, regTranslation, 0, 0, 10, beyond, search_extend);
    }
  else if (my_round_type == 3)
    {
      // Special, change at will
      distances = closestPointCalculations(pts, structure, regRotation, regTranslation, 0, 100, beyond, beyond, search_extend);
      if (old_also)
	distances = closestVectorsOld(pts, structure, regRotation, regTranslation, 0, 0, 100, 20, search_extend);
    }
  else if (my_round_type == 4)
    {
      // Compare different calls
      cout << "Running distance calculations" << endl;
      distances = closestDistances(pts, structure, regRotation, regTranslation);
      cout << "Running signed distance calculations" << endl;
      vector<float> signedDistances = closestSignedDistances(pts, structure, regRotation, regTranslation);
      cout << "Running points calculations" << endl;
      vector<float> clPoints = closestPoints(pts, structure, regRotation, regTranslation);

      if (3 * distances.size() != pts.size())
	{
	  cout << "3 x Size of distances is " << (3 * distances.size()) << " but should be " << pts.size() << endl;
	  exit(1);
	}

      if (3 * signedDistances.size() != pts.size())
	{
	  cout << "3 x Size of signedDistances is " << (3 * signedDistances.size()) << " but should be " << pts.size() << endl;
	  exit(1);
	}

      if (clPoints.size() != pts.size())
	{
	  cout << "Size of clPoints is " << clPoints.size() << " but should be " << pts.size() << endl;
	  exit(1);
	}

      double tol = 1.0e-4;

      for (int i = 0; i < distances.size(); ++i)
	if (abs(distances[i] - abs(signedDistances[i])) > tol)
	  {
	    cout << "Error: distances[" << i << "] = " << distances[i] << " while signedDistances[" << i << "] = " << signedDistances[i] << endl;
	    exit(1);
	  }

      for (int i = 0; i < distances.size(); ++i)
	{
	  double dist2 = 0.0;
	  for (int j = 0; j < 3; ++j)
	    {
	      double diff = pts[3*i + j] - clPoints[3*i + j];
	      dist2 += diff * diff;
	    }
	  double dist = sqrt(dist2);
	  if (abs(distances[i] - dist) > tol)
	    {
	      cout << "Error: distances[" << i << "] = " << distances[i] << endl;
	      cout << "       pts[" << i << "] = (";
	      for (int j = 0; j < 3; ++j)
		{
		  if (j > 0)
		    cout << ", ";
		  cout << pts[3*i + j];
		}
	      cout << ")" << endl;
	      cout << "       clPoints[" << i << "] = (";
	      for (int j = 0; j < 3; ++j)
		{
		  if (j > 0)
		    cout << ", ";
		  cout << clPoints[3*i + j];
		}
	      cout << ")" << endl;
	      cout << "       pts[" << i << "]-clPoints[" << i << "] = (";
	      for (int j = 0; j < 3; ++j)
		{
		  if (j > 0)
		    cout << ", ";
		  cout << (pts[3*i + j] - clPoints[3*i + j]);
		}
	      cout << ")" << endl;
	      cout << "       having length " << dist << " diffrence is to distances[" << i << "] is " << abs(distances[i] - dist) << endl;
	      exit(1);
	    }
	}

      cout << "All tests sucessfull. Now writing to file" << endl;

      vector<vector<int> > inout_pts(4);
      for (int i = 0; i < signedDistances.size(); ++i)
	inout_pts[((signedDistances[i] >= 0.0) ? 0 : 1) + ((distances[i] >= 3.0) ? 2 : 0)].push_back(i);

      /*
      for (int j = 1; j <= 10; ++j)
	{
	  int dist_cnt;
	  double tol_plus = (double)j;
	  double tol_minus = -tol_plus;
	  dist_cnt = 0;
	  for (int i = 0; i < signedDistances.size(); ++i)
	    if (signedDistances[i] > tol_plus)
	      ++dist_cnt;
	  cout << "Signed distance over " << tol_plus << " gave count = " << dist_cnt << endl;
	  dist_cnt = 0;
	  for (int i = 0; i < signedDistances.size(); ++i)
	    if (signedDistances[i] < tol_minus)
	      ++dist_cnt;
	  cout << "Signed distance under " << tol_minus << " gave count = " << dist_cnt << endl;
	}
      */
      ofstream inout_str("in_out.g2");

      for (int i = 0; i < 4; ++i)
	{
	  if (i == 0)
	    inout_str << "400 1 0 4 0 255 0 255" << endl;
	  else if (i == 1)
	    inout_str << "400 1 0 4 255 0 0 255" << endl;
	  else if (i == 2)
	    inout_str << "400 1 0 4 0 255 128 255" << endl;
	  else
	    inout_str << "400 1 0 4 255 128 0 255" << endl;
	  inout_str << inout_pts[i].size() << endl;
	  for (int j = 0; j < inout_pts[i].size(); ++j)
	    {
		int pos = 3 * inout_pts[i][j];
		inout_str << pts[pos] << " " << pts[pos + 1] << " " << pts[pos + 2] << endl;
	    }
	}
      inout_str.close();

      cout << "Completed, " << (inout_pts[0].size()) << " points are inside, " << (inout_pts[1].size()) << " points are outside" << endl;
    }
  else if (my_round_type == 5)
    {
      // Drop signed distances to file
      cout << "Running signed distance calculations" << endl;
      vector<float> signedDistances = closestSignedDistances(pts, structure, regRotation, regTranslation);
      cout << "Signed distances calculations completed" << endl;

      ofstream sig_dist_str("signed_dist.txt");

      for (int i = 0; i < signedDistances.size(); ++i)
	sig_dist_str << signedDistances[i] << endl;
      sig_dist_str.close();
    }
  else
    cout << "No closestVector call" << endl;

  if (distances.size() > 0 && argc == 6)
    {
      makeColors();
      int cols = r_vals.size();
      vector<vector<int> > color_group(cols);

      double worst_distance = -1.0;
      for (int i = 0; i < distances.size(); ++i)
	{
	  double dist = sqrt(distances[i]);
	  if (worst_distance < dist)
	    worst_distance = dist;
	}
      double factor = 1.33 * (double)cols / worst_distance;

      for (int i = 0; i < distances.size(); ++i)
	{
	  double dist = sqrt(distances[i]);
	  int group = (int)(dist * factor);
	  if (group < 0)
        group = 0;
	  if (group >= cols)
	    group = cols - 1;
	  color_group[group].push_back(i);
	}

      for (int i = 0; i < cols; ++i)
	cout << "Group " << i << " has " << (color_group[i].size()) << " points" << endl;

      ofstream out_pts(argv[5]);
      for (int i = 0; i < cols; ++i)
	if (color_group[i].size() > 0)
	  {
	    out_pts << "400 1 0 4 " << r_vals[i] << " " << g_vals[i] << " " << b_vals[i] << " 255" << endl;
	    out_pts << color_group[i].size() << endl;
	    for (int j = 0; j < color_group[i].size(); ++j)
	      {
		int pos = 3 * color_group[i][j];
		out_pts << pts[pos] << " " << pts[pos + 1] << " " << pts[pos + 2] << endl;
	      }
	  }
      out_pts.close();
    }
}
