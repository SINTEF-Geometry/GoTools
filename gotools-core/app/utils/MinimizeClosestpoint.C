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
//#include "GoTools/geometry/GeomObject.h"
#include "GoTools/geometry/Factory.h"
#include "GoTools/geometry/Utils.h"
#include "GoTools/utils/RegistrationUtils.h"
#include "GoTools/utils/ClosestPointUtils.h"



using namespace Go;
using namespace std;

// int show_pts = 10;
int show_pts = 0;


typedef pair<vector<vector<double> >, Point>  transformation_type;

transformation_type currentTransformation;

transformation_type combine(const transformation_type& first_transf, const transformation_type& second_transf)
{
  vector<vector<double> > rot(3);
  for (int i = 0; i < 3; ++i)
    {
      rot[i].resize(3);
      for (int j = 0; j < 3; ++j)
	{
	  double sum = 0.0;
	  for (int k = 0; k < 3; ++k)
	    sum += second_transf.first[i][k] * first_transf.first[k][j];
	  rot[i][j] = sum;
	}
    }

  Point transl(3);
  for (int i = 0; i < 3; ++i)
    {
      double sum = second_transf.second[i];
      for (int j = 0; j < 3; ++j)
	sum += second_transf.first[i][j] * first_transf.second[j];
      transl[i] = sum;
    }

  return transformation_type(rot, transl);
}


vector<Point> floatsToPoints(const vector<float>& pts_in)
{
  vector<Point> result(pts_in.size() / 3);
  for (int i = 0, idx = 0; i < pts_in.size(); i += 3, ++idx)
    result[idx] = Point(pts_in[i], pts_in[i + 1], pts_in[i + 2]);
  return result;
}


vector<Point> floatsToPoints(const vector<float>& pts_in, const transformation_type& transformation)
{
  vector<Point> result(pts_in.size() / 3);
  vector<vector<double> > rotation = transformation.first;
  Point translation = transformation.second;
  for (int i = 0, idx = 0; i < pts_in.size(); i += 3, ++idx)
    {
      if (idx < show_pts)
	{
	  cout << "  *** Applying tansformation on point " << idx << endl;
	  cout << "    before =";
	  for (int j = 0; j < 3; ++j)
	    cout << " " << pts_in[i+j];
	  cout << endl;
	}
      Point p(3);
      for (int j = 0; j < 3; ++j)
	{
	  double sum = translation[j];
	  for (int k = 0; k < 3; ++k)
	    sum += rotation[j][k] * pts_in[i + k];
	  p[j] = sum;
	}
      result[idx] = p;
      if (idx < show_pts)
	{
	  cout << "    after =";
	  for (int j = 0; j < 3; ++j)
	    cout << " " << p[j];
	  cout << endl;
	}
    }
  return result;
}


double transformationL2(const transformation_type& transformation)
{
  vector<vector<double> > rotation = transformation.first;
  Point translation = transformation.second;
  double sum2 = translation.length2();
  for (int i = 0; i < 3; ++i)
    {
      int next_i = (i + 1) % 3;
      double term = 0.5 * (rotation[i][next_i] - rotation[next_i][i]);
      sum2 += term * term;
    }
  return sum2;
}


double avgDist(const vector<float>& pts1, const vector<float>& pts2, const transformation_type& transformation)
{
  double sum2 = 0.0;
  vector<vector<double> > rotation = transformation.first;
  Point translation = transformation.second;
  int nmb_pts = 0;
  for (int i = 0; i < pts1.size(); i += 3, ++nmb_pts)
    {
      double prev_sum = sum2;
      if (nmb_pts < show_pts)
	{
	  cout << "  *** Dist calculations at " << nmb_pts << endl;
	  cout << "    clp =";
	  for (int j = 0; j < 3; ++j)
	    cout << " " << pts1[i+j];
	  cout << endl;
	  cout << "    pt =";
	  for (int j = 0; j < 3; ++j)
	    cout << " " << pts2[i+j];
	  cout << endl;
	  cout << "    T(pt) =";
	  for (int j = 0; j < 3; ++j)
	    {
	      double sum = translation[j];
	      for (int k = 0; k < 3; ++k)
		sum += rotation[j][k] * pts2[i + k];
	      cout << " " << sum;
	    }
	  cout << endl;
	  cout << "    T(pt) - clp =";
	}
      for (int j = 0; j < 3; ++j)
	{
	  double sum = translation[j] - pts1[i + j];
	  for (int k = 0; k < 3; ++k)
	    sum += rotation[j][k] * pts2[i + k];
	  sum2 += sum * sum;
	  if (nmb_pts < show_pts)
	    cout << " " << sum;
	}
      if (nmb_pts < show_pts)
	cout << endl << "    L2 is " << (sum2 - prev_sum) << endl;
    }
  return sum2 / (double)nmb_pts;
}


void drop_final(const vector<float>& pts1, const vector<float>& pts2, const transformation_type& transformation)
{
  ofstream out_str_1("transformed_points.txt");
  ofstream out_str_2("clp_distances.txt");

  vector<vector<double> > rotation = transformation.first;
  Point translation = transformation.second;

  for (int i = 0; i < pts1.size(); i += 3)
    {
      double sum2 = 0.0;
      for (int j = 0; j < 3; ++j)
	{
	  double sum = translation[j];
	  for (int k = 0; k < 3; ++k)
	    sum += rotation[j][k] * pts2[i + k];
	  out_str_1 << sum;
	  if (j < 2)
	    out_str_1 << " ";
	  else
	    out_str_1 << endl;
	  sum -= pts1[i + j];
	  sum2 += sum * sum;
	}
      out_str_2 << sqrt(sum2) << endl;
    }
  out_str_1.close();
  out_str_2.close();
}


void dropTransformation(const transformation_type& transformation, string s)
{
  // cout << "Dropping transformation:" << endl;
  cout << s << endl;
  vector<vector<double> > rotation = transformation.first;
  Point translation = transformation.second;
  for (int i = 0; i < 3; ++i)
    {
      for (int j = 0; j < 3; ++j)
	cout << "\t" << rotation[i][j];
      cout << "\t\t" << translation[i] << endl;
    }
  // cout << "Done dropping" << endl;
}


void registrationIteration(const vector<float>& pts, const shared_ptr<boxStructuring::BoundingBoxStructure>& structure, double changeL2tol)
{
  int nmb_pts = pts.size() / 3;
  for (int i = 0; i < 10000; ++i)
    {
      // cout << "Running iteration " << i << " for point set of size " << nmb_pts << endl;

      vector<float> clp = closestPoints(pts, structure, currentTransformation.first, currentTransformation.second);
      double avg_dist1 = avgDist(clp, pts, currentTransformation);

      RegistrationInput regParameters;
      vector<Point> clp_p = floatsToPoints(clp);
      vector<Point> pts_p = floatsToPoints(pts, currentTransformation);
      int max_newton_iterations = regParameters.max_newton_iterations_;
      RegistrationResult regResult = fineRegistration(clp_p, pts_p, false, regParameters);

      transformation_type changeTransformation = transformation_type(regResult.rotation_matrix_, regResult.translation_);
      double changeL2 = transformationL2(changeTransformation);
      currentTransformation = combine(currentTransformation, changeTransformation);
      double avg_dist2 = avgDist(clp, pts, currentTransformation);

      int newton_iterations = regResult.last_newton_iteration_;
      bool reg_OK = (regResult.result_type_ == RegistrationOK) && (newton_iterations < max_newton_iterations);

      cout << "Nmb pts: " << nmb_pts << " Iter: " << i;
      cout << " Avg clp dist: " << avg_dist1 << " => " << avg_dist2;
      cout << " Transf L2-chg: " << changeL2;
      cout << " Nmb Nwt it: " << (regResult.last_newton_iteration_ + 1);
      cout << " Transl:";
      for (int j = 0; j < 3; ++j)
	cout << " " << currentTransformation.second[j];
      cout << endl;
      if (!reg_OK)
	{
	  cout << endl << "******* REGISTRATION FAILED!!! *****" << endl;
	  cout << "  Newton method failing reason: ";
	  switch (regResult.result_type_)
	    {
	    case RegistrationOK:
	      cout << "Did not get close enough in maximum number of iterations (" << max_newton_iterations << ")" << endl;
	      break;
	    case TooFewPoints:
	      cout << "To few input points (must be at least 3, was " << clp.size() << ")" << endl;
	      break;
	    case PointSetSizeDiff:
	      cout << "Different size of input point sets (" << clp.size() << " vs " << pts.size() << ")" << endl;
	      break;
	    case SolveFailed:
	      cout << "Solving linear system failed" << endl;
	      break;
	    default:
	      cout << "Unknown reason (should not happen)" << endl;
	    }
	  exit(1);
	}

      if (changeL2 < changeL2tol)
	{
	  if (nmb_pts > 600000)
	    {
	      cout << "Dropping results to files" << endl;
	      drop_final(clp, pts, currentTransformation);
	    }
	  break;
	}

    }

}


int main( int argc, char* argv[] )
{
  GoTools::init();

  if (argc < 4 || argc > 5)
    {
      cout << "Usage:  " << argv[0] << " surfaceFile pointFile use_id_reg [run_reg? default=1]" << endl;
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

  vector<int> reduce_factors;
  vector<double> tolerances;

  bool use_id_reg = atoi(argv[3]) == 1;

  vector<vector<double> > startRotation;
  Point startTranslation;

  if (use_id_reg)
    {
      reduce_factors.push_back(100);
      tolerances.push_back(1.0e-5);
      reduce_factors.push_back(1);
      tolerances.push_back(1.0e-3);
      startTranslation = Point(0.0, 0.0, 0.0);
      startRotation.resize(3);
      for(int i = 0; i < 3; i++)
	{
	  startRotation[i].resize(3);
	  for (int j = 0; j < 3; ++j)
	    startRotation[i][j] = i==j;
	}
    }
  else
    {
      reduce_factors.push_back(1000);
      tolerances.push_back(1.0e-6);
      reduce_factors.push_back(100);
      tolerances.push_back(1.0e-5);
      reduce_factors.push_back(1);
      tolerances.push_back(1.0e-3);
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

      startRotation = regResult.rotation_matrix_;
      startTranslation = regResult.translation_;
    }

  currentTransformation = transformation_type(startRotation, startTranslation);
  dropTransformation(currentTransformation, "  Input rotation and transformation:");

  cout << "Preprocessing surface model..." << endl;
  shared_ptr<boxStructuring::BoundingBoxStructure> structure = preProcessClosestVectors(surfaces, 200.0);
  cout << "...done" << endl;

  bool perform_registration = (argc == 4 || atoi(argv[4]) == 1);

  if (perform_registration)
    {
      for (int i = 0; i < reduce_factors.size(); ++i)
	{
	  int red_fact = reduce_factors[i];
	  if (red_fact == 1)
	    registrationIteration(pts, structure, tolerances[i]);
	  else
	    {
	      vector<float> few_pts;
	      for (int j = 0, idx = 0; j < pts.size(); j += 3, ++idx)
		if ((idx % red_fact) == 0)
		  for (int k = 0; k < 3; ++k)
		    few_pts.push_back(pts[j + k]);
	      registrationIteration(few_pts, structure, tolerances[i]);
	    }
	}
      dropTransformation(currentTransformation, "  Final rotation and transformation:");
      cout << "Main diagonal entries diff from 1.0:";
      for (int i = 0; i < 3; ++i)
	cout << " " << (currentTransformation.first[i][i] - 1.0);
      cout << endl;
    }
  else
    {
      vector<float> clp = closestPoints(pts, structure, currentTransformation.first, currentTransformation.second);
      cout << "Dropping closest points for initial transformation to files" << endl;
      drop_final(clp, pts, currentTransformation);
    }

}
