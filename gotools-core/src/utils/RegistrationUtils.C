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

using namespace std;

#define LOGLINE cout<<"Line number " << __LINE__ << endl;

namespace Go
{

  namespace matrix3DUtils
  {

    typedef vector<vector<double> > matrix3D;


    void addInMatrix(matrix3D& base, const matrix3D& add)
    {
      for (int i = 0; i < 3; ++i)
	for (int j = 0; j < 3; ++j)
	  base[i][j] += add[i][j];
    }

    void multiplyInScalar(matrix3D& base, double scalar)
    {
      for (int i = 0; i < 3; ++i)
	for (int j = 0; j < 3; ++j)
	  base[i][j] *= scalar;
    }

    void addInMatrixScalar(matrix3D& base, const matrix3D& add, double scalar)
    {
      for (int i = 0; i < 3; ++i)
	for (int j = 0; j < 3; ++j)
	  base[i][j] += add[i][j] * scalar;
    }

    matrix3D multiply(const matrix3D& m1, const matrix3D& m2)
    {
      matrix3D result(3);
      for (int i = 0; i < 3; ++i)
	{
	  result[i].resize(3);
	  for (int j = 0; j < 3; ++j)
	    {
	      double sum = 0.0;
	      for (int k = 0; k < 3; ++k)
		sum += m1[i][k] * m2[k][j];
	      result[i][j] = sum;
	    }
	}
      return result;
    }

    /// Return the 3x3 matrix of the outer product of p and q
    matrix3D tensorProduct(Point p, Point q)
    {
      matrix3D result(3);
      for (int i = 0; i < 3; ++i)
	{
	  result[i].resize(3);
	  for (int j = 0; j < 3; ++j)
	    result[i][j] = p[i] * q[j];
	}
      return result;
    }

    /// Return pTq + qTp where T is the tensor product
    matrix3D symmetricTensorProduct(Point p, Point q)
    {
      matrix3D mat1 = tensorProduct(p, q);
      matrix3D mat2 = tensorProduct(q, p);
      addInMatrix(mat1, mat2);
      return mat1;
    }

    /// Return the cross product matrix of a vector
    matrix3D crossProductMatrix(Point p)
    {
      matrix3D result(3);
      for (int i = 0; i < 3; ++i)
	result[i].resize(3);
      result[0][0] = result[1][1] = result[2][2] = 0.0;
      result[0][1] = -p[2];
      result[0][2] = p[1];
      result[1][0] = p[2];
      result[1][2] = -p[0];
      result[2][0] = -p[1];
      result[2][1] = p[0];
      return result;
    }

    /// Return the 3x3 identity matrix
    matrix3D identity3D()
    {
      matrix3D result(3);
      for (int i = 0; i < 3; ++i)
	result[i].resize(3);
      result[0][0] = result[1][1] = result[2][2] = 1.0;
      return result;
    }

    /// Apply a 3x3 matrix on a point
    Point apply(const matrix3D& mat, const Point p)
    {
      Point q(3);
      for (int i = 0; i < 3; ++i)
	{
	  q[i] = 0.0;
	  for (int j = 0; j < 3; ++j)
	    q[i] += p[j] * mat[i][j];
	}
      return q;
    }

    /// REMOVE LATER
    void dropMatrix(const matrix3D& mat)
    {
      for (int i = 0; i < 3; ++i)
	{
	  if (i == 0) cout << "["; else cout << " ";
	  cout << "[" << mat[i][0] << " " << mat[i][1] << " " << mat[i][2] << "]";
	  if (i == 2) cout << "]";
	  cout << endl;
	}
    }

  }   // end namespace Go::matrix3DUtils


  using namespace Go::matrix3DUtils;


//===========================================================================
  bool registration(const vector<Point>& points_fixed, const vector<Point>& points_transform,
		    vector<vector<double> >& rotation_matrix, Point& translate)
//===========================================================================
  {
    int n_pts = (int)points_fixed.size();
    if (n_pts < 3 || n_pts != (int)points_transform.size())
      return false;

    // Find best tripple of point pairs for raw registration
    int best_idx1, best_idx2, best_idx3;
    Point norm_fix, norm_tr;
    double best_sq_area = -1.0;

    vector<Point> vec_fix(n_pts - 1);
    vector<Point> vec_tr(n_pts - 1);

    for (int i1 = 0; i1 < n_pts - 2; ++i1)
      {
	for (int idx = 0, j = i1 + 1; j < n_pts; ++j, ++idx)
	  {
	    vec_fix[idx] = points_fixed[j] - points_fixed[0];
	    vec_tr[idx] = points_transform[j] - points_transform[0];
	  }
	for (int i2 = i1 + 1; i2 < n_pts - 1; ++i2)
	  for (int i3 = i2 + 1; i3 < n_pts; ++i3)
	    {
	      Point n_fix = vec_fix[i3 - (i1+1)] % vec_fix[i2 - (i1+1)];
	      Point n_tr = vec_tr[i3 - (i1+1)] % vec_tr[i2 - (i1+1)];
	      double area_prod = n_fix.length2() * n_tr.length2();
	      if (area_prod > best_sq_area)
		{
		  best_idx1 = i1;
		  best_idx2 = i2;
		  best_idx3 = i3;
		  best_sq_area = area_prod;
		  norm_fix = n_fix;
		  norm_tr = n_tr;
		}
	    }
      }

    // Some tripple must be non-colinear on both sets if we shall make a registration
    if (best_sq_area <= 0.0)
      return false;

    // Get point tripples
    Point q1 = points_fixed[best_idx1];
    Point q2 = points_fixed[best_idx2];
    Point q3 = points_fixed[best_idx3];
    Point p1 = points_transform[best_idx1];
    Point p2 = points_transform[best_idx2];
    Point p3 = points_transform[best_idx3];

    // Calculate rescaling factor without applying it on the tripple
    Point norm_q = (q2 - q1) % (q3 - q1);
    Point norm_p = (p2 - p1) % (p3 - p1);
    double raw_scale = sqrt(sqrt(norm_q.length2() / norm_p.length2()));

    // Translation vector for both tripples, to have center in origo
    Point center_q = (q1 + q2 + q3) / 3.0;
    Point center_p = (p1 + p2 + p3) / 3.0;

    // Rotation matrix rot_1, assures the two normal vectors coincide after the p_i are sent to rot_1*p_i
    
    matrix3D id = identity3D();
    matrix3D rot_1;

    norm_q.normalize();
    norm_p.normalize();
    double np_d_nq = norm_p * norm_q;

    if (np_d_nq == -1.0)
      {
	// The normal vectors point in opposite directions.
	Point best_unit(0.0, 0.0, 0.0);
	if (abs(norm_q[0]) <= abs(norm_q[1]))
	  {
	    if (abs(norm_q[0]) <= abs(norm_q[2]))
	      best_unit[0] = 1.0;
	    else
	      best_unit[2] = 1.0;
	  }
	else
	  {
	    if (abs(norm_q[1]) <= abs(norm_q[2]))
	      best_unit[1] = 1.0;
	    else
	      best_unit[2] = 1.0;
	  }
	Point n = norm_q % best_unit;
	n.normalize();
	rot_1 = tensorProduct(n, n);
	multiplyInScalar(rot_1, 2.0);
	addInMatrixScalar(rot_1, id, -1.0);
      }

    else

      {

	Point np_x_nq = norm_p % norm_q;

	rot_1 = tensorProduct(np_x_nq, np_x_nq);
	multiplyInScalar(rot_1, 1.0 / (1.0 + np_d_nq));
	addInMatrixScalar(rot_1, id, np_d_nq);
	addInMatrix(rot_1, crossProductMatrix(np_x_nq));
      }

    // Rotation matrix rot_2, assures that p1 and q1 point in same direction after rot_2 is applied to p1

    q1 -= center_q;
    p1 -= center_p;
    p1 = apply(rot_1, p1);

    double inv_norms_p1_q1 = 1.0 / (p1.length() * q1.length());
    double cos_alpha = (p1 * q1) * inv_norms_p1_q1;
    double sin_alpha = norm_q * (p1 % q1) * inv_norms_p1_q1;

    matrix3D rot_2 = tensorProduct(norm_q, norm_q);
    multiplyInScalar(rot_2, 1.0 - cos_alpha);
    addInMatrixScalar(rot_2, id, cos_alpha);
    addInMatrixScalar(rot_2, crossProductMatrix(norm_q), sin_alpha);

    // Returning the raw registration only
    rotation_matrix = multiply(rot_2, rot_1);
    multiplyInScalar(rotation_matrix, raw_scale);
    translate = center_q - apply(rotation_matrix, center_p);
    return true;
  }

}   // end namespace Go
