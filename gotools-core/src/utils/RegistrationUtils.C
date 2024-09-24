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
#include "GoTools/creators/SolveCG.h"
#include "GoTools/geometry/ParamSurface.h"
#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/geometry/BoundedSurface.h"
#include "GoTools/utils/omp.h"

using namespace std;
using namespace Go;


namespace Go
{

  namespace matrix3DUtils
  {

    typedef vector<vector<double> > matrix3D;


    /// Return the 3x3 zero matrix
    matrix3D zeroMatrix()
    {
      matrix3D result(3);
      for (int i = 0; i < 3; ++i)
	result[i].resize(3);
      return result;
    }

    /// Return the 3x3 identity matrix
    matrix3D identity3D()
    {
      matrix3D result = zeroMatrix();
      result[0][0] = result[1][1] = result[2][2] = 1.0;
      return result;
    }

    /// Add the 3x3 matrix 'add' into the 3x3 matrix 'base'. This will change 'base'
    void addInMatrix(matrix3D& base, const matrix3D& add)
    {
      for (int i = 0; i < 3; ++i)
	for (int j = 0; j < 3; ++j)
	  base[i][j] += add[i][j];
    }

    /// Add a scalar into the 3x3 matrix 'base'. This will change 'base'
    void multiplyInScalar(matrix3D& base, double scalar)
    {
      for (int i = 0; i < 3; ++i)
	for (int j = 0; j < 3; ++j)
	  base[i][j] *= scalar;
    }

    /// Add the 3x3 matrix 'add' multiplied by a scalar into the 3x3 matrix 'base'. This will change 'base'
    void addInMatrixScalar(matrix3D& base, const matrix3D& add, double scalar)
    {
      for (int i = 0; i < 3; ++i)
	for (int j = 0; j < 3; ++j)
	  base[i][j] += add[i][j] * scalar;
    }

    /// Multiply two 3x3 matrices and return the result. The input matrices are not changed
    matrix3D multiply(const matrix3D& m1, const matrix3D& m2)
    {
      matrix3D result = zeroMatrix();
      for (int i = 0; i < 3; ++i)
	for (int j = 0; j < 3; ++j)
	  {
	    double sum = 0.0;
	    for (int k = 0; k < 3; ++k)
	      sum += m1[i][k] * m2[k][j];
	    result[i][j] = sum;
	  }
      return result;
    }

    /// Return the 3x3 matrix of the outer product of the 3-dimensional vectors p and q
    matrix3D tensorProduct(Point p, Point q)
    {
      matrix3D result = zeroMatrix();
      for (int i = 0; i < 3; ++i)
	for (int j = 0; j < 3; ++j)
	  result[i][j] = p[i] * q[j];
      return result;
    }

    /// Return pTq + qTp where T is the tensor product and p and q are 3-dimensional vectors
    matrix3D symmetricTensorProduct(Point p, Point q)
    {
      matrix3D mat1 = tensorProduct(p, q);
      matrix3D mat2 = tensorProduct(q, p);
      addInMatrix(mat1, mat2);
      return mat1;
    }

    /// Return the cross product matrix of a 3-dimensional vector
    matrix3D crossProductMatrix(Point p)
    {
      matrix3D result = zeroMatrix();
      result[1][2] = -p[0];
      result[2][0] = -p[1];
      result[0][1] = -p[2];
      result[2][1] = p[0];
      result[0][2] = p[1];
      result[1][0] = p[2];
      return result;
    }

    /// Return the rotation matrix representing rotation of a given angle around a given normal vector
    matrix3D rotationMatrix(Point normal_vector, double angle)
    {
      matrix3D result = zeroMatrix();
      double cos_angle = cos(angle);
      double sin_angle = sin(angle);
      addInMatrixScalar(result, identity3D(), cos_angle);
      addInMatrixScalar(result, tensorProduct(normal_vector, normal_vector), 1.0 - cos_angle);
      addInMatrixScalar(result, crossProductMatrix(normal_vector), sin_angle);
      return result;
    }

    /// Return the rotation matrix from a rotation vector R, where length of R is rotation angle, and direction of R is normal vector of rotation
    matrix3D rotationMatrix(Point r)
    {
      double r2 = r.length2();
      if (r2 == 0.0)
	return identity3D();
      else
	{
	  double angle = sqrt(r2);
	  return rotationMatrix(r / angle, angle);
	}
    }

    /// Return the result of applying a 3x3 matrix to a point
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

  }   // end namespace Go::matrix3DUtils


  using namespace Go::matrix3DUtils;


//===========================================================================
  RegistrationResult rawRegistration(const vector<Point>& points_fixed, const vector<Point>& points_transform, bool allow_rescaling, RegistrationInput params)
//===========================================================================
  {
    RegistrationResult result;

    int n_pts = (int)points_fixed.size();
    if (n_pts < 3 || n_pts != (int)points_transform.size())
      {
	if (n_pts < 3)
	  result.result_type_ = TooFewPoints;
	else
	  result.result_type_ = PointSetSizeDiff;
	return result;
      }

    // Try to find a good tripple of point pairs for raw registration
    int best_idx1, best_idx2, best_idx3;
    double best_sq_area = -1.0;

    vector<Point> vec_fix(n_pts - 1);
    vector<Point> vec_tr(n_pts - 1);

    if (n_pts < 150)
      {
	// Test all point tripples, this is O(n*n*n)
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
		    }
		}
	  }
      }
    else
      {
	// Too many points, use a O(n) method instead, not guaranteed to give best (or useful) tripple, but should be good enough

	// First point is the one furthest away from mass center of points
	Point mass_center(0.0, 0.0, 0.0);
	for (int i = 0; i < n_pts; ++i)
	  mass_center += points_fixed[i];
	mass_center /= (double)n_pts;

	double best_sq_dist = -1.0;
	for (int i = 0; i < n_pts; ++i)
	  {
	    double sq_dist = mass_center.dist2(points_fixed[i]);
	    if (sq_dist > best_sq_dist)
	      {
		best_sq_dist = sq_dist;
		best_idx1 = i;
	      }
	  }
	Point first_pt = points_fixed[best_idx1];

	// Second point is the one furthest away from the first
	best_sq_dist = -1.0;
	for (int i = 0; i < n_pts; ++i)
	  if (i != best_idx1)
	    {
	      double sq_dist = first_pt.dist2(points_fixed[i]);
	      if (sq_dist > best_sq_dist)
		{
		  best_sq_dist = sq_dist;
		  best_idx2 = i;
		}
	    }

	// Third point is the one spanning the biggest triangle together with the first two points,
	// both for the fixed points and the corresponding points in the changable point set
	Point first_pt_tr = points_transform[best_idx1];
	Point vec_fix = points_fixed[best_idx2] - first_pt;
	Point vec_tr = points_transform[best_idx2] - first_pt_tr;
	for (int i = 0; i < n_pts; ++i)
	  if (i != best_idx1 && i != best_idx2)
	    {
	      Point n_fix = (points_fixed[i] - first_pt) % vec_fix;
	      Point n_tr = (points_transform[i]- first_pt_tr) % vec_tr;
	      double area_prod = n_fix.length2() * n_tr.length2();
	      if (area_prod > best_sq_area)
		{
		  best_sq_area = area_prod;
		  best_idx3 = i;
		}
	    }
      }

    // It is required that the tripples are non-colinear on both sets
    if (best_sq_area <= 0.0)
      {
	result.result_type_ = AreaTooSmall;
	return result;
      }

    // Get point tripples
    Point q1 = points_fixed[best_idx1];
    Point q2 = points_fixed[best_idx2];
    Point q3 = points_fixed[best_idx3];
    Point p1 = points_transform[best_idx1];
    Point p2 = points_transform[best_idx2];
    Point p3 = points_transform[best_idx3];

    // Get normal vectors and test if points are far enough from being colinear
    Point edge_q21 = q2 - q1;
    Point edge_q31 = q3 - q1;
    Point edge_q32 = q3 - q2;
    Point edge_p21 = p2 - p1;
    Point edge_p31 = p3 - p1;
    Point edge_p32 = p3 - p2;

    double edge_q21_l2 = edge_q21.length2();
    double edge_q31_l2 = edge_q31.length2();
    double edge_q32_l2 = edge_q32.length2();
    double edge_p21_l2 = edge_p21.length2();
    double edge_p31_l2 = edge_p31.length2();
    double edge_p32_l2 = edge_p32.length2();

    double max_q_len_prod = max(edge_q21_l2 * edge_q31_l2, edge_q32_l2 * max(edge_q21_l2, edge_q31_l2));
    double max_p_len_prod = max(edge_p21_l2 * edge_p31_l2, edge_p32_l2 * max(edge_p21_l2, edge_p31_l2));

    Point norm_q = edge_q21 % edge_q31;
    Point norm_p = edge_p21 % edge_p31;

    if (norm_q.length2() < params.area_tolerance_sq_ * max_q_len_prod ||
	norm_p.length2() < params.area_tolerance_sq_ * max_p_len_prod)
      {
	result.result_type_ = AreaTooSmall;
	return result;
      }

    // Calculate rescaling factor without applying it on the tripple
    double raw_scale = 1.0;
    if (allow_rescaling)
      raw_scale = sqrt(sqrt(norm_q.length2() / norm_p.length2()));

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

    // Put together to create final raw registration
    result.rotation_matrix_ = multiply(rot_2, rot_1);
    result.rescaling_ = raw_scale;
    result.translation_ = center_q - apply(result.rotation_matrix_, center_p);
    result.result_type_ = RegistrationOK;
    return result;
  }


  void addToLinearSystem(int pt_idx, const vector<Point>& points_fixed, const vector<Point>& points_transform, bool allow_rescaling,
			 const vector<vector<double> >& id, const Point& fine_R, const Point& fine_T, double fine_s,
			 const vector<vector<double> >& m_rot_R, double s2, double R2, bool zero_R,
			 vector<vector<vector<double> > >& all_lhs_matrix,
			 vector<vector<double> >& all_rhs_matrix)
  {
#ifdef _OPENMP
    int thread_id = omp_get_thread_num();
#else
    int thread_id = 0;
#endif

    // Variables holding the contributions to the linear system
    Point dEdR(0.0, 0.0, 0.0);
    Point dEdT(0.0, 0.0, 0.0);
    double dEds = 0.0;
    matrix3D ddEdRdR = zeroMatrix();
    matrix3D ddEdRdT = zeroMatrix();
    matrix3D ddEdTdT = zeroMatrix();
    Point ddEdRds(0.0, 0.0, 0.0);
    Point ddEdTds(0.0, 0.0, 0.0);
    double ddEdsds = 0.0;

    // Some matrices, vectors and operations on these, used both for R = 0 and R != 0
    Point v_p = points_transform[pt_idx];
    Point v_q = points_fixed[pt_idx];

    Point v_t = apply(m_rot_R, v_p) + fine_T - v_q;
    Point v_U = v_t * 2.0 + v_q - fine_T;
    Point v_pxt = v_p % v_t;
    Point v_pxU = v_p % v_U;

    matrix3D p_sym_t = symmetricTensorProduct(v_p, v_t);
    matrix3D p_ten_p = tensorProduct(v_p, v_p);
    matrix3D p_cross = crossProductMatrix(v_p);

    double p2 = v_p * v_p;
    double p_dot_t = v_p * v_t;

    if (zero_R)
      {
	// Current fine rotation vector is zero
	dEdR += v_pxt * fine_s;
	dEdT += v_t;
	dEds += p_dot_t;

	addInMatrixScalar(ddEdRdR, p_sym_t, 0.5 * fine_s);
	addInMatrixScalar(ddEdRdR, p_ten_p, -s2);
	addInMatrixScalar(ddEdRdR, id, s2 * p2 - fine_s * p_dot_t);
	addInMatrixScalar(ddEdRdT, p_cross, fine_s);
	addInMatrix(ddEdTdT, id);
	ddEdRds += v_pxU;
	ddEdTds += v_p;
	ddEdsds += p2;
      }

    else
      {
	// Current fine rotation vector is non-zero

	// Calculate A, B and C-values
	double len_R = sqrt(R2);
	double sin_R = sin(len_R);
	double cos_R = cos(len_R);
	double one_min_cos = 1.0 - cos_R;

	double inv_R1 = 1.0 / len_R;
	double inv_R2 = 1.0 / R2;
	double inv_R3 = inv_R1 * inv_R2;
	double inv_R4 = inv_R2 * inv_R2;
	double inv_R5 = inv_R2 * inv_R3;
	double inv_R6 = inv_R2 * inv_R4;

	double A0 = cos_R;
	double B0 = one_min_cos * inv_R2;
	double C0 = sin_R * inv_R1;
	double A1 = -sin_R * inv_R1;
	double B1 = (len_R * sin_R - 2.0 * one_min_cos) * inv_R4;
	double C1 = (len_R * cos_R - sin_R) * inv_R3;
	double A2 = (sin_R - len_R * cos_R) * inv_R3;
	double B2 = (R2 * cos_R - 5.0 * len_R * sin_R + 8.0 * one_min_cos) * inv_R6;
	double C2 = (3.0 * sin_R - 3.0 * len_R * cos_R - R2 * sin_R) * inv_R5;

	// Some matrices, vectors and operations on these, used only for R != 0
	Point v_pxR = v_p % fine_R;

	matrix3D R_ten_R = tensorProduct(fine_R, fine_R);
	matrix3D R_ten_p = tensorProduct(fine_R, v_p);
	matrix3D p_ten_R = tensorProduct(v_p, fine_R);
	matrix3D R_ten_Rxp = tensorProduct(fine_R, -v_pxR);
	matrix3D p_sym_R = symmetricTensorProduct(v_p, fine_R);
	matrix3D pxR_sym_R = symmetricTensorProduct(v_pxR, fine_R);
	matrix3D p_sym_pxR = symmetricTensorProduct(v_p, v_pxR);
	matrix3D R_sym_t = symmetricTensorProduct(fine_R, v_t);
	matrix3D p_sym_t = symmetricTensorProduct(v_p, v_t);
	matrix3D pxt_sym_R = symmetricTensorProduct(v_pxt, fine_R);

	double p_dot_R = v_p * fine_R;
	double R_dot_t = fine_R * v_t;
	double p_dot_U = v_p * v_U;
	double R_dot_U = fine_R * v_U;
	double pR2 = p_dot_R * p_dot_R;
	double pxR2 = v_pxR * v_pxR;
	double det_tRp = v_pxt * fine_R;
	double det_URp = v_pxU * fine_R;

	// Calculate dEdR
	double coef;
	coef = fine_s * (A1 * p_dot_t + B1 * p_dot_R * R_dot_t + C1 * det_tRp);
	dEdR += fine_R * coef;
	dEdR += v_p * (fine_s * B0 * R_dot_t);
	dEdR += v_t * (fine_s * B0 * p_dot_R);
	dEdR += v_pxt * (fine_s * C0);

	// Calculate dEdT
	dEdT += v_t;

	// Calculate dEds
	dEds += A0 * p_dot_t;
	dEds += B0 * p_dot_R * R_dot_t;
	dEds += C0 * det_tRp;

	// Calculate ddEdRdR

	// Id contribution to ddEdRdR
	coef = B0 * B0 * pR2;
	coef += C0 * C0 * p2;
	coef *= fine_s;
	coef += A1 * p_dot_t;
	coef += B1 * p_dot_R * R_dot_t;
	coef += C1 * det_tRp;
	coef *= fine_s;
	addInMatrixScalar(ddEdRdR, id, coef);

	// R_tensor_R contribution to ddEdRdR
	coef = A1 * A1 * p2;
	coef += 2.0 * A1 * B1 * pR2;
	coef += B1 * B1 * pR2 * R2;
	coef += 2.0 * B0 * B1 * pR2;
	coef += C1 * C1 * pxR2;
	coef += 2.0 * C0 * C1 * p2;
	coef *= fine_s;
	coef += A2 * p_dot_t;
	coef += B2 * p_dot_R * R_dot_t;
	coef += C2 * det_tRp;
	coef *= fine_s;
	addInMatrixScalar(ddEdRdR, R_ten_R, coef);

	// R_sym_p contribution to ddEdRdR
	coef = 2.0 * A1 * B0;
	coef += B0 * B1 * R2;
	coef += B0 * B0;
	coef -= C0 * C1;
	coef *= fine_s * p_dot_R;
	coef += B1 * R_dot_t;
	coef *= fine_s;
	addInMatrixScalar(ddEdRdR, p_sym_R, coef);

	// Other contributions to ddEdRdR
	addInMatrixScalar(ddEdRdR, p_ten_p, s2 * (B0 * B0 * R2 - C0 * C0));
	addInMatrixScalar(ddEdRdR, pxR_sym_R, s2 * p_dot_R * (B1 * C0 - B0 * C1));
	addInMatrixScalar(ddEdRdR, p_sym_pxR, s2 * B0 * C0);
	addInMatrixScalar(ddEdRdR, R_sym_t, fine_s * B1 * p_dot_R);
	addInMatrixScalar(ddEdRdR, p_sym_t, fine_s * B0);
	addInMatrixScalar(ddEdRdR, pxt_sym_R, fine_s * C1);

	// Calculate ddEdRdT
	addInMatrixScalar(ddEdRdT, id, fine_s * B0 * p_dot_R);
	addInMatrixScalar(ddEdRdT, R_ten_p, fine_s * A1);
	addInMatrixScalar(ddEdRdT, R_ten_R, fine_s * B1 * p_dot_R);
	addInMatrixScalar(ddEdRdT, p_ten_R, fine_s * B0);
	addInMatrixScalar(ddEdRdT, R_ten_Rxp, fine_s * C1);
	addInMatrixScalar(ddEdRdT, p_cross, fine_s * C0);

	// Calculate ddEdTdT
	addInMatrix(ddEdTdT, id);

	// Calculate ddEdRds
	ddEdRds += fine_R * (A1 * p_dot_U + B1 * p_dot_R * R_dot_U + C1 * det_URp);
	ddEdRds += v_p * (B0 * R_dot_U);
	ddEdRds += v_U * (B0 * p_dot_R);
	ddEdRds += v_pxU * C0;

	// Calculate ddEdTds
	ddEdTds += v_p * A0;
	ddEdTds += fine_R * (B0 * p_dot_R);
	ddEdTds -= v_pxR * C0;

	// Calculate ddEdsds
	ddEdsds += A0 * A0 * p2;
	ddEdsds += B0 * B0 * pR2 * R2;
	ddEdsds += C0 * C0 * pxR2;
	ddEdsds += 2.0 * A0 * B0 * pR2;
      }

    // Add contributions to the coefficients of the linear system
    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j)
	all_lhs_matrix[thread_id][i][j] += ddEdRdR[i][j];

    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j)
	{
	  all_lhs_matrix[thread_id][i][j+3] += ddEdRdT[i][j];
	  all_lhs_matrix[thread_id][j+3][i] += ddEdRdT[i][j];
	}

    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j)
	all_lhs_matrix[thread_id][i+3][j+3] += ddEdTdT[i][j];

    for (int i = 0; i < 3; ++i)
      {
	all_lhs_matrix[thread_id][i][6] += ddEdRds[i];
	all_lhs_matrix[thread_id][6][i] += ddEdRds[i];
      }

    for (int i = 0; i < 3; ++i)
      {
	all_lhs_matrix[thread_id][i+3][6] += ddEdTds[i];
	all_lhs_matrix[thread_id][6][i+3] += ddEdTds[i];
      }

    all_lhs_matrix[thread_id][6][6] += ddEdsds;

    for (int i = 0; i < 3; ++i)
      all_rhs_matrix[thread_id][i] += dEdR[i];

    for (int i = 0; i < 3; ++i)
      all_rhs_matrix[thread_id][i+3] += dEdT[i];

    all_rhs_matrix[thread_id][6] += dEds;
   }


//===========================================================================
  RegistrationResult fineRegistration(const vector<Point>& points_fixed, const vector<Point>& points_transform, bool allow_rescaling, RegistrationInput params)
//===========================================================================
  {
    RegistrationResult result;

#ifdef _OPENMP
    int max_threads = omp_get_max_threads();
#else
    int max_threads = 1;
#endif

    int n_pts = (int)points_fixed.size();
    if (n_pts < 3 || n_pts != (int)points_transform.size())
      {
	if (n_pts < 3)
	  result.result_type_ = TooFewPoints;
	else
	  result.result_type_ = PointSetSizeDiff;
	return result;
      }

    double tol_2 = params.newton_tolerance_;
    int max_iterations = params.max_newton_iterations_;

    // Description of fine registration. Originally the identity operation
    // For description of the formulas, see own document
    Point fine_R(0.0, 0.0, 0.0);
    Point fine_T(0.0, 0.0, 0.0);
    double fine_s = 1.0;

    result.last_change_ = 0.0;
    matrix3D id = identity3D();

    for (int iteration = 0; iteration < max_iterations && (iteration == 0 || result.last_change_ >= tol_2); ++iteration)
      {
	result.last_newton_iteration_ = iteration;

	// Coefficients for linear system
	vector<vector<vector<double> > > all_lhs_matrix(max_threads);
	vector<vector<double> > all_rhs_matrix(max_threads);
	for (int th = 0; th < max_threads; ++th)
	  {
	    all_lhs_matrix[th].resize(7);
	    all_rhs_matrix[th].resize(7);
	    for (int i = 0; i < 7; ++i)
	      all_lhs_matrix[th][i].resize(7);
	  }

	// Calculations independent of each pair of points
	matrix3D m_rot_R = rotationMatrix(fine_R);
	multiplyInScalar(m_rot_R, fine_s);
	double s2 = fine_s * fine_s;
	double R2 = fine_R.length2();
	bool zero_R = R2 == 0.0;

#ifdef _OPENMP
	if (params.multi_core_)
	  {
	    // Run linear system calculations in multicore, because params.multi_core_=true and OPENMP is included
	    int pt_idx;
#pragma omp parallel \
  default(none)	\
  private(pt_idx) \
  shared(n_pts, points_fixed, points_transform, allow_rescaling, id, fine_R, fine_T, fine_s, m_rot_R, s2, R2, zero_R, all_lhs_matrix, all_rhs_matrix)
#pragma omp for OMP_SCHEDULE_AUTO
	    for (pt_idx = 0; pt_idx < n_pts; ++pt_idx)
	      addToLinearSystem(pt_idx, points_fixed, points_transform, allow_rescaling,
			    id, fine_R, fine_T, fine_s,
			    m_rot_R, s2, R2, zero_R,
			    all_lhs_matrix,
			    all_rhs_matrix);
	  }

	else

	  {
	    // Run linear system calculations in one single thread, because params.multi_core_=false
	    for (int pt_idx = 0; pt_idx < n_pts; ++pt_idx)
	      addToLinearSystem(pt_idx, points_fixed, points_transform, allow_rescaling,
			    id, fine_R, fine_T, fine_s,
			    m_rot_R, s2, R2, zero_R,
			    all_lhs_matrix,
			    all_rhs_matrix);
	  }

#else   // #ifdef _OPENMP

	// Run linear system calculations in one single thread, because OPENMP is not included
	for (int pt_idx = 0; pt_idx < n_pts; ++pt_idx)
	  addToLinearSystem(pt_idx, points_fixed, points_transform, allow_rescaling,
	  id, fine_R, fine_T, fine_s,
	  m_rot_R, s2, R2, zero_R,
	  all_lhs_matrix,
	  all_rhs_matrix);

#endif   // #ifdef _OPENMP

	vector<vector<double> > lhs_matrix(7);
	vector<double> rhs_matrix(7);
	for (int i = 0; i < 7; ++i)
	  {
	    lhs_matrix[i].resize(7);
	    for (int j = 0; j < 7; ++j)
	      {
		for (int th = 0; th < max_threads; ++th)
		  lhs_matrix[i][j] += all_lhs_matrix[th][i][j];
	      }
	    for (int th = 0; th < max_threads; ++th)
	      rhs_matrix[i] += all_rhs_matrix[th][i];
	  }

	// Make flat representation of the 6x6 (rescaling not allowed) or 7x7 (rescaling allowed) Hessian matrix
	vector<double> lhs_matrix_flat;
	int n_rows = allow_rescaling ? 7 : 6;
	for (int i = 0; i < n_rows; ++i)
	  for (int j = 0; j < n_rows; ++j)
	    lhs_matrix_flat.push_back(lhs_matrix[i][j]);

	// Solve the equation system by Conjugate Gradient Method.
	SolveCG solveCg;
	solveCg.attachMatrix(&lhs_matrix_flat[0], n_rows);
	solveCg.setTolerance(params.solve_tolerance_);
	solveCg.setMaxIterations(params.max_solve_iterations_);
	solveCg.precondRILU(0.1);

	// Fetch solution
	vector<double> change(n_rows, 0.0);
	int solve_res = solveCg.solve(&change[0], &rhs_matrix[0], n_rows);
	if (solve_res < 0 || solve_res == 1)
	  {
	    result.result_type_ = SolveFailed;
	    result.solve_result_ = solve_res;
	    return result;
	  }

	// Calculate weigths for tolerance calculations
	Point change_R = Point(change[0], change[1], change[2]);
	Point change_T = Point(change[3], change[4], change[5]);
	if (iteration == 0 && params.calculate_tolerance_weights_)
	  {
	    double len_R = change_R.length2();
	    if (len_R == 0.0)
	      len_R = 1.0;
	    double len_T = change_T.length2();
	    if (len_T == 0.0)
	      len_T = 1.0;
	    params.tolerance_weight_translation_ = 1.0;
	    params.tolerance_weight_rotation_ = len_T / len_R;
	    if (allow_rescaling)
	      {
		double len_s = change[6] * change[6];
		if (len_s == 0.0)
		  len_s = 1.0;
	      params.tolerance_weight_rescale_ = len_T / len_s;
	      }
	    else
	      params.tolerance_weight_rescale_ = 1.0;
	  }

	// Apply result of solution of equation system
	fine_R -= change_R;
	fine_T -= change_T;
	if (allow_rescaling)
	  fine_s -= change[6];

	result.last_change_ = change_R.length2() * params.tolerance_weight_rotation_ + change_T.length2() * params.tolerance_weight_translation_;
	if (allow_rescaling)
	  result.last_change_ += change[6] * change[6] * params.tolerance_weight_rescale_;

      }   // End Newton iteration

    if (result.last_change_ >= tol_2)
      result.last_newton_iteration_ = max_iterations;

    result.rotation_matrix_ = rotationMatrix(fine_R);
    result.rescaling_ = fine_s;
    result.translation_ = fine_T;
    result.result_type_ = RegistrationOK;

    return result;
  }


//===========================================================================
  RegistrationResult registration(const vector<Point>& points_fixed, const vector<Point>& points_transform, bool allow_rescaling, RegistrationInput params)
//===========================================================================
  {
    // First create raw registration
    RegistrationResult raw_result = rawRegistration(points_fixed, points_transform, allow_rescaling, params);
    if (raw_result.result_type_ != RegistrationOK)
      return raw_result;

    // Apply raw registration on second point set
    int n_pts = points_transform.size();
    vector<Point> points_raw_transformed(0);
    matrix3D rotation_and_scalar = zeroMatrix();
    addInMatrixScalar(rotation_and_scalar, raw_result.rotation_matrix_, raw_result.rescaling_);
    for (int i = 0; i < n_pts; ++i)
      points_raw_transformed.push_back(apply(rotation_and_scalar, points_transform[i]) + raw_result.translation_);

    // Create fine registration
    RegistrationResult fine_result = fineRegistration(points_fixed, points_raw_transformed, allow_rescaling, params);
    if (fine_result.result_type_ != RegistrationOK)
      return fine_result;

    // Finally, put together the two registrations
    RegistrationResult result;
    result.rotation_matrix_ = multiply(fine_result.rotation_matrix_, raw_result.rotation_matrix_);
    result.rescaling_ = fine_result.rescaling_ * raw_result.rescaling_;
    result.translation_ = apply(fine_result.rotation_matrix_, raw_result.translation_) * fine_result.rescaling_ + fine_result.translation_;
    result.last_newton_iteration_ = fine_result.last_newton_iteration_;
    result.last_change_ = fine_result.last_change_;
    result.result_type_ = fine_result.result_type_;
    return result;
  }


}   // end namespace Go
