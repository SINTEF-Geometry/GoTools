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

#include "GoTools/geometry/Utils.h"
#include "GoTools/geometry/SplineUtils.h"
#include "GoTools/utils/BaryCoordSystemTriangle3D.h"
#include <vector>
#include <fstream>
#include <assert.h>

using namespace std;

namespace Go
{


//===========================================================================
void SplineUtils::transpose_array(int dim, int num_old_rows, int num_old_cols, 
		     double* array_start)
//===========================================================================
{
    std::vector<double> c(array_start,
			  array_start + dim*num_old_cols*num_old_rows);
    for (int new_row = 0; new_row < num_old_cols; ++new_row) {
	for (int new_col = 0; new_col < num_old_rows; ++new_col) {
	    for (int dd = 0; dd < dim; ++dd) {
		array_start[new_row*num_old_rows*dim + new_col*dim + dd]
		    = c[new_col*num_old_cols*dim + new_row*dim + dd];
	    }
	}
    }
}



//===========================================================================
int SplineUtils::closest_in_array(const double* pt, const double* array, int n, int dim)
//===========================================================================
{
    double best_dist2 = Utils::distance_squared(pt, pt+dim, array);
    int best_index = 0;
    double dist2;
    for (int i = 1; i < n; ++i) {
        dist2 = Utils::distance_squared(pt, pt+dim, array + i*dim);
        if (dist2 < best_dist2) {
            best_dist2 = dist2;
            best_index = i;
        }
    }
    return best_index;
}


//===========================================================================
void SplineUtils::closest_on_rectgrid(const double* pt, const double* array,
			 int m, int n,
			 double& clo_u, double& clo_v)
//===========================================================================
{
    // m is number of columns
    // n is number of rows
    // in array, the column index runs fastest.

    // For each (not necessarily planar) quadrangle, split into 2
    // triangles. Then find distance to triangle, and closest point
    // in triangle.

    Vector3D pnt(pt);
    Vector3D p[4];
    Vector3D tri[3];

    Vector3D umask1(0, 1, 0);
    Vector3D umask2(1, 1, 0);
    Vector3D vmask1(0, 0, 1);
    Vector3D vmask2(0, 1, 1);

    double bestdist2 = 1e100;
    double best_u = 0;
    double best_v = 0;

    for (int col = 0; col < m-1; ++col) {
	for (int row = 0; row < n-1; ++row) {
	    // Pick the corner points counterclockwise
	    p[0].setValue(array + (row*m + col)*3);
	    p[1].setValue(array + (row*m + col+1)*3);
	    p[2].setValue(array + ((row+1)*m + col+1)*3);
	    p[3].setValue(array + ((row+1)*m + col)*3);
	    // Lower triangle, points 0, 1, 3.
	    tri[0] = p[0];
	    tri[1] = p[1];
	    tri[2] = p[3];
	    double clo_dist2;
	    Vector3D cltri = SplineUtils::closest_on_triangle(pnt, tri, clo_dist2);
	    if (clo_dist2 < bestdist2) {
		best_u = col + umask1*cltri;
		best_v = row + vmask1*cltri;
		bestdist2 = clo_dist2;
	    }
	    // Upper triangle, points 1, 2, 3.
	    tri[0] = p[1];
	    tri[1] = p[2];
	    tri[2] = p[3];
	    cltri = SplineUtils::closest_on_triangle(pnt, tri, clo_dist2);
	    if (clo_dist2 < bestdist2) {
		best_u = col + umask2*cltri;
		best_v = row + vmask2*cltri;
		bestdist2 = clo_dist2;
	    }
	}
    }
    clo_u = best_u;
    clo_v = best_v;
    int i = int(floor(clo_u));
    if (i >= m-1)
	i = m-2;
    int j = int(floor(clo_v));
    if (j >= n-1)
	j = n-2;
    p[0].setValue(array + (j*m + i)*3);
    p[1].setValue(array + (j*m + i+1)*3);
    p[2].setValue(array + ((j+1)*m + i+1)*3);
    p[3].setValue(array + ((j+1)*m + i)*3);
}

//===========================================================================
void SplineUtils::closest_on_rectgrid(const double* pt, const double* array,
			 int u_min, int u_max, int v_min, int v_max,
			 int nmb_coefs_u,
			 double& clo_u, double& clo_v)
//===========================================================================
{
    // m is number of columns
    // n is number of rows
    // in array, the column index runs fastest.

    // For each (not necessarily planar) quadrangle, split into 2
    // triangles. Then find distance to triangle, and closest point
    // in triangle.

    Vector3D pnt(pt);
    Vector3D p[4];
    Vector3D tri[3];

    Vector3D umask1(0, 1, 0);
    Vector3D umask2(1, 1, 0);
    Vector3D vmask1(0, 0, 1);
    Vector3D vmask2(0, 1, 1);

    double bestdist2 = 1e100;
    double best_u = 0;
    double best_v = 0;

    for (int i = u_min; i < u_max; ++i) {
	for (int j = v_min; j < v_max; ++j) {
	    // Pick the corner points counterclockwise
	    p[0].setValue(array + (j*nmb_coefs_u + i)*3);
	    p[1].setValue(array + (j*nmb_coefs_u + i+1)*3);
	    p[2].setValue(array + ((j+1)*nmb_coefs_u + i+1)*3);
	    p[3].setValue(array + ((j+1)*nmb_coefs_u + i)*3);
	    // Lower triangle, points 0, 1, 3.
	    tri[0] = p[0];
	    tri[1] = p[1];
	    tri[2] = p[3];
	    double clo_dist2;
	    Vector3D cltri = SplineUtils::closest_on_triangle(pnt, tri, clo_dist2);
	    if (clo_dist2 < bestdist2) {
		best_u = i + umask1*cltri;
		best_v = j + vmask1*cltri;
		bestdist2 = clo_dist2;
	    }
	    // Upper triangle, points 1, 2, 3.
	    tri[0] = p[1];
	    tri[1] = p[2];
	    tri[2] = p[3];
	    cltri = SplineUtils::closest_on_triangle(pnt, tri, clo_dist2);
	    if (clo_dist2 < bestdist2) {
		best_u = i + umask2*cltri;
		best_v = j + vmask2*cltri;
		bestdist2 = clo_dist2;
	    }
	}
    }
    clo_u = best_u;
    clo_v = best_v;
    int ii = min(int(floor(clo_u)), u_max - 1);
    int jj = min(int(floor(clo_v)), v_max - 1);
    if (ii < 0)
      ii = floor(0.5*u_max);
    if (jj < 0)
      jj = floor(0.5*v_max);
    p[0].setValue(array + (jj*nmb_coefs_u + ii)*3);
    p[1].setValue(array + (jj*nmb_coefs_u + ii+1)*3);
    p[2].setValue(array + ((jj+1)*nmb_coefs_u + ii+1)*3);
    p[3].setValue(array + ((jj+1)*nmb_coefs_u + ii)*3);
}

//===========================================================================
Vector3D SplineUtils::closest_on_triangle(const Vector3D& pt,
			     const Vector3D tri[3],
			     double& clo_dist2)
//===========================================================================
{
    // Project the point to the triangle plane.
    Vector3D norm = (tri[1]-tri[0]) % (tri[2]-tri[0]);
    norm.normalize();
    Vector3D diff = pt-tri[0];
    Vector3D proj_pt = pt - (diff*norm)*norm;
    double proj_dist2 = (diff*norm)*(diff*norm);
    BaryCoordSystemTriangle3D cs(tri);
    Vector3D bc = cs.cartToBary(proj_pt);
    bool allpos = true;
    double bestdist2 = 1e100;
    Vector3D bestb;
    for (int i = 0; i < 3; ++i) {
	if (bc[i] < 0) {
	    allpos = false;
	    // Closest point may be on the segment given
	    // by the two other points.
	    double lp = SplineUtils::closest_on_line_segment(proj_pt,
						tri[(i+1)%3],
						tri[(i+2)%3]);
	    Vector3D b;
	    b[i] = 0;
	    b[(i+1)%3] = 1.0 - lp;
	    b[(i+2)%3] = lp;
	    Vector3D realpt = cs.baryToCart(b);
	    double dist2 = proj_pt.dist2(realpt);
	    if (dist2 < bestdist2) {
		bestb = b;
		bestdist2 = dist2;
	    }
	}
    }
    if (allpos) {
	bestb = bc;
	bestdist2 = 0;
    }
    clo_dist2 = bestdist2 + proj_dist2;
    return bestb;
}

//===========================================================================
double SplineUtils::closest_on_line_segment(const Vector3D& pt,
			       const Vector3D& beg,
			       const Vector3D& end)
//===========================================================================
{
    double t = ((pt-beg)*(end-beg))/beg.dist2(end);
    if (t < 0) return 0;
    if (t > 1) return 1;
    return t;
}

//===========================================================================
void SplineUtils::make_coef_array_from_rational_coefs(const double* from,
					 double* to,
					 int num_coefs,
					 int dim)
//===========================================================================
{
    for (int i = 0; i < num_coefs; ++i) {
	double w = from[dim];
	if (w != 1.0) {
	    for (int dd = 0; dd < dim; ++dd) {
		to[dd] = from[dd]/w;
	    }
	} else {
	    for (int dd = 0; dd < dim; ++dd) {
		to[dd] = from[dd];
	    }
	}
	from += dim+1;
	to += dim;
    }
}


//===========================================================================
void SplineUtils::splineToBezierTransfMat(const double* knots,
					  vector<double>& transf_mat)
//===========================================================================
{
    // We implement the algorithm given in:
    // "A generalized conversion matrix between non-uniform B-spline
    // and Bezier representations with applications in CAGD",
    // by Giulio Casciola, Lucia Romani (2004).
    // s^(m) is a (m+1)x(m+1) matrix. Given degree n, we start with
    // the 1x1 matrix identity matrix s^(0), then proceed to build a
    // (n+1)x(n+1) matrix. Each s^(m) matrix is built row by row.

    // Assuming cubic degree, i.e. the order is 4, giving matrix size
    // 16.
    if (transf_mat.size()!= 16)
	transf_mat.resize(16);

    // We loop through the the steps for constructing the elements.
    const int deg = 3; // Cubic case, i.e. degree 3.
    int kk = 3; // The pointer to the interval (tmin = knots[kk]).
    double tmin = knots[kk];
    double tmax = knots[kk+1];

    // We store the values used when building the row.
    vector<double> mat(16), mat_prev(16);
    mat_prev[0] = 1; // The initial value for s^(0).

    // We construct each row in the refinement matrix.
    // We're using the (13) formulation from the article (row procedure).
    double tpar, frac1, frac2;
    size_t ki, kj, kn;
    int ind_ki, ind1, ind2;
    // We compute the matrices s^(i) for i = 1, .., n.
    for (kn = 1; kn < deg + 1; ++kn) // s^1, ..., s^n.
    {   // Row based: Keeping each row ki fixed, varying kj.
	// Building a (kn+1)x(kn+1) matrix.
	for (ki = 0; ki < kn + 1; ++ki) // The rows.
	{
	    for (kj = 0; kj < kn + 1; ++kj) // The entries in the row.
	    {
		tpar = (ki == kn) ? tmax : tmin;
		frac1 = 1.0/(knots[kk+kj]-knots[kk+kj-kn]);
		frac2 = 1.0/(knots[kk+kj+1]-knots[kk+kj+1-kn]);
		ind_ki = (ki < kn) ? ki : kn - 1;
		ind1 = 4*ind_ki+kj-1;
		ind2 = 4*ind_ki+kj;
		mat[4*ki+kj] = (kj == 0) ?
		    (knots[kk+kj+1]-tpar)*frac2*mat_prev[ind2] :
		    ((kj == kn) ?
		     (tpar-knots[kk+kj-kn])*frac1*mat_prev[ind1] :
		     (tpar-knots[kk+kj-kn])*frac1*mat_prev[ind1] +
		     (knots[kk+kj+1]-tpar)*frac2*mat_prev[ind2]);
	    }
	}
	// @@sbr201206 Copying more than we need. but otherwise we
	// must juggle different matrix sizes and indexing.
	copy(mat.begin(), mat.end(), mat_prev.begin());
    }
    // We then copy the matrix values to the output matrix.
    copy(mat.begin(), mat.end(), transf_mat.begin());

    return;
}


//===========================================================================
void SplineUtils::extractBezierCoefs(const double* coefs,
				     const int num_coefs_u, const int num_coefs_v,
				     const int ind_u_min, const int ind_v_min,
				     const std::vector<double>& transf_mat_u,
				     const std::vector<double>& transf_mat_v,
				     std::vector<double>& bezier_coefs)
//===========================================================================
{
    const int dim = 3;
    const int cv_dim = dim*num_coefs_u; // Geometric dim 3, order_v = 4.
    const int order_u = 4;
    const int order_v = 4;
    int coef_ind_u = ind_u_min - order_u + 1;
    int coef_ind_v = ind_v_min - order_v + 1;
    size_t ki, kj, kk, kh;
    const double* coef_ptr = &coefs[dim*(coef_ind_v*num_coefs_u+coef_ind_u)];
    const double* ref_mat_ptr = &transf_mat_v[0];
    vector<double> ref_coefs_v(48, 0.0);
    for (ki = 0; ki < 4; ++ki) // Summing over the rows in the refinement matrix, in v-dir.
	for (kj = 0; kj < 4; ++kj) // The coefs in dir v in ref mat.
	    for (kk = 0; kk < 4; ++kk) // The dimension in dir u.
		for (kh = 0; kh < dim; ++kh)
		{
		    ref_coefs_v[ki*12+kk*dim+kh] += ref_mat_ptr[4*ki+kj]*coef_ptr[kj*cv_dim+kk*dim+kh];
		}

    // Then the u-dir (considering the sf a curve, with dim = dim*num_coefs_u).
    coef_ptr = &ref_coefs_v[0];
    ref_mat_ptr = &transf_mat_u[0];
    for (ki = 0; ki < 4; ++ki) // Summing over refinement matrix. Going in u-dir.
	for (kj = 0; kj < 4; ++kj) // The ref mat, going in the u-dir.
	    for (kk = 0; kk < 4; ++kk) // The coefs in dir v.
		for (kh = 0; kh < dim; ++kh) // Geometry dimension.
		{
		    bezier_coefs[12*kk+ki*dim+kh] += coef_ptr[kk*12+kj*dim+kh]*ref_mat_ptr[4*ki+kj];
		}


}


//===========================================================================
void SplineUtils::refinedBezierCoefsCubic(SplineSurface& spline_sf,
					  int ind_u_min, int ind_v_min,
					  vector<double>& bez_coefs)
//===========================================================================
{
    assert(!spline_sf.rational());

    if (bez_coefs.size() != 48)
	bez_coefs.resize(48);
    std::fill(bez_coefs.begin(), bez_coefs.end(), 0.0);

    // Values for inpute spline surface.
    int dim = spline_sf.dimension();
    int order_u = spline_sf.order_u();
    int order_v = spline_sf.order_u();
    int num_coefs_u = spline_sf.numCoefs_u();
    int num_coefs_v = spline_sf.numCoefs_v();

    // Checking that input index is within range.
    assert(ind_u_min >= order_u - 1 && ind_u_min < num_coefs_u);
    assert(ind_v_min >= order_v - 1 && ind_v_min < num_coefs_v);

    BsplineBasis& basis_u = spline_sf.basis_u();
    BsplineBasis& basis_v = spline_sf.basis_v();
    double* knot_u = &basis_u.begin()[0];
    double* knot_v = &basis_v.begin()[0];

    // We expect the knot index to refer to the last occurence.
    assert(knot_u[ind_u_min] != knot_u[ind_u_min+1]);
    assert(knot_v[ind_v_min] != knot_v[ind_v_min+1]);

    // We expect knot mult to be 1 or 4.
    int knot_mult_umin = (knot_u[ind_u_min-1] == knot_u[ind_u_min]) ? 4 : 1;
    int knot_mult_umax = (knot_u[ind_u_min+1] == knot_u[ind_u_min+2]) ? 4 : 1;
    int knot_mult_vmin = (knot_v[ind_v_min-1] == knot_v[ind_v_min]) ? 4 : 1;
    int knot_mult_vmax = (knot_v[ind_v_min+1] == knot_v[ind_v_min+2]) ? 4 : 1;

    bool kreg_at_ustart = (knot_mult_umin == 4);
    bool kreg_at_uend = (knot_mult_umax == 4);
    vector<double> transf_mat_u(16, 0.0);
    // if (!kreg_at_ustart && !kreg_at_uend)
    splineToBezierTransfMat(knot_u + ind_u_min - 3,
			    transf_mat_u);

#ifndef NDEBUG
    std::cout << "\ntransf_mat_u=" << std::endl;
    for (size_t kj = 0; kj < 4; ++kj)
    {
	for (size_t ki = 0; ki < 4; ++ki)
	    std::cout << transf_mat_u[kj*4+ki] << " ";
	std::cout << std::endl;
    }
    std::cout << std::endl;
#endif // NDEBUG

    // else
    // 	cubicTransfMat(knot_u + ind_u_min - 3,
    // 		       kreg_at_ustart, kreg_at_uend,
    // 		       transf_mat_u);

    bool kreg_at_vstart = (knot_mult_vmin == 4);
    bool kreg_at_vend = (knot_mult_vmax == 4);
    vector<double> transf_mat_v(16, 0.0);
    // if (!kreg_at_ustart && !kreg_at_uend)
    splineToBezierTransfMat(knot_v + ind_v_min - 3,
			    transf_mat_v);

#ifndef NDEBUG
    std::cout << "\ntransf_mat_v=" << std::endl;
    for (size_t kj = 0; kj < 4; ++kj)
    {
	for (size_t ki = 0; ki < 4; ++ki)
	    std::cout << transf_mat_v[kj*4+ki] << " ";
	std::cout << std::endl;
    }
    std::cout << std::endl;
#endif // NDEBUG

    extractBezierCoefs(&spline_sf.coefs_begin()[0],
		       num_coefs_u, num_coefs_v,
		       ind_u_min, ind_v_min,
		       transf_mat_u, transf_mat_v,
		       bez_coefs);

    return;
}


//===========================================================================
shared_ptr<SplineSurface> SplineUtils::refineToBezier(const SplineSurface& spline_sf)
//===========================================================================
{
    shared_ptr<SplineSurface> bez_sf;

    const BsplineBasis& bas_u = spline_sf.basis_u();
    const BsplineBasis& bas_v = spline_sf.basis_v();
    const int order_u = bas_u.order();
    const int order_v = bas_v.order();

    // We extract the unique knots.
    vector<double> new_knots_u, new_knots_v;
    // vector<double> ref_knots_u, ref_knots_v;
    vector<double>::const_iterator iter = bas_u.begin();
    while (iter != bas_u.end() - order_u)
    {
	if (iter[0] != iter[1])
	{
	    int knot_mult = bas_u.knotMultiplicity(iter[0]);
	    int num_insert = order_u - knot_mult;
	    if (num_insert > 0)
	    {
		new_knots_u.insert(new_knots_u.end(), num_insert, iter[0]);
	    }
	}
	++iter;
    }

    iter = bas_v.begin();
    while (iter != bas_v.end() - order_v)
    {
	if (iter[0] != iter[1])
	{
	    int knot_mult = bas_v.knotMultiplicity(iter[0]);
	    int num_insert = order_v - knot_mult;
	    if (num_insert > 0)
	    {
		new_knots_v.insert(new_knots_v.end(), num_insert, iter[0]);
	    }
	}
	++iter;
    }

    bez_sf = insertKnots(spline_sf,
			 new_knots_u, new_knots_v);

    return bez_sf;
}


//===========================================================================
shared_ptr<SplineSurface> GO_API SplineUtils::insertKnots(const Go::SplineSurface& spline_sf,
							  const std::vector<double> new_knots_u,
							  const std::vector<double> new_knots_v)
//===========================================================================
{
    shared_ptr<SplineSurface> ref_sf;
    int ki, kj, kk, kh, kl;

    const BsplineBasis& bas_u = spline_sf.basis_u();
    const BsplineBasis& bas_v = spline_sf.basis_v();
    const double* knots_u = &bas_u.begin()[0];
    const double* knots_v = &bas_v.begin()[0];
    const int order_u = bas_u.order();
    const int order_v = bas_v.order();
    const int num_coefs_u = bas_u.numCoefs();
    const int num_coefs_v = bas_v.numCoefs();

    vector<double> ref_knots_u(num_coefs_u + order_u + new_knots_u.size());
    copy(knots_u, knots_u + num_coefs_u + order_u, ref_knots_u.begin());
    copy(new_knots_u.begin(), new_knots_u.end(), ref_knots_u.begin() + num_coefs_u + order_u);
    sort(ref_knots_u.begin(), ref_knots_u.end());

    vector<double> ref_knots_v(num_coefs_v + order_v + new_knots_v.size());
    copy(knots_v, knots_v + num_coefs_v + order_v, ref_knots_v.begin());
    copy(new_knots_v.begin(), new_knots_v.end(), ref_knots_v.begin() + num_coefs_v + order_v);
    sort(ref_knots_v.begin(), ref_knots_v.end());

    const int num_coefs_u_ref = ref_knots_u.size() - order_u;
    const int num_coefs_v_ref = ref_knots_v.size() - order_v;

    vector<int> first_u(num_coefs_u_ref, -1), last_u(num_coefs_u_ref, -1);
    vector<double> ref_mat_u(num_coefs_u_ref*order_u);
    refmatrix(&ref_knots_u[0], num_coefs_u_ref, order_u, 
	      knots_u, num_coefs_u,
	      &ref_mat_u[0], &first_u[0], &last_u[0]);

    vector<int> first_v(num_coefs_v_ref, -1), last_v(num_coefs_v_ref, -1);
    vector<double> ref_mat_v(num_coefs_v_ref*order_v);
    refmatrix(&ref_knots_v[0], num_coefs_v_ref, order_v, 
	      knots_v, num_coefs_v,
	      &ref_mat_v[0], &first_v[0], &last_v[0]);

    // We refine line by line.
    const int dim = spline_sf.dimension();
    vector<double> coefs_ref_v(num_coefs_u*num_coefs_v_ref*dim, 0.0); // After refinement in the 2nd dir.
    // First in the v-dir.
    const int cv_dim = dim*num_coefs_u;
    const double* coef_ptr = &spline_sf.coefs_begin()[0];
    for (ki = 0; ki < num_coefs_v_ref; ++ki) // Summing over the rows in the refinement matrix, in v-dir.
	// The coefs in dir v in ref mat. order_v values (at most).
	for (kj = first_v[ki], kk = order_v - 1 - last_v[ki] + first_v[ki]; kj < last_v[ki] + 1; ++kj, ++kk)
	    for (kh = 0; kh < cv_dim; ++kh) // Considering the sf a curve, going in dir v.
	    {
		coefs_ref_v[ki*cv_dim+kh] += ref_mat_v[order_v*ki+kk]*coef_ptr[kj*cv_dim+kh];
	    }

#ifndef NDEBUG
    std::ofstream debug_out("tmp/sf_ref_v.g2");
    spline_sf.writeStandardHeader(debug_out);
    spline_sf.write(debug_out);
    vector<double> knot_vec_u(knots_u, knots_u + num_coefs_u + order_u);
    SplineSurface sf_ref_v(num_coefs_u, num_coefs_v_ref,
			   order_u, order_v,
			   knot_vec_u.begin(), ref_knots_v.begin(),
			   coefs_ref_v.begin(),
			   dim);
    sf_ref_v.writeStandardHeader(debug_out);
    sf_ref_v.write(debug_out);
#endif

    // Then in the u-dir.
    vector<double> coefs_ref_uv(num_coefs_u_ref*num_coefs_v_ref*dim, 0.0); // After refinement in both dirs.
    coef_ptr = &coefs_ref_v.begin()[0];
    for (ki = 0; ki < num_coefs_u_ref; ++ki) // Summing over the rows in the refinement matrix, in u-dir.
	// The coefs in dir u in ref mat. order_u values (at most).
	for (kj = first_u[ki], kk = order_u - 1 - last_u[ki] + first_u[ki]; kj < last_u[ki] + 1; ++kj, ++kk)
	    for (kh = 0; kh < num_coefs_v_ref; ++kh) // Considering the sf a curve, going in dir v.
		for (kl = 0; kl < dim; ++kl)
		{
		    coefs_ref_uv[(kh*num_coefs_u_ref+ki)*dim+kl] += ref_mat_u[order_u*ki+kk]*coef_ptr[(kh*num_coefs_u+kj)*dim+kl];
		}

    ref_sf = shared_ptr<SplineSurface>
	(new SplineSurface(num_coefs_u_ref, num_coefs_v_ref,
			   order_u, order_v,
			   ref_knots_u.begin(), ref_knots_v.begin(),
			   coefs_ref_uv.begin(),
			   dim));


#ifndef NDEBUG
    std::ofstream debug_out2("tmp/sf_ref_uv.g2");
    spline_sf.writeStandardHeader(debug_out2);
    spline_sf.write(debug_out2);
    ref_sf->writeStandardHeader(debug_out2);
    ref_sf->write(debug_out2);
#endif

    return ref_sf;

}


} // namespace Go

