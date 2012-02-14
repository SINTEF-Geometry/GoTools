//===========================================================================
//                                                                           
// File: SplineUtils.C                                                           
//                                                                           
// Created: Tue Oct 17 16:05:26 2000                                         
//                                                                           
// Author: Atgeirr F Rasmussen <atgeirr@sintef.no>
//                                                                           
// Revision: $Id: Utils.C,v 1.22 2005-10-03 08:15:53 sbr Exp $
//                                                                           
// Description:
//                                                                           
//===========================================================================


#include "GoTools/geometry/Utils.h"
#include "GoTools/geometry/SplineUtils.h"
#include "GoTools/utils/BaryCoordSystemTriangle3D.h"
#include <vector>
#include <fstream>

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


} // namespace Go

