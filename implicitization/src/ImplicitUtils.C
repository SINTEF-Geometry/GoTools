//===========================================================================
//                                                                           
// File: ImplicitUtils.C                                                
//                                                                           
// Created: Tue Apr 24 15:19:41 2001                                         
//                                                                           
// Authors: Atgeirr F Rasmussen <atgeirr@sintef.no>
//          Jan B. Thomassen <jbt@math.sintef.no>
//                                                                           
// Revision: $Id: ImplicitUtils.C,v 1.68 2006-09-19 09:23:14 jbt Exp $
//                                                                           
// Description:
//                                                                           
//===========================================================================


#include "GoTools/implicitization/ImplicitUtils.h"
#include "GoTools/implicitization/BernsteinPoly.h"
#include "GoTools/implicitization/BernsteinMulti.h"
#include "GoTools/implicitization/BernsteinTriangularPoly.h"
#include "GoTools/implicitization/BernsteinTetrahedralPoly.h"
#include "GoTools/implicitization/BernsteinUtils.h"
#include "GoTools/geometry/SplineCurve.h"
#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/utils/BaryCoordSystem.h"
#include "GoTools/utils/errormacros.h"
#include "newmat.h"
#include "newmatap.h"
//#include "newmatio.h"


using namespace std;
using namespace NEWMAT;


namespace Go {


//==========================================================================
void create_bary_coord_system2D(const SplineCurve& curve,
				BaryCoordSystem2D& bc)
//==========================================================================
{
    // We make a coordinate system based on a "regular" triangle. Its
    // properties are: 1) One corner is located at the lower left
    // corner of the bounding box, 2) The two edges going out from
    // this corner are axis parallell, and 3) It encloses the bounding
    // box.

    Point low = curve.boundingBox().low();
    Point high = curve.boundingBox().high();
    double len = (high[0] - low[0]) + (high[1] - low[1]);

    Vector2D corners[3];
    corners[0] = Vector2D(low[0], low[1]);
    corners[1] = Vector2D(low[0] + len, low[1]);
    corners[2] = Vector2D(low[0], low[1] + len);

    bc = BaryCoordSystem2D(corners);
    return;
}


//==========================================================================
void create_bary_coord_system3D(const SplineCurve& curve,
				BaryCoordSystem3D& bc)
//==========================================================================
{
    BoundingBox box = curve.boundingBox();
    create_bary_coord_system3D(box, bc);
    return;
}


//==========================================================================
void create_bary_coord_system3D(const SplineSurface& surface,
				BaryCoordSystem3D& bc)
//==========================================================================
{
    BoundingBox box = surface.boundingBox();
    create_bary_coord_system3D(box, bc);
    return;
}


//==========================================================================
void create_bary_coord_system3D(const PointCloud3D& cloud,
				BaryCoordSystem3D& bc)
//==========================================================================
{
    BoundingBox box = cloud.boundingBox();
    create_bary_coord_system3D(box, bc);
    return;
}


//==========================================================================
void create_bary_coord_system3D(const BoundingBox& box,
				BaryCoordSystem3D& bc)
//==========================================================================
{
    // We make a coordinate system based on a "regular"
    // tetrahedron. Its properties are: 1) One corner is located at
    // the "lower" corner of the bounding box, 2) The three edges
    // going out from this corner are axis parallell, and 3) It
    // encloses the bounding box.

    Point low = box.low();
    Point high = box.high();
    double len = (high[0] - low[0])
	+ (high[1] - low[1]) + (high[2] - low[2]);

    Vector3D corners[4];
    corners[0] = Vector3D(low[0], low[1], low[2]);
    corners[1] = Vector3D(low[0] + len, low[1], low[2]);
    corners[2] = Vector3D(low[0], low[1] + len, low[2]);
    corners[3] = Vector3D(low[0], low[1], low[2] + len);

    bc = BaryCoordSystem3D(corners);
    return;
}


//==========================================================================
void cart_to_bary(const SplineCurve& cv, const BaryCoordSystem2D& bc,
		  SplineCurve& cv_bc)
//==========================================================================
{
    ALWAYS_ERROR_IF(cv.dimension() != 2, "Dimension must be 2.");

    int n = cv.numCoefs();
    Vector2D cart;
    Vector3D bary;
    vector<double> new_coefs;
    if (!cv.rational()) {
	new_coefs.resize(3 * n);
	for (int i = 0; i < n; ++i) {
	    cart = Vector2D(cv.coefs_begin() + 2*i);
	    bary = bc.cartToBary(cart);
	    for (int j = 0; j < 3; ++j) {
		new_coefs[3*i + j] = bary[j];
	    }
	}
    } else {
	new_coefs.resize(4 * n);
	for (int i = 0; i < n; ++i) {
	    cart = Vector2D(cv.coefs_begin() + 2*i);
	    bary = bc.cartToBary(cart);
	    double w = cv.rcoefs_begin()[3*i + 2];
	    for (int j = 0; j < 3; ++j) {
		new_coefs[4*i + j] = bary[j] * w;
	    }
	    new_coefs[4*i + 3] = w;
	}
    }
    cv_bc = SplineCurve(n, cv.order(), cv.basis().begin(),
			new_coefs.begin(), 3, cv.rational());
    return;


}


//==========================================================================
void cart_to_bary(const SplineCurve& cv, const BaryCoordSystem3D& bc,
		  SplineCurve& cv_bc)
//==========================================================================
{
    ALWAYS_ERROR_IF(cv.dimension() != 3, "Dimension must be 3.");


    int n = cv.numCoefs();
    Vector3D cart;
    Vector4D bary;
    vector<double> new_coefs;
    if (!cv.rational()) {
	new_coefs.resize(4 * n);
	for (int i = 0; i < n; ++i) {
	    cart = Vector3D(cv.coefs_begin() + 3*i);
	    bary = bc.cartToBary(cart);
	    for (int j = 0; j < 4; ++j) {
		new_coefs[4*i + j] = bary[j];
	    }
	}
    } else {
	new_coefs.resize(5 * n);
	for (int i = 0; i < n; ++i) {
	    cart = Vector3D(cv.coefs_begin() + 3*i);
	    bary = bc.cartToBary(cart);
	    double w = cv.rcoefs_begin()[4*i + 3];
	    for (int j = 0; j < 4; ++j) {
		new_coefs[5*i + j] = bary[j] * w;
	    }
	    new_coefs[5*i + 4] = w;
	}
    }
    cv_bc = SplineCurve(n, cv.order(), cv.basis().begin(),
			new_coefs.begin(), 4, cv.rational());
    return;


}


//==========================================================================
void cart_to_bary(const SplineSurface& sf, const BaryCoordSystem3D& bc,
		  SplineSurface& sf_bc)
//==========================================================================
{
    ALWAYS_ERROR_IF(sf.dimension() != 3, "Dimension must be 3.");


    int nu = sf.numCoefs_u();
    int nv = sf.numCoefs_v();
    Vector3D cart;
    Vector4D bary;
    vector<double> new_coefs;
    if (!sf.rational()) {
	new_coefs.resize(4 * nu * nv);
	for (int iv = 0; iv < nv; ++iv) {
	    for (int iu = 0; iu < nu; ++iu) {
		int offset = nu * iv + iu;
		cart = Vector3D(sf.coefs_begin() + 3 * offset);
		bary = bc.cartToBary(cart);
		for (int j = 0; j < 4; ++j) {
		    new_coefs[4*offset + j] = bary[j];
		}
	    }
	}
    } else {
	new_coefs.resize(5 * nu * nv);
	for (int iv = 0; iv < nv; ++iv) {
	    for (int iu = 0; iu < nu; ++iu) {
		int offset = nu * iv + iu;
		cart = Vector3D(sf.coefs_begin() + 3 * offset);
		bary = bc.cartToBary(cart);
		double w = sf.rcoefs_begin()[4*offset + 3];
		for (int j = 0; j < 4; ++j) {
		    new_coefs[5*offset + j] = bary[j] * w;
		}
		new_coefs[5*offset + 4] = w;
	    }
	}
    }
    sf_bc = SplineSurface(nu, nv, sf.order_u(), sf.order_v(),
			  sf.basis_u().begin(), sf.basis_v().begin(),
			  new_coefs.begin(), 4, sf.rational());
    return;	
}


//==========================================================================
void cart_to_bary(const PointCloud3D& cloud, const BaryCoordSystem3D& bc,
		  PointCloud4D& cloud_bc)
//==========================================================================
{
    int num = cloud.numPoints();
    vector<Array<double, 4> > points_bc(num);
    Array<double, 3> cart;
    for (int i = 0; i < num; ++i) {
	cart = cloud.point(i);
	points_bc[i] = bc.cartToBary(cart);
    }

    cloud_bc = PointCloud4D(points_bc);

    return;	
}


//==========================================================================
void make_matrix(const SplineCurve& curve, int deg,
		 vector<vector<double> >& mat)
//==========================================================================
{
    // Create BernsteinPoly. In the rational case the weights are
    // included in an "extra" coordinate.
    int dim = curve.dimension();
    bool rational = curve.rational();
    vector<BernsteinPoly> beta;
    spline_to_bernstein(curve, beta);

    // Make vector of basis functions (with the curve plugged in) by
    // using recursion
    int num = (deg+1) * (deg+2) / 2;
    vector<BernsteinPoly> basis(num);
    vector<BernsteinPoly> tmp(num);
    basis[0] = BernsteinPoly(1.0);
    BernsteinPoly zero = BernsteinPoly(0.0);
    for (int r = 1; r <= deg; ++r) {
	int m = -1;
	int tmp_num = (r + 1) * (r + 2) / 2;
	fill(tmp.begin(), tmp.begin() + tmp_num, zero);
	for (int i = 0; i < r; ++i) {
	    for (int l = 0; l <= i; ++l) {
		++m;
		tmp[m] += beta[0] * basis[m];
		tmp[m + 1 + i] += beta[1] * basis[m];
		tmp[m + 2 + i] += beta[2] * basis[m];
	    }
	}
	basis.swap(tmp);
    }

    // Fill up the matrix mat
    int degt = curve.order() - 1;
    int numbas = deg * degt + 1;
    mat.resize(numbas);
    for (int row = 0; row < numbas; ++row) {
	mat[row].resize(num);
	for (int col = 0; col < num; ++col) {
	    mat[row][col] = basis[col][row];
	}
    }

    // If rational, include diagonal scaling matrix. Dividing the
    // D-matrix by the weights has the effect of multiplying the basis
    // with the same weights. (Included for numerical reasons only -
    // it makes the basis a partition of unity.)
    if (rational) {
        BernsteinPoly weights = BernsteinPoly(1.0);
	for (int i = 1; i <= deg; ++i)
	    weights *= beta[dim];
	for (int row = 0; row < numbas; ++row) {
	    double scaling = 1.0 / weights[row];
	    for (int col = 0; col < num; ++col) {
		mat[row][col] *= scaling;
	    }
	}
    }

//     // Check Frobenius norm
//     double norm = 0.0;
//     for (int irow = 0; irow < numbas; ++irow) {
//  	for (int icol = 0; icol < num; ++icol) {
//  	    norm += mat[irow][icol] * mat[irow][icol];
//  	}
//     }
//     norm = sqrt(norm);
//     cout << "Frobenius norm = " << norm << endl;

    return;
}


//==========================================================================
void make_matrix(const SplineSurface& surf, int deg,
		 vector<vector<double> >& mat)
//==========================================================================
{
    // Create BernsteinMulti. In the rational case the weights are
    // included in an "extra" coordinate.
    int dim = surf.dimension();
    bool rational = surf.rational();
    vector<BernsteinMulti> beta;
    spline_to_bernstein(surf, beta);

    // Make vector of basis functions (with the surface plugged in) by
    // using recursion
    int num = (deg+1) * (deg+2) * (deg+3) / 6;
    vector<BernsteinMulti> basis(num);
    vector<BernsteinMulti> tmp(num);
    basis[0] = BernsteinMulti(1.0);
    BernsteinMulti zero_multi = BernsteinMulti(0.0);
    for (int r = 1; r <= deg; ++r) {
	int m = -1;
	int tmp_num = (r + 1) * (r + 2) * (r + 3) / 6;
	fill(tmp.begin(), tmp.begin() + tmp_num, zero_multi);
	for (int i = 0; i < r; ++i) {
	    int k = (i + 1) * (i + 2) / 2;
	    for (int j = 0; j <= i; ++j) {
		for (int l = 0; l <= j; ++l) {
		    ++m;
		    tmp[m] += beta[0] * basis[m];
		    tmp[m + k] += beta[1] * basis[m];
		    tmp[m + 1 + j + k] += beta[2] * basis[m];
		    tmp[m + 2 + j + k] += beta[3] * basis[m];
		}
	    }
	}
	basis.swap(tmp);
    }

    // Fill up the matrix mat
    int deg_u = surf.order_u() - 1;
    int deg_v = surf.order_v() - 1;
    int numbas = (deg * deg_u + 1) * (deg * deg_v + 1);
    mat.resize(numbas);
    for (int row = 0; row < numbas; ++row) {
	mat[row].resize(num);
	for (int col = 0; col < num; ++col) {
	    mat[row][col] = basis[col][row];
	}
    }

    // If rational, include diagonal scaling matrix. Dividing the
    // D-matrix by the weights has the same effect as multiplying the
    // basis with the same weights. (Included for numerical reasons only -
    // it makes the basis a partition of unity.)
    if (rational) {
        BernsteinMulti weights = BernsteinMulti(1.0);
	for (int i = 1; i <= deg; ++i)
	    weights *= beta[dim];
	for (int row = 0; row < numbas; ++row) {
	    double scaling = 1.0 / weights[row];
	    for (int col = 0; col < num; ++col) {
		mat[row][col] *= scaling;
	    }
	}
    }

//     // Check Frobenius norm
//     double norm = 0.0;
//     for (int irow = 0; irow < numbas; ++irow) {
//  	for (int icol = 0; icol < num; ++icol) {
//  	    norm += mat[irow][icol] * mat[irow][icol];
//  	}
//     }
//     norm = sqrt(norm);
//     cout << "Frobenius norm = " << norm << endl;

    return;
}


//==========================================================================
void make_matrix(const PointCloud4D& cloud, int deg,
		 vector<vector<double> >& mat)
//==========================================================================
{
    // The matrix mat has the form mat_ij = B_{j,d}(p_i), where p_i is
    // the i'th point and B_{j,d} is the j'th triangluar Berstein
    // polynomial of degree d

    int numpts = cloud.numPoints();
    int numbas = (deg+1) * (deg+2) * (deg+3) / 6;
    mat.resize(numpts);

    // For each row - i.e. point - we make the Bernstein polynomials
    // by recursion. This we fill into mat.
    vector<double> basis(numbas);
    vector<double> tmp(numbas);
    Array<double, 4> pt;
    for (int i = 0; i < numpts; ++i) {
	pt = cloud.point(i);
	basis[0] = 1.0;
	for (int r = 1; r <= deg; ++r) {
	    int m = 0;
	    int tmp_num = (r + 1) * (r + 2) * (r + 3) / 6;
	    fill(tmp.begin(), tmp.begin() + tmp_num, 0.0);
	    for (int i = 0; i < r; ++i) {
		int k = (i + 1) * (i + 2) / 2;
		for (int j = 0; j <= i; ++j) {
		    for (int l = 0; l <= j; ++l) {
			tmp[m] += pt[0] * basis[m];
			tmp[m + k] += pt[1] * basis[m];
			tmp[m + 1 + j + k] += pt[2] * basis[m];
			tmp[m + 2 + j + k] += pt[3] * basis[m];
			++m;
		    }
		}
	    }
	    basis.swap(tmp);
	}
	mat[i].resize(numbas);
	for (int col = 0; col < numbas; ++col)
	    mat[i][col] = basis[col];
    }

    return;
}


//==========================================================================
void make_implicit_svd(vector<vector<double> >& mat,
		       vector<double>& b, double& sigma_min)
//==========================================================================
{
    int rows = (int)mat.size();
    int cols = (int)mat[0].size();
//     cout << "Rows = " << rows << endl
// 	 << "Cols = " << cols << endl;

    Matrix nmat;
    nmat.ReSize(rows, cols);
    for (int i = 0; i < rows; ++i) {
	for (int j = 0; j < cols; ++j) {
	    nmat.element(i, j) = mat[i][j];
	}
    }

    // Check if mat has enough rows. If not fill out with zeros.
    if (rows < cols) {
	RowVector zero(cols);
	zero = 0.0; // Initializes zero to a null-vector.
	for (int i = rows; i < cols; ++i) {
	    nmat &= zero; // & means horizontal concatenation in newmat
	}
    }

    // Perform SVD.
//     cout << "Running SVD..." << endl;
    static DiagonalMatrix diag;
    static Matrix V;
    Try {
	SVD(nmat, diag, nmat, V);
    } CatchAll {
	cout << Exception::what() << endl;
	b = vector<double>(cols, 0.0);
	sigma_min = -1.0;
	return;
    }

//     // Write out singular values.
//     cout << "Singular values:" << endl;
//     for (int ik = 0; ik < cols; ik++)
//  	cout << ik << "\t" << diag.element(ik, ik) << endl;

//     // Write out info about singular values
//     double s_min = diag.element(cols-1, cols-1);
//     double s_max = diag.element(0, 0);
//     cout << "Implicitization:" << endl
// 	 << "s_min = " << s_min << endl
// 	 << "s_max = " << s_max << endl
// 	 << "Ratio of s_min/s_max = " << s_min/s_max << endl;

//     // Find square sum of singular values
//     double sum = 0.0;
//     for (int i = 0; i < cols; i++)
//  	sum += diag.element(i, i) * diag.element(i, i);
//     sum = sqrt(sum);
//     cout << "Square sum = " << sum << endl;

    // Get the appropriate null-vector and corresponding singular value
    const double eps = 1.0e-15;
    double tol = cols * fabs(diag.element(0, 0)) * eps;
    int nullvec = 0;
    for (int i = 0; i < cols-1; ++i) {
	if (fabs(diag.element(i, i)) > tol) {
	    ++nullvec;
	}
    }
    sigma_min = diag.element(nullvec, nullvec);
//     cout << "Null-vector: " << nullvec << endl
// 	 << "sigma_min = " << sigma_min << endl;

    // Set the coefficients
    b.resize(cols);
    for (int jk = 0; jk < cols; ++jk)
	b[jk] = V.element(jk, nullvec);

    return;
}


//==========================================================================
void make_implicit_gauss(vector<vector<double> >& mat, vector<double>& b)
//==========================================================================
{
    int rows = (int)mat.size();
    int cols = (int)mat[0].size();
    int mindim = (rows < cols ? rows : cols);
//     cout << "Rows = " << rows << endl
// 	 << "Cols = " << cols << endl;

    // Gaussian elimination with complete pivoting. Algorithm 3.4.2 in
    // Golub and van Loan.
//     cout << "Gaussian elimination with complete pivoting..." << endl;
    vector<int> q(cols);
    for (int i = 0; i < cols; ++i)
	q[i] = i;
    for (int k = 0; k < mindim; ++k) {
	int pivi = k;
	int pivj = k;
	double big = 0.0;
	for (int i = k; i < rows; ++i) {
	    for (int j = k; j < cols; ++j) {
		if (fabs(mat[i][j]) > big) {
		    pivi = i;
		    pivj = j;
		    big = fabs(mat[i][j]);
		}
	    }
	}
	if (pivi != k) {
	    for (int l = 0; l < cols; ++l) {
		double dummy = mat[k][l];
		mat[k][l] = mat[pivi][l];
		mat[pivi][l] = dummy;
	    }
	}
	if (pivj != k) {
	    for (int l = 0; l < rows; ++l) {
		double dummy = mat[l][k];
		mat[l][k] = mat[l][pivj];
		mat[l][pivj] = dummy;
	    }
	    q[k] = pivj;
	}
	if (mat[k][k] != 0.0) {
	    for (int i = k + 1; i < rows; ++i) {
		double z = mat[i][k] / mat[k][k];
		for (int j = k + 1; j < cols; ++j) {
		    mat[i][j] -= z * mat[k][j];
		}
	    }
	}
    }

//     // Write out diagonal values.
//     for (int ik = 0; ik < mindim; ik++) {
//  	cerr << ik << "\t" << mat[ik][ik] << endl;
//     }

    // Find smallest pivots and prepare for backsubstitution
//     cout << "Find smallest pivot..." << endl;
//     const double eps = 1.0e-6;
    int pivmin = 0;
    double min = fabs(mat[0][0]);
    for (int i = 1; i < mindim; ++i) {
	// We test also for 0.0 because there may be artifacts in the
	// form of columns of zeros
// 	if (fabs(mat[i][i]) < min && fabs(mat[i][i]) > eps * min) {
	if (fabs(mat[i][i]) < min && fabs(mat[i][i]) > 0.0) {
	    pivmin = i;
	    min = fabs(mat[i][i]);
      }
    }
//     cout << "Smallest pivot: " << pivmin << " -> " << min << endl;
    b = vector<double>(cols, 0.0);
    b[pivmin] = 1.0;
//     for (int i = pivmin; i < cols; ++i)
// 	b[i] = 1.0;

    // Solve by backsubstitution
//     cout << "Backsubstitution..." << endl;
    for (int i = pivmin - 1; i >= 0; --i) {
	double sum = 0.0;
	for (int j = i + 1; j < cols; ++j)
	    sum -= mat[i][j] * b[j];
	b[i] = sum / mat[i][i];
    }

    // Unscramble when using complete pivoting
    for (int i = cols - 1; i >= 0; --i) {
	double dummy = b[q[i]];
	b[q[i]] = b[i];
	b[i] = dummy;
    }

//     // Write info about b
//     double norm = fabs(b[0]);
//     double bmax = fabs(b[0]);
//     double bmin = fabs(b[0]);
//     for (int i = 1; i < cols; ++i) {
// 	norm += fabs(b[i]);
// 	if (fabs(b[i]) < bmin)
// 	    bmin = fabs(b[i]);
// 	if (fabs(b[i]) > bmax)
// 	    bmax = fabs(b[i]);
//     }
//     norm /= cols;
//     cout << "Max-norm of b = " << norm << endl
// 	 << "bmin / bmax = " << bmin / bmax << endl;

    return;
}


//==========================================================================


} // namespace Go
