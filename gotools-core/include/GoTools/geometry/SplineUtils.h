#ifndef _SPLINEUTILS_H
#define _SPLINEUTILS_H


#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/geometry/SplineCurve.h"
#include "GoTools/utils/errormacros.h"
#include <math.h>
#include <algorithm>
#include <ctype.h>
#include "GoTools/utils/config.h"


namespace Go {

/// Utility functionality related to spline curves and surfaces or
/// functions operating on spline curves and surfaces
namespace SplineUtils {

  /// Utility functionality related to spline curves and surfaces or
  /// functions operating on spline curves and surfaces

    /// Transpose an (m x n) matrix of dim-dimensional points
    /// stored as an array in row-major order
    /// (i.e. all elements of a single row are stored together).
    /// \param dim the dimension of each of the points in the array
    /// \param m number of lines in the matrix
    /// \param n number of columns in the matrix 
    /// \param array_start pointer to the start of array where matrix
    ///                    is stored
    void GO_API transpose_array(int dim, int m, int n, double* array_start);

    /// Find the point in an array closest to the given base point
    /// (measured in the usual, Euclidean distance)    
    /// \param pt pointer to the base point (of dimension 'n')
    /// \param array pointer to the array in which we will search 
    /// \param n number of elements in the array
    /// \param dim dimension of the elements/base point
    int GO_API closest_in_array(const double* pt,
                                const double* array,
                                int n, 
                                int dim);

    /// Find the closest point on a triangle, in barycentric coordinates.
    /// \param pt the point for which we want to find the closest point on the triangle
    /// \param tri the triangle, specified by its three corners
    /// \param clo_dist2 the squared Euclidean distance from 'pt' to the closest
    ///                  point in the triangle.
    /// \return the barycentric coordinates of the closest point on the triangle
    Vector3D GO_API closest_on_triangle(const Vector3D& pt,
                                        const Vector3D tri[3],
                                        double& clo_dist2);

    /// Find the closest point on a line segment to a given point.
    /// \param pt the point for which we want to find the closest point on the 
    ///           line segment
    /// \param beg the start point of the line segment
    /// \param end the end point of the line segment
    /// \return the barycentric coordinate of the closest point on the line segment
    ///         (number between 0 and 1, where 0 is the start point and 1 is the 
    ///         end point)
    double GO_API closest_on_line_segment(const Vector3D& pt,
                                          const Vector3D& beg,
                                          const Vector3D& end);


    /// Find an approximate closest point on a rectangular grid of 3-dimensional points. 
    /// The method used is to 'triangulate' the grid, and locate the closest point on this
    /// triangulation.  We imagine that the point grid is uniformly parametrized from 0 to
    /// m-1 along each row and from 0 to n-1 along each column, and we refer to these 
    /// parameters as 'u' and 'v'.
    /// \param pt the point we want to search the closest point for
    /// \param array pointer to the array of 3D-points; coordinates should be stored 
    ///              consecutively in a row-wise fashion.
    /// \param m number of columns in the point array
    /// \param n number of rows in the point array
    /// \param clo_u upon function return; the u-parameter of the closest point, using the
    ///              parametrization explained above.
    /// \param clo_v upon function return; the v-parameter of the closest point, using the
    ///              parametrization explained above.
    void GO_API closest_on_rectgrid(const double* pt,
                                    const double* array,
                                    int m, int n,
                                    double& clo_u,
                                    double& clo_v);

    /// Find an approximate closest point on a sub-grid of a rectangular grid of 
    /// 3-dimensional points.  The method used is to 'triangulate' the grid, and locate
    /// the closest point on this triangulation.  We imagine that the parametrization of the
    /// sub-grid is from 'u_min' to 'u_max' along each row and from 'v_min' to 'v_max'
    /// along each column.
    /// \param pt the point we want to search the closest point for
    /// \param array pointer to the array of 3D-points; coordinates should be stored
    ///              consecutively in a row-wise fashion.
    /// \param u_min lowest column index of the subgrid we want to search
    /// \param u_max highest column index of the subgrid we want to search
    /// \param v_min lowest row index of the subgrid we want to search
    /// \param v_max highest row index of the subgrid we want to search
    /// \param nmb_coefs_u total number of columns in the complete grid (necessary to know
    ///                    in order to determine how far to jump to get to the next row
    ///                    in the subgrid.
    /// \param clo_u upon function return; the u-parameter of the closest point, using
    ///              the parametrization explained above.
    /// \param clo_v upon function return; the v-parameter of the closest point, using
    ///              the parametrization explained above.
    void GO_API closest_on_rectgrid(const double* pt, const double* array,
                                    int u_min, int u_max, int v_min, int v_max,
                                    int nmb_coefs_u,
                                    double& clo_u, double& clo_v);


    /// convert an array of rational coefficients to an array of nonrational coefficients
    /// \param rationals pointer to the array of rational coefficients
    /// \param coefs pointer to the array where the nonrational coefficients are to be
    ///              written
    /// \param num_coefs total number of coefficients
    /// \param dim dimension of coefficients (not counting the rational component).
    void GO_API make_coef_array_from_rational_coefs(const double* rationals,
                                                    double* coefs,
                                                    int num_coefs,
                                                    int dim);

    /// This function takes as input the position and a certain number of derivatives
    /// in homogenous space of a point on a rational curve.  It outputs the position 
    /// and derivatives of this point in non-homogenous ("ordinary") space.
    /// (Corresponds to s6ratder in SISL)
    /// \param eder pointer to the array where the point's position and derivatives in
    ///             homogenous coordinates (input) are consecutively stored
    /// \param idim the dimension of the non-homogenous space
    /// \param ider the number of derivatives sought (0 means that we only want to convert
    ///             the \em position to non-homogenous coordinates).
    /// \param gder pointer to the array where the result will be written (in same
    ///             order as 'eder').
    void GO_API curve_ratder(double const eder[],int idim,int ider,double gder[]);

    /// This function takes as input the position and a certain number of derivatives
    /// in homogenous space of a point on a rational surface.  It outputs the 
    /// positioin and derivatives of this point in non-homogenous ("ordinary") space.
    /// (Corresponds to s6strider in SISL)
    /// \param  eder pointer to double array of dimenson [(ider+1)*(ider+2)*(idim+1)/2]
    ///              containing the position and the derivative vectors
    ///              of the homogeneous surface at the point with parameter value
    ///              (epar[0],epar[1]).
    ///              (idim+1 is the number of components of each B-spline
    ///              coefficient, i.e. the dimension of the homogemous
    ///              space in which the surface lies.)
    ///              These vectors are stored in the following order:
    ///              First the idim+1 components of the position vector,
    ///              then the idim+1 components of the D(1,0) vector,
    ///              then the idim+1 components of the D(0,1) vector,
    ///              then the idim+1 components of the D(2,0) vector,
    ///              followed by D(1,1), D(0,2)
    ///              and so on up to the idim+1 components of the D(0,ider).
    /// \param idim The dimension of the non homogenous space
    /// \param ider The number of input derivatives with respect to both parameter directions.
    /// \param gder pointer to the array where the result will be written (in same order
    ///             as 'eder').
    void GO_API surface_ratder(double const eder[],int idim,int ider,double gder[]);

    /// Corresponds to s1701 in SISL
    /// This function computes in a compact format a line in the discrete
    /// B-spline matrix converting between an orginal basis
    /// "etau" and a new basis "et".
    /// \param ij The index of the new vertice
    /// \param imy An index on etau, where the input value are to be
    ///            etau(imy) <= et(ij) < etau(imy + 1).
    /// \param ik The order of the B-spline.
    /// \param in The number of the orginal vertices.
    /// \param et The new knot vector.
    /// \param etau The old knot vector.
    /// \param et An array ep(ik) to local use. Such that we
    ///           do not need to allocate the array locally after
    ///           each call.
    /// \param jpl The negativ difference between the index in galfa
    ///            and the real knot inserten matrix.
    /// \param jfi The index of the first element in the line j in the
    ///            the real knot inserten matrix whice is not zero.
    ///            The element with the index (jfi+jpl) in galfa
    ///            is the same as the element with index jfi in
    ///            the real line j in the knot inserten matrix.
    /// \param jla The index of the last element in the line j in the
    ///            real knot inserten matrix whice is not zero.
    ///            The element with the index (jla+jpl) in galfa
    ///            is the same as the element with index jla in
    ///            the real line j in the knot inserten matrix.
    /// \param galfa A compressed line in the knot inserten matrix.
    void GO_API osloalg(int ij,int imy,int ik,int in,int *jpl,int *jfi,int *jla,
                        const double *et, const double *etau,double *galfa);

    /// Corresponds to sh1922 in SISL
    /// Computes the B-spline refinement transformation matrix
    /// from the spline space generated by the knot vector etau
    /// to the refined spline space generated by the refined knot
    /// vector et.
    /// \param et Real array of length (im+ik) containing the refined
    ///           knot vector.
    /// \param im The dimension of the spline space corresponding to et.
    /// \param ik The order of the spline space.
    /// \param etau Real array of length (in+ik) containing the original 
    ///            knot vector.
    /// \param in The dimension of the spline space corresponding
    ///           to etau.
    /// \param ea Real array of dimension (im*ik) containing 
    ///           the B-spline refinement matrix from the knot vector
    ///           etau to the knot vector et. This matrix has
    ///           dimension im*in but since at most
    ///           ik entries are nonzero in each row, it can
    ///           be stored in a im*ik array together
    ///           with two integer arrays indicating the position
    ///           of the first and last nonzero elements in each
    ///           row.
    /// \param nfirst Integer array of dimension (im) containing 
    ///               pointers to the first nonzero element of each row 
    ///               of the B-spline refinement matrix from etau to et.
    /// \param nlast Integer array of dimension (im) containing 
    ///              pointers to the last nonzero element of each row 
    ///              of the B-spline refinement matrix from etau to et.
    void GO_API refmatrix(const double *et, int im, int ik, 
                          const double *etau, int in,
                          double *ea, int *nfirst,int *nlast);

    /// Assuming basis is cubic (i.e. order 4).
    /// Create the transformation matrix which extract the bezier coefs
    /// for the interval (knots[3], knots[4]).
    void GO_API splineToBezierTransfMat(const double* knots,
					std::vector<double>& transf_mat);

    /// Assuming surface is bi-cubic (i.e. order 4).
    /// Extract the bezier patch corr to the domain
    /// (knots_u[ind_u_min], knots_u[int_u_min+1])x(knots_v[ind_v_min], knots_v[int_v_min+1]).
    /// \param ind_u_min Basis pointer corresponding to umin.
    /// \param ind_v_min Basis pointer corresponding to vmin.
    void GO_API extractBezierCoefs(const double* coefs,
				   const int num_coefs_u, const int num_coefs_v,
				   const int ind_u_min, const int ind_v_min,
				   const std::vector<double>& transf_mat_u,
				   const std::vector<double>& transf_mat_v,
				   std::vector<double>& bezier_coefs);

    /// Method expecting bi-cubic input.
    void GO_API refinedBezierCoefsCubic(Go::SplineSurface& spline_sf,
					int ind_u_min, int ind_v_min,
					std::vector<double>& bez_coefs);

    // We insert knots so that all inner knots are of mult 'order'.
    shared_ptr<SplineSurface> GO_API refineToBezier(const Go::SplineSurface& spline_sf);

} // End of namespace SplineUtils

} // End of namespace Go


#endif

