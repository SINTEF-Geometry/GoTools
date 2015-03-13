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

#ifndef _LRSURFAPPROX_H_
#define _LRSURFAPPROX_H_

#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/lrsplines2D/LRSplineSurface.h"
#include "GoTools/lrsplines2D/LRSurfSmoothLS.h"
#include <vector>



//class LRSurfSmoothLS;

namespace Go
{
/// This class can generate a LR B-spline surface that approximates
/// a set of points for a given accuracy.

class LRSurfApprox
{
 public:

  /// Constructor given a parameterized point set
  /// \param points Parameterized point set given as (u1,v1,x1,y1,z1, u2, v2, ...)
  ///               The length of the array is (2+dim)x(the number of points)
  /// \param dim    The dimension of the geometry space
  /// \param epsge  Requested approximation accuracy
  /// \param closest_dist Check accuracy in closest point or in corresponding 
  ///                     parameter value
  /// \param repar Perform reparameterization during iterations
  LRSurfApprox(std::vector<double>& points, 
	       int dim, double epsge, bool init_mba=false, 
	       double mba_level = 0.0,
	       bool closest_dist=true, bool repar=false);

  /// Constructor given a parameterized point set and an initial spline surface
  /// \param srf    Given spline surface
  /// \param points Parameterized point set given as (u1,v1,x1,y1,z1, u2, v2, ...)
  ///               The length of the array is (2+dim)x(the number of points)
  ///               Note that the sequence of the points can be changed and this class
  ///               references the initial points (no copy)
  /// \param dim    The dimension of the geometry space
  /// \param epsge  Requested approximation accuracy
  /// \param closest_dist Check accuracy in closest point or in corresponding 
  ///                     parameter value
  /// \param repar Perform reparameterization during iterations
  LRSurfApprox(shared_ptr<SplineSurface>& srf,
	       std::vector<double>& points, 
	       double epsge, bool init_mba=false, 
	       double mba_level = 0.0, bool closest_dist=true,
	       bool repar=false);

  /// Constructor given a parameterized point set and an initial LR B-spline surface
  /// \param srf    Given LR B-spline surface
  /// \param points Parameterized point set given as (u1,v1,x1,y1,z1, u2, v2, ...)
  ///               The length of the array is (2+dim)x(the number of points)
  ///               Note that the sequence of the points can be changed and this class
  ///               references the initial points (no copy)
  /// \param dim    The dimension of the geometry space
  /// \param epsge  Requested approximation accuracy
  /// \param closest_dist Check accuracy in closest point or in corresponding 
  ///                     parameter value
  /// \param repar Perform reparameterization during iterations
  LRSurfApprox(shared_ptr<LRSplineSurface>& srf,
	       std::vector<double>& points, 
	       double epsge, bool init_mba=false, 
	       double mba_level = 0.0, bool closest_dist=true,
	       bool repar=false, bool check_init_accuracy=false);

  /// Constructor given a parameterized point set and the size of an initial
  /// spline space
  /// \param ncoef_u Number of coefficients in the 1. parameter direction
  /// \param order_u Order in the 1. parameter direction
  /// \param ncoef_v Number of coefficients in the 2. parameter direction
  /// \param order_v Order in the 2. parameter direction
  /// \param points Parameterized point set given as (u1,v1,x1,y1,z1, u2, v2, ...)
  ///               The length of the array is (2+dim)x(the number of points)
  ///               Note that the sequence of the points can be changed and this class
  ///               references the initial points (no copy)
  /// \param dim    The dimension of the geometry space
  /// \param epsge  Requested approximation accuracy
  /// \param closest_dist Check accuracy in closest point or in corresponding 
  ///                     parameter value
  /// \param repar Perform reparameterization during iterations
  LRSurfApprox(int ncoef_u, int order_u, int ncoef_v, int order_v,
	       std::vector<double>& points, 
	       int dim, double epsge, bool init_mba=false, 
	       double mba_level = 0.0,
	       bool closest_dist=true, bool repar=false);
  /// Constructor given a parameterized point set and an initial
  /// spline space
  /// \param order_u Order in the 1. parameter direction
  /// \param knots_u Knot vector in the 1. parameter direction
  /// \param order_v Order in the 2. parameter direction
  /// \param knots_v Knot vector in the 2. parameter direction
  /// \param points Parameterized point set given as (u1,v1,x1,y1,z1, u2, v2, ...)
  ///               The length of the array is (2+dim)x(the number of points)
  ///               Note that the sequence of the points can be changed and this class
  ///               references the initial points (no copy)
  /// \param dim    The dimension of the geometry space
  /// \param epsge  Requested approximation accuracy
  /// \param closest_dist Check accuracy in closest point or in corresponding 
  ///                     parameter value
  /// \param repar Perform reparameterization during iterations
  LRSurfApprox(int order_u, std::vector<double>& knots_u, 
	       int order_v, std::vector<double>& knots_v,
	       std::vector<double>& points, 
	       int dim, double epsge, bool init_mba=false, 
	       double mba_level = 0.0,
	       bool closest_dist=true, bool repar=false);

  /// Constructor given a parameterized point set and the size of an initial
  /// spline space
  /// \param ncoef_u Number of coefficients in the 1. parameter direction
  /// \param order_u Order in the 1. parameter direction
  /// \param ncoef_v Number of coefficients in the 2. parameter direction
  /// \param order_v Order in the 2. parameter direction
  /// \param points Parameterized point set given as (u1,v1,x1,y1,z1, u2, v2, ...)
  ///               The length of the array is (2+dim)x(the number of points)
  ///               Note that the sequence of the points can be changed and this class
  ///               references the initial points (no copy)
  /// \param dim    The dimension of the geometry space
  /// \param domain A possibly extende parameter domain (umin, umax, vmin, vmax)
  /// \param epsge  Requested approximation accuracy
  /// \param closest_dist Check accuracy in closest point or in corresponding 
  ///                     parameter value
  /// \param repar Perform reparameterization during iterations
  LRSurfApprox(int ncoef_u, int order_u, int ncoef_v, int order_v,
	       std::vector<double>& points, int dim, 
	       double domain[4], double epsge, bool init_mba=false, 
	       double mba_level = 0.0,
	       bool closest_dist=true, bool repar=false);

  /// Destructor
  ~LRSurfApprox();

    /// Sets the smoothing weight to something other than the default (1e-9).
    /// The value should lie in the unit interval, typically close to 0.
    /// \param smooth the new smoothing weight.
    void setSmoothingWeight(double smooth)
	{
	    ASSERT(smooth >= 0.0 && smooth <= 1.0);
	    smoothweight_ = smooth;
	}

    void setSmoothBoundary(bool smoothbd)
    {
      smoothbd_ = smoothbd;
    }

    /// Decide whether or not the total bondary of the surface should be kept fixed
    /// (i.e. unchanged by approximation process).   Default is true.  Cross derivatives
    /// will not be kept fixed.  (If you want to keep cross derivatives fixed, use
    /// the edgeFix() member function instead).
    /// \param fix_boundary if 'true' the boundary of the surface will not be modified,
    ///                     if 'false' it will be open to modification.
    void setFixBoundary(bool fix_boundary)
    {
      int fix = (fix_boundary) ? 1 : 0;
      edge_derivs_[0] = edge_derivs_[1] = edge_derivs_[2] = edge_derivs_[3] = fix;
      setCoefKnown();
    }

    void setFixCorner(bool fix_corner)
    {
      fix_corner_ = fix_corner;
    }

    /// Decide whether specific edges of the surface's boundary should be kept fixed
    /// (i.e. unchanged by approximation process), as well as a certain number of cross-
    /// derivatives across these curves.
    /// \param edge_fix pointer to an array of four integers specifying to which extent
    ///                 each surface edge should be kept fixed during the approximation 
    ///                 process.
    ///                 A value of 0 means that it will not be kept fixed, 1 means that
    ///                 its position will be kept fixed, derivatives at the boundary
    ///                 are currently NOT kept fixed
    ///                 The integers are associated with the surface edges starting 
    ///                 with the edge 'v=vmin' and moving counterclockwise.
    void edgeFix(int edge_fix[])  // CCV
    {
      for (int ki=0; ki<4; ki++)
	edge_derivs_[ki] = std::min(1, edge_fix[ki]);
      setCoefKnown();
    }

    void setTurn3D(int iter)
    {
      to3D_ = iter;
    }

    void setGridInfo(double grid_start[], double cell_size[])
    {
      grid_ = true;
      grid_start_[0] = grid_start[0];
      grid_start_[1] = grid_start[1];
      cell_size_[0] = cell_size[0];
      cell_size_[1] = cell_size[1];
    }

    /// Decide if ghost points should be constructed
    void setMakeGhostPoints(bool make_ghost_points)
    {
      make_ghost_points_ = make_ghost_points;
    }

    /// Decide if the MBA algorithm should be used to approximate
    /// points (default == false)
    void setUseMBA(bool useMBA)
    {
      useMBA_ = useMBA;
    }

    /// Add lower constraint. Only functional (1D surface)
    void addLowerConstraint(double minval)
    {
      has_min_constraint_ = true;
      minval_ = minval;
    }

    /// Add upper constraint. Only functional (1D surface)
    void addUpperConstraint(double maxval)
    {
      has_max_constraint_ = true;
      maxval_ = maxval;
    }

    /// Set local constraint on coefficient. Only functional (1D surface)
    void setLocalConstraint(double constraint_factor)
    {
      has_local_constraint_ = true;
      constraint_fac_ = constraint_factor;
    }

    /// Set minimum element size
    void setMinimumElementSize(double usize_min, double vsize_min)
    {
      usize_min_ = usize_min;
      vsize_min_ = vsize_min;
    }

    /// Iteration count on which to swith to MBA (if LS)
    void setSwitchToMBA(int iter)
    {
      toMBA_ = iter;
    }

    /// Set if an initial MBA surface should be used in LS approximation
    /// Default is false
    void setInitMBA(bool initMBA, double start_coef=0.0)
    {
      initMBA_ = initMBA;
      initMBA_coef_ = start_coef;
    }

    /// Whether or not intermediate information should be written to
    /// standard output (default is not)
    void setVerbose(bool verbose)
    {
      verbose_ = verbose;
    }

    /// When everything else is set, this function can be used to run the 
    /// approximation process and fetch the approximated surface.
    /// \retval maxdist report the maximum distance between the approximated 
    ///                 surface and the data points
    /// \retval avdist report the average distance between the approximated 
    ///                surface and the datapoints
    /// \param nmb_out_eps report the number of data points that were found to
    ///                    lie outside the geometric tolerance (as specified by 
    ///                    the 'aepsge' argument to the ApproxSurf constructor.
    /// \param max_iter maximum number of iterations
    /// \return a shared pointer to the generated SplineCurve, approximating 
    ///         the points as specified.
    shared_ptr<LRSplineSurface> getApproxSurf(double& maxdist, 
					      double& avdist_all,
					      double& avdist,
					      int& nmb_out_eps, 
					      int max_iter=4);

 private:
    shared_ptr<LRSplineSurface> srf_;
    int nmb_pts_;
    std::vector<double>& points_;  // Reference to input points and parameter values
    std::vector<int> coef_known_;
    shared_ptr<LRSplineSurface> prev_;  // Previous surface, no point information
    // in elements
    
    bool useMBA_;    // Only LR-MBA
    int toMBA_;      // Start with LR-MBA at the given iteration step
    bool initMBA_;   // The initial surface is made using LR-MBA
    double initMBA_coef_;  // Initial hight of constant surface

    int edge_derivs_[4];
    double maxdist_;
    double maxdist_prev_;
    double avdist_;
    double avdist_all_;
    double avdist_all_prev_;
    int outsideeps_;
    double aepsge_;
    double smoothweight_;
    bool smoothbd_;
    bool repar_;
    bool check_close_;
    bool fix_corner_;
    int to3D_;
    bool grid_;
    bool check_init_accuracy_;
    double grid_start_[2];
    double cell_size_[2];
    bool initial_surface_;

    bool has_min_constraint_;
    bool has_max_constraint_;
    double minval_;
    double maxval_;
    bool has_local_constraint_;
    double constraint_fac_;
    bool verbose_;
    double usize_min_;  // Minimum element size in u direction, negative 
    // if not set
    double vsize_min_;  // Minimum element size in v direction, negative 
    // if not set

    bool fix_boundary_;
    bool make_ghost_points_;

    /// Define free and fixed coefficients
    void setCoefKnown();
    void unsetCoefKnown();
    void updateCoefKnown();
    void setCoefKnown(Direction2D d, Direction2D d2, bool atstart, int fixcoef);

    /// Perform least squares approximation with a smoothing term
    void performSmooth(LRSurfSmoothLS *LSapprox);

    void computeAccuracy(std::vector<Element2D*>& ghost_elems);
    // The same as the above, but with OpenMP support (if flag is turned on).
    void computeAccuracy_omp(std::vector<Element2D*>& ghost_elems);
    void computeAccuracyElement(std::vector<double>& points, int nmb, int del,
				RectDomain& rd, const Element2D* elem);
    // The same as the above, but with OpenMP support (if flag is turned on).
    void computeAccuracyElement_omp(std::vector<double>& points, int nmb, int del,
				    RectDomain& rd, const Element2D* elem);
    /// Refine surface
    int refineSurf();
    void refineSurf2();

    /// Create initial LR B-spline surface
    void makeInitSurf(int dim);
    void makeInitSurf(shared_ptr<SplineSurface> surf);
    void makeInitSurf(int dim, int ncoef_u, int order_u, int ncoef_v, int order_v,
		      double domain[4]);
    void makeInitSurf(int dim, int ncoef_u, int order_u, int ncoef_v, 
		      int order_v);
    void makeInitSurf(int dim, int ncoef_u, int order_u, int ncoef_v, 
		      int order_v, double *knots_u, double *knots_v);

    /// Create spline surface
    shared_ptr<SplineSurface> createSurf(double* points, int nmb_pts,
					 int dim, int ncoef_u, int order_u, 
					 int ncoef_v, int order_v, 
					 double *knots_u, double *knots_v,
					 double smoothweight,
					 double& maxdist, double& avdist,
					 int& nmb_outside);

    /// Parameter domain surrounding the parameter values of all data points
    void computeParDomain(int dim, double& umin, double& umax, double& vmin, double& vmax);

    void defineRefs(LRBSpline2D* bspline,
		    std::vector<LRSplineSurface::Refinement2D>& refs,
		    int choice);

    void checkFeasibleRef(Element2D* elem, 
			  std::vector<LRSplineSurface::Refinement2D>& refs,
			  std::vector<Element2D*>& affected);

    void constructGhostPoints(std::vector<double>& ghost_points);
    void constructLocalGhostPts(double *startpt, int kn2,
				int dim, double u1, double u2,
				double v1, double v2, 
				int nmb_u, int nmb_v,
				std::vector<double>& ghostpts);

    void constructInnerGhostPoints();

    void updateGhostElems(std::vector<Element2D*>& elems, bool enable_omp = false);
    void updateGhostPoints(std::vector<Element2D*>& elems);

    void addConstraintGhostPoints();

    void adaptSurfaceToConstraints();

    // Turn function into a 3D surface
    void turnTo3D();
};

}  // namespace Go

#endif
