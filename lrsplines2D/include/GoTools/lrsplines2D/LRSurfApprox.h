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
#include "GoTools/geometry/RectDomain.h"
#include "GoTools/creators/Eval1D3DSurf.h"
#include "GoTools/lrsplines2D/LRSplineSurface.h"
#include "GoTools/lrsplines2D/LRSurfSmoothLS.h"
#include <vector>



//class LRSurfSmoothLS;

namespace Go
{
/// This class can generate a LR B-spline surface that approximates
/// a set of points with a given tolerance. The resulting accuracy depends
  /// also on the characteristics of the point cloud and other parameters.

class LRSurfApprox
{
 public:
  /// Storage for variable tolerances depending on domain
  struct TolBox {
    RectDomain box;
    double tol;

    void setVal(double umin, double umax, double vmin, double vmax,
		double tolerance)
    {
      box = RectDomain(Vector2D(umin, vmin), Vector2D(umax, vmax));
      tol = tolerance;
    }

    void setTol(double tolerance)
    {
      tol = tolerance;
    }

    void translateBox(double udel, double vdel)
    {
      box.move(Vector2D(udel,vdel));
    }

    bool contains(double uval, double vval)
    {
      return box.isInDomain(Vector2D(uval,vval), 0.0);
    }
  };
  
  /// Constructor given a parameterized point set
  /// \param points Parameterized point set given as (u1,v1,x1,y1,z1, u2, v2, ...)
  ///               The length of the array is (2+dim)x(the number of points)
  /// \param dim    The dimension of the geometry space
  /// \param epsge  Requested approximation accuracy
  /// \param init_mba If true, initiate surface with mba_level
  /// \param mba_level Value of all coefficients in initial surface if mba_level is set
  /// \param closest_dist Check accuracy in closest point or in corresponding 
  ///                     parameter value
  /// \param repar Perform reparameterization during iterations
  LRSurfApprox(std::vector<double>& points, 
	       int dim, double epsge, bool init_mba=false, 
	       double mba_level = 0.0,
	       bool closest_dist=true, bool repar=false);

  /// Constructor given an initial spline surface
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
	       double epsge, bool closest_dist=true,
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
	       double epsge, bool closest_dist=true,
	       bool repar=false, bool check_init_accuracy=false);

  /// Constructor given an initial LR B-spline surface where the points already 
  /// are distributed to the elements. 
  /// \param srf    Given surface
  /// \param epsge  Requested approximation accuracy
  /// \param init_mba If true, initiate surface with mba_level
  /// \param mba_level Value of all coefficients in initial surface if mba_level is set
  /// \param closest_dist Check accuracy in closest point or in corresponding 
  ///                     parameter value
  /// \param repar Perform reparameterization during iterations
  LRSurfApprox(shared_ptr<LRSplineSurface>& srf,
	       std::vector<double>& points, 
	       double epsge, bool init_mba=true, double mba_level = 0.0,
	       bool repar=false, bool closest_dist=true);

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

  /// Add significant points to whom special attention must be paid
  void addSignificantPoints(std::vector<double>& sign_points, 
			    double sign_epsge)
  {
    sign_points_.insert(sign_points_.end(), sign_points.begin(),
			sign_points.end());
    nmb_sign_ = (int)sign_points_.size()/(2+srf_->dimension());
    sign_aepsge_ = sign_epsge;
  }

  /// Set factor for approximation of significant points (default == 5)
  void setSignificantFactor(double significant_fac)
  {
    significant_fac_ = significant_fac;
  }

  /// Get results of significant point approximation
  void getSignificantPointInfo(double& maxdist_sign, double& avdist_sign,
			       int& outside_sign)
  {
    maxdist_sign  = maxdist_sign_;
    avdist_sign = avdist_sign_;
    outside_sign = outsideeps_sign_;
  }

    /// Sets the smoothing weight to something other than the default (1e-3).
    /// The value should lie in the unit interval, typically close to 0.
    /// \param smooth the new smoothing weight.
    void setSmoothingWeight(double smooth)
	{
	    ASSERT(smooth >= 0.0 && smooth <= 1.0);
	    smoothweight_ = smooth;
	}

  /// Additional smoothing at boundaries. Default is false.
    void setSmoothBoundary(bool smoothbd)
    {
      smoothbd_ = smoothbd;
    }

    /// Decide whether or not the total bondary of the surface should be kept fixed
    /// (i.e. unchanged by approximation process).   Default is false.  Cross derivatives
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

  /// Decide of the surface corners should be kept fixed during the approximation
  /// process. Default is false.
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
  /// Default is 0 for all edges
    void edgeFix(int edge_fix[])  // CCV
    {
      for (int ki=0; ki<4; ki++)
	edge_derivs_[ki] = std::min(1, edge_fix[ki]);
      setCoefKnown();
    }

  /// Change 1D surface to 3D by letting the parameterization define the two
  /// first coordinates. Default is not (iter = -1)
    void setTurn3D(int iter)
    {
      to3D_ = iter;
    }

  /// Special treatment of gridded point cloud. Default is not. Does not give
  /// added accuracy
    void setGridInfo(double grid_start[], double cell_size[])
    {
      grid_ = true;
      grid_start_[0] = grid_start[0];
      grid_start_[1] = grid_start[1];
      cell_size_[0] = cell_size[0];
      cell_size_[1] = cell_size[1];
    }

    /// Decide if ghost points should be constructed. Default is false
    void setMakeGhostPoints(bool make_ghost_points)
    {
      make_ghost_points_ = make_ghost_points;
    }

    /// Decide if the MBA algorithm should be used to approximate
    /// points during the entire computation (default == false)
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

    /// Set number of iterations for each knot insertion using mba
    /// Default = 2
    void setMBAiter(int nmb_iter)
    {
      nmb_mba_iter_ = nmb_iter;
    }

    /// Set if only points with approximation errors of a given sign should
    /// be included in the computations (only 1D surface, only MBA)
    void setMBASign(int sgn)
    {
      mba_sgn_ = sgn;
    }

    /// Add information about initial knot values used in refinement
    void setInitialKnots(std::vector<double>& init_knots_u,
			 std::vector<double>& init_knots_v)
    {
      init_knots_u_ = init_knots_u;
      init_knots_v_ = init_knots_v;
    }

    /// Set flag on whether or not a search for outliers should be performed
    void setOutlierFlag(bool outlier_detection)
    {
      outlier_detection_ = outlier_detection;
    }

    /// Fetch number of registered outliers (may be non-zero only if
    /// outlier flag is true
    int getNmbOutliers()
    {
      return nmb_outliers_;
    }

    /// Get outlier points (only after call to getApproxSurf)
    void getOutlierPts(std::vector<double>& outliers, int& nmb_outliers);

    /// Get non-outlier points (only after call to getApproxSurf)
    void getRegularPts(std::vector<double>& regular, int& nmb_regular);

    /// Get classified points (only after call to getApproxSurf)
    void getClassifiedPts(std::vector<double>& outliers, int& nmb_outliers,
			  std::vector<double>& regular, int& nmb_regular);

    /// Set information about variable tolerance depending on elevation
  /// (default not triggered)
    void setVarTol(double fac_pos, double fac_neg, bool var_tol_sign = false)
    {
      has_var_tol_ = true;
      var_fac_pos_ = fac_pos;
      var_fac_neg_ = fac_neg;
      has_var_tol_sign_ = var_tol_sign;
    }

  /// Set minimum tolerance. Used if a variable tolerance depending on elevation is applied
    void setMinTol(double mintol)
    {
      mintol_ = mintol;
    }

  /// Unset variable tolerance depending on elevation
    void unsetVarTol()
    {
      has_var_tol_ = false;
      has_var_tol_sign_ = false;
      var_fac_pos_ = var_fac_neg_ = 1.0;
    }

  // Set variable tolerance depending on domain
    void setVarTolBox(std::vector<TolBox> tolerances)
    {
      tolerances_ = tolerances;
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

    /// Fetch additional accuracy information relevant for variable
    /// tolerance. To be called after getApproxSurf
    void fetchOutsideTolInfo(double& max_outside,
			     double& average_outside)
    {
      max_outside = maxout_;
      average_outside = avout_;
    }

    /// Refinement strategy
    /// category1: 1 = full span, 2 = minimum span l, 3 = minimum span u,
    ///            4 = minimum span c, 5 = structured mesh, 6 = restricted mesh
    ///            7 = restricted mesh with element extension
    /// alter: 0 = refine in both parameter directions, 1 = refine in alternating directions
    /// threshold1: 0 = no, 1 = td, 2 = tn, 3 = tk, 4 = td+tk (not all combinations with
    ///            category exist, if non-existing threshold is set to 0
    /// swap: Change category to category2 when approximation efficiency is below this level
    /// category2, threshold2: after swap, as above
    void setRefinementStrategy(int category1, int alter, int threshold1,
			       double swap, int category2, int threshold2)
    {
      category1_ = category1;
      category2_ = category2;
      alter_ = alter;
      threshold1_ = threshold1;
      threshold2_ = threshold2;
      swap_ = swap;
    }

  /// Return informaion on applied refinement strategy
    void getRefinementStrategy(int& category1, int& alter, int& threshold1,
			       double& swap, int& category2, int& threshold2)
    {
      category1 = category1_;
      category2 = category2_;
      threshold1 = threshold1_;
      threshold2 = threshold2_;
      alter = alter_;
      swap = swap_;
    }

     /// Compute feature output in ncell x ncell grid
  /// The result is stored in files called "cellinfox.txt" where x is replaced
  /// by the iteration level, see \link LRFeatureUtils.h \endlink for implented
  /// features
    void setFeatureOut(int ncell)
    {
      write_feature_ = true;
      ncell_ = ncell;
    }

  /// Unset request to compute feature output
    void unsetFeatureOut()
    {
      write_feature_ = false;
    }

  void setAIC(bool AIC)
  {
    compute_AIC_ = AIC;
  }

  void getAICInfo(std::vector<double>& AIC_info, std::vector<int>& ncoef)
  {
    AIC_info.clear();
    ncoef.clear();
    AIC_info.insert(AIC_info.end(), AIC_.begin(), AIC_.end());
    ncoef.insert(ncoef.end(), ncoef_.begin(), ncoef_.end());
  }

private:
    shared_ptr<LRSplineSurface> srf_;
    shared_ptr<Eval1D3DSurf> evalsrf_;
    int nmb_pts_;
    int nmb_sign_;
    int nmb_outliers_;
    std::vector<double>& points_;  // Reference to input points and parameter values
    std::vector<double> sign_points_;

    std::vector<int> coef_known_;
    shared_ptr<LRSplineSurface> prev_;  // Previous surface, no point information
    // in elements

    bool useMBA_;    // Only LR-MBA
    int toMBA_;      // Start with LR-MBA at the given iteration step
    bool initMBA_;   // The initial surface is made using LR-MBA
    double initMBA_coef_;  // Initial hight of constant surface

    std::vector<double> init_knots_u_; // Initial knots to select for refinement
    std::vector<double> init_knots_v_; // Initial knots to select for refinement

    int edge_derivs_[4];
    double maxdist_;
    double maxdist_prev_;
    double maxdist_sign_;
    double avdist_;
    double avdist_all_;
    double avdist_all_prev_;
    double avdist_sign_;
    int outsideeps_;
    int outsideeps_sign_;
    double maxout_;
    double avout_;
    double aepsge_;
    double sign_aepsge_;
    double smoothweight_;
    double significant_fac_;
    int maxLScoef_;
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
    int nmb_mba_iter_;
    int mba_sgn_;
    bool verbose_;
    double usize_min_;  // Minimum element size in u direction, negative 
    // if not set
    double vsize_min_;  // Minimum element size in v direction, negative 
    // if not set

    double prev_el_out_;
    double prev_thresh_;
    
    bool fix_boundary_;
    bool make_ghost_points_;
    bool outlier_detection_;

    // Variable tolerance (linear with distane from zero)
    bool has_var_tol_;
    double var_fac_pos_;
    double var_fac_neg_;
    double mintol_;
    bool has_var_tol_sign_;
    std::vector<TolBox> tolerances_;

    // Refinement strategy
    int category1_, category2_, alter_, threshold1_, threshold2_;
    double swap_;
    
    // Features output
    bool write_feature_;
    int ncell_;

  // AIC
  bool compute_AIC_;
  std::vector<double> AIC_;
  std::vector<int> ncoef_;

    void initDefaultParams();

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
				RectDomain& rd, const Element2D* elem,
				std::vector<double>& prev_points_dist);
    // The same as the above, but with OpenMP support (if flag is turned on).
    void computeAccuracyElement_omp(std::vector<double>& points, int nmb, int del,
				    RectDomain& rd, const Element2D* elem,
				    std::vector<double>& prev_points_dist);

    void runMBAUpdate(bool computed_accuracy);

    int defineOutlierPts(Element2D* element, 
			 std::vector<double>& prev_dist, double lim,
			 double rad);

    //double density);
    /// Refine surface
    int refineSurf(int iter, int& dir, double threshold);
    int refineSurf3(int iter, int& dir, double threshold);
    int refineSurf4(int& dir, double threshold);
    void getRefineExtension(Element2D *elem, Direction2D fixdir,
			    int strategy, double& ppar, double& pmin, double& pmax,
			    std::set<std::pair<Element2D*,std::pair<Direction2D,double> > >& unbalanced_elem);

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
					 int& nmb_outside,
					 double *points2=0, int nmb_pts2=0);

    /// Parameter domain surrounding the parameter values of all data points
    void computeParDomain(int dim, double& umin, double& umax, double& vmin, double& vmax);

    void defineRefs(LRBSpline2D* bspline, double average_out, int dir,
		    std::vector<LRSplineSurface::Refinement2D>& refs_x,
		    std::vector<LRSplineSurface::Refinement2D>& refs_y,
		    int choice,
		    std::vector<Element2D*>& elem_div);

    void checkFeasibleRef(Element2D* elem, int dir, int iter,
			  std::vector<LRSplineSurface::Refinement2D>& refs_x,
			  std::vector<LRSplineSurface::Refinement2D>& refs_y,
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

    void getCandElements(double x, double y, double rad, 
			 Element2D* start_elem,
			 std::vector<Element2D*>& elems);

    void appendRef(std::vector<LRSplineSurface::Refinement2D>& refs,
		   LRSplineSurface::Refinement2D& curr_ref, double tol);

    // Turn function into a 3D surface
    void turnTo3D();
};

}  // namespace Go

#endif
