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

#ifndef _LRVOLAPPROX_H_
#define _LRVOLAPPROX_H_

#include "GoTools/trivariate/SplineVolume.h"
#include "GoTools/trivariate/RectTriDomain.h"
#include "GoTools/lrsplines3D/LRSplineVolume.h"
#include <vector>


namespace Go
{
/// This class can generate a LR B-spline volume that approximates
/// a set of points for a given accuracy.

class LRVolApprox
{
 public:

  /// Storage of variable tolerances
  struct TolBox {
    RectTriDomain box;
    double tol;

    void setVal(double umin, double umax, double vmin, double vmax,
		 double wmin, double wmax,double tolerance)
    {
      box = RectTriDomain(Vector3D(umin, vmin, wmin), Vector3D(umax, vmax, wmax));
      tol = tolerance;
    }

    void setTol(double tolerance)
    {
      tol = tolerance;
    }

    void translateBox(double udel, double vdel, double wdel)
    {
      box.move(Vector3D(udel,vdel,wdel));
    }

    bool contains(double uval, double vval, double wval)
    {
      return box.isInDomain(Vector3D(uval,vval,wval), 0.0);
    }
  };
  
  /// Constructor given a parameterized point set
  /// \param points Parameterized point set given as (u1,v1,w1,x1,y1,z1, u2,v2,w2,...)
  ///               The length of the array is (3+dim)x(the number of points)
  /// \param dim    The dimension of the geometry space
  /// \param epsge  Requested approximation accuracy
  /// \param closest_dist Check accuracy in closest point or in corresponding 
  ///                     parameter value
  /// \param repar Perform reparameterization during iterations
  LRVolApprox(std::vector<double>& points, int dim, double epsge, 
              double mba_level = 0.0,
              bool closest_dist=true, 
              bool repar=false);

  /// Constructor given a parameterized point set and an initial spline volume
  /// \param vol    Given spline volume
  /// \param points Parameterized point set given as (u1,v1,w1,x1,y1,z1, u2,v2,w2,...)
  ///               The length of the array is (3+dim)x(the number of points)
  ///               Note that the sequence of the points can be changed and this class
  ///               references the initial points (no copy)
  /// \param dim    The dimension of the geometry space
  /// \param epsge  Requested approximation accuracy
  /// \param closest_dist Check accuracy in closest point or in corresponding 
  ///                     parameter value
  /// \param repar Perform reparameterization during iterations
  LRVolApprox(shared_ptr<SplineVolume>& vol, std::vector<double>& points, 
              double epsge, 
              bool closest_dist=true,
              bool repar=false);

  /// Constructor given a parameterized point set and an initial LR B-spline volume
  /// \param vol    Given LR B-spline volume
  /// \param points Parameterized point set given as (u1,v1,w1,x1,y1,z1, u2,v2,w2,...)
  ///               The length of the array is (3+dim)x(the number of points)
  ///               Note that the sequence of the points can be changed and this class
  ///               references the initial points (no copy)
  /// \param dim    The dimension of the geometry space
  /// \param epsge  Requested approximation accuracy
  /// \param closest_dist Check accuracy in closest point or in corresponding 
  ///                     parameter value
  /// \param repar Perform reparameterization during iterations
  LRVolApprox(shared_ptr<LRSplineVolume>& vol, std::vector<double>& points, 
	       double epsge, 
               bool closest_dist=true,
	       bool repar=false, 
               bool check_init_accuracy=false);

  /// Constructor given a parameterized point set and the size of an initial
  /// spline space
  /// \param ncoef_u Number of coefficients in the 1. parameter direction
  /// \param order_u Order in the 1. parameter direction
  /// \param ncoef_v Number of coefficients in the 2. parameter direction
  /// \param order_v Order in the 2. parameter direction
  /// \param ncoef_w Number of coefficients in the 3. parameter direction
  /// \param order_w Order in the 3. parameter direction
  /// \param points Parameterized point set given as (u1,v1,w1,x1,y1,z1, u2,v2,w2,...)
  ///               The length of the array is (3+dim)x(the number of points)
  ///               Note that the sequence of the points can be changed and this class
  ///               references the initial points (no copy)
  /// \param dim    The dimension of the geometry space
  /// \param epsge  Requested approximation accuracy
  /// \param closest_dist Check accuracy in closest point or in corresponding 
  ///                     parameter value
  /// \param repar Perform reparameterization during iterations
  LRVolApprox(int ncoef_u, int order_u, int ncoef_v, int order_v, 
	      int ncoef_w, int order_w,
              std::vector<double>& points, 
              int dim, double domain[], double epsge, 
              //bool init_mba=false, 
              double mba_level = 0.0,
              bool closest_dist=true, bool repar=false);

  /// Constructor given a parameterized point set and an initial
  /// spline space
  /// \param order_u Order in the 1. parameter direction
  /// \param knots_u Knot vector in the 1. parameter direction
  /// \param order_v Order in the 2. parameter direction
  /// \param knots_v Knot vector in the 2. parameter direction
  /// \param order_w Order in the 3. parameter direction
  /// \param knots_w Knot vector in the 3. parameter direction
  /// \param points Parameterized point set given as (u1,v1,w1,x1,y1,z1, u2,v2,w2,...)
  ///               The length of the array is (3+dim)x(the number of points)
  ///               Note that the sequence of the points can be changed and this class
  ///               references the initial points (no copy)
  /// \param dim    The dimension of the geometry space
  /// \param epsge  Requested approximation accuracy
  /// \param closest_dist Check accuracy in closest point or in corresponding 
  ///                     parameter value
  /// \param repar Perform reparameterization during iterations
  LRVolApprox(int order_u, std::vector<double>& knots_u, 
              int order_v, std::vector<double>& knots_v,
              int order_w, std::vector<double>& knots_w,
              std::vector<double>& points, 
              int dim, double epsge, 
              //bool init_mba=false, 
              double mba_level = 0.0,
              bool closest_dist=true, bool repar=false);
  
  /// Constructor given a parameterized point set and the size of an initial
  /// spline space
  /// \param ncoef_u Number of coefficients in the 1. parameter direction
  /// \param order_u Order in the 1. parameter direction
  /// \param ncoef_v Number of coefficients in the 2. parameter direction
  /// \param order_v Order in the 2. parameter direction
  /// \param ncoef_w Number of coefficients in the 3. parameter direction
  /// \param order_w Order in the 3. parameter direction
  /// \param points Parameterized point set given as (u1,v1,w1,x1,y1,z1, u2,v2,w2,...)
  ///               The length of the array is (3+dim)x(the number of points)
  ///               Note that the sequence of the points can be changed and this class
  ///               references the initial points (no copy)
  /// \param dim    The dimension of the geometry space
  /// \param domain A possibly extended parameter domain (umin, umax, vmin, vmax, wmin, wmax)
  /// \param epsge  Requested approximation accuracy
  /// \param closest_dist Check accuracy in closest point or in corresponding 
  ///                     parameter value
  /// \param repar Perform reparameterization during iterations
  LRVolApprox(int ncoef_u, int order_u, int ncoef_v, int order_v, int ncoef_w, int order_w,
              std::vector<double>& points, int dim, 
              double epsge, 
              //bool init_mba=false, 
              double mba_level = 0.0,
              bool closest_dist=true, bool repar=false);

  /// Destructor
  ~LRVolApprox();
  
  /// Decide whether or not the total boundary of the volume should be kept fixed
  /// (i.e. unchanged by approximation process). Default is true. Cross derivatives
  /// will not be kept fixed. (If you want to keep cross derivatives fixed, use
  /// the faceFix() member function instead).
  /// \param fix_boundary if 'true' the boundary of the volume will not be modified,
  ///                     if 'false' it will be open to modification.
  void setFixBoundary(bool fix_boundary)
  {
    int fix = (fix_boundary) ? 1 : 0;
    face_derivs_[0] = face_derivs_[1] 
                    = face_derivs_[2]
                    = face_derivs_[3]
                    = face_derivs_[4]
                    = face_derivs_[5] = fix;
    setCoefKnown();
  }
  
  void setFixCorner(bool fix_corner)
  {
    fix_corner_ = fix_corner;
  }

  /// Decide whether specific faces of the volume's boundary should be kept fixed
  /// (i.e. unchanged by approximation process), as well as a certain number of cross-
  /// derivatives across these volumes. (does this make sense??)
  /// \param face_fix pointer to an array of four integers specifying to which extent
  ///                 each volume face should be kept fixed during the approximation 
  ///                 process.
  ///                 A value of 0 means that it will not be kept fixed, 1 means that
  ///                 its position will be kept fixed, derivatives at the boundary
  ///                 are currently NOT kept fixed
  ///                 The integers are associated with the volume faces starting 
  ///                 with the face 'w=wmin' and moving counterclockwise.
    void faceFix(int face_fix[])  // CCV
    {
      for (int ki=0; ki<4; ki++)
	face_derivs_[ki] = std::min(1, face_fix[ki]);
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
      grid_start_[2] = grid_start[2];
      cell_size_[0] = cell_size[0];
      cell_size_[1] = cell_size[1];
      cell_size_[2] = cell_size[2];
    }

    /// Decide if ghost points should be constructed
    //void setMakeGhostPoints(bool make_ghost_points)
    //{
      //  make_ghost_points_ = make_ghost_points;
    //}

    /// Decide if the MBA algorithm should be used to approximate
    /// points (default == false)
    void setUseMBA(bool useMBA)
    {
      useMBA_ = useMBA;
    }

    void setInitMBA(bool initMBA)
    {
      initMBA_ = initMBA;
    }

    /// Add lower constraint. Only functional (1D volume)
    void addLowerConstraint(double minval)
    {
      has_min_constraint_ = true;
      minval_ = minval;
    }

    /// Add upper constraint. Only functional (1D volume)
    void addUpperConstraint(double maxval)
    {
      has_max_constraint_ = true;
      maxval_ = maxval;
    }

    /// Set local constraint on coefficient. Only functional (1D volume)
    void setLocalConstraint(double constraint_factor)
    {
      has_local_constraint_ = true;
      constraint_fac_ = constraint_factor;
    }

    /// Set minimum element size
    void setMinimumElementSize(double usize_min, double vsize_min, double wsize_min)
    {
      usize_min_ = usize_min;
      vsize_min_ = vsize_min;
      wsize_min_ = wsize_min;
    }

    /// Iteration count on which to switch to MBA (if LS)
    //void setSwitchToMBA(int iter)
    //{
    //  toMBA_ = iter;
    //}

    /// Set if an initial MBA volume should be used in LS approximation
    /// Default is false
    //void setInitMBA(bool initMBA, double start_coef=0.0)
    //{
    // initMBA_ = initMBA;
    //  initMBA_coef_ = start_coef;
    //}

    /// Set information about variable tolerance (default not triggered)
    void setVarTolBox(std::vector<TolBox> tolerances)
    {
      tolerances_ = tolerances;
    }
    
    /// Local measure for when the fraction of points outside the tolerance should
    /// lead to volume splitting
    void setOutFraction(double outfrac)
    {
      outfrac_ = outfrac;
    }
    
   /// Whether or not intermediate information should be written to
    /// standard output (default is not)
    void setVerbose(bool verbose)
    {
      verbose_ = verbose;
    }

    /// When everything else is set, this function can be used to run the 
    /// approximation process and fetch the approximated volume.
    /// \retval maxdist report the maximum distance between the approximated 
    ///                 volume and the data points
    /// \retval avdist report the average distance between the approximated 
    ///                volume and the datapoints
    /// \param nmb_out_eps report the number of data points that were found to
    ///                    lie outside the geometric tolerance (as specified by 
    ///                    the 'aepsge' argument to the ApproxVol constructor.
    /// \param max_iter maximum number of iterations
    /// \return a shared pointer to the generated LRSplinVolume, approximating 
    ///         the points as specified.
    shared_ptr<LRSplineVolume> getApproxVol(double& maxdist,
                                            double& avdist_all,
                                            double& avdist,
                                            int& nmb_out_eps,
                                            int& max_iter);

     /// Fetch additional accuracy information relevant for variable
    /// tolerance. To be called after getApproxVol
    void fetchOutsideTolInfo(double& max_outside,
			     double& average_outside)
    {
      max_outside = maxout_;
      average_outside = avout_;
    }

    
     /// Feature output
    void setFeatureOut(int ncell1, int ncell2, int ncell3)
    {
      write_feature_ = true;
      ncell1_ = ncell1;
      ncell2_ = ncell2;
      ncell3_ = ncell3;
    }

    void unsetFeatureOut()
    {
      write_feature_ = false;
    }

    void setFeatureLevel(std::vector<int>& levels)
    {
      feature_levels_ = levels;
    }
    
private:
    shared_ptr<LRSplineVolume> vol_;
    int nmb_pts_;
    std::vector<double>& points_;  // Reference to input points and parameter values
    std::vector<int> coef_known_;
    shared_ptr<LRSplineVolume> prev_;  // Previous volume, no point information
    // in elements
    
    bool useMBA_;         // Only LR-MBA
    //int toMBA_;         // Start with LR-MBA at the given iteration step
    bool initMBA_;      // The initial volume is made using LR-MBA
    double initMBA_coef_; // Initial height of constant volume

    int face_derivs_[6];
    double maxdist_;
    double maxdist_prev_;
    double avdist_;
    double avdist_all_;
    double avdist_all_prev_;
    int outsideeps_;
    int outsideeps_prev_;
    double maxout_;
    double avout_;
    double aepsge_;
    double smoothweight_;
    //bool smoothbd_;
    bool repar_;
    bool check_close_;
    bool fix_corner_;
    bool to3D_;
    bool grid_;
    bool check_init_accuracy_;
    double grid_start_[3];
    double cell_size_[3];
    bool initial_volume_;

    bool has_min_constraint_;
    bool has_max_constraint_;
    double minval_;
    double maxval_;
    bool has_local_constraint_;
    double constraint_fac_;
    double outfrac_;
    bool verbose_;
    double usize_min_;  // Minimum element size in u direction, negative 
    // if not set
    double vsize_min_;  // Minimum element size in v direction, negative 
    // if not set
    double wsize_min_;  // Minimum element size in w direction, negative 
    // if not set
    
    bool fix_boundary_;
    //bool make_ghost_points_;

    // Variable tolerance 
    std::vector<TolBox> tolerances_;

    // Features output
    bool write_feature_;
    int ncell1_, ncell2_, ncell3_;
    std::vector<int> feature_levels_;
    
    // DEBUG
    int ref_x_, ref_y_, ref_z_;
    int nmb1_, nmb2_, nmb3_;
    
    /// Define free and fixed coefficients
    void setCoefKnown();
    void unsetCoefKnown();
    void updateCoefKnown();
    void setCoefKnown(Direction3D d, Direction3D d2, bool atstart, int fixcoef);

    void computeAccuracy(std::vector<Element3D*>& ghost_elems); // @obar: why ghost???
    // The same as the above, but with OpenMP support (if flag is turned on).
    void computeAccuracy_omp(std::vector<Element3D*>& ghost_elems);
    void computeAccuracyElement(std::vector<double>& points, 
                                int nmb, int del, const Element3D* elem);
    //// The same as the above, but with OpenMP support (if flag is turned on).
    void computeAccuracyElement_omp(std::vector<double>& points,
                                   int nmb, int del, const Element3D* elem);

    /// Refine volume
    int refineVol(double threshold);
    void refineVol2();
    
    /// Create initial LR B-spline volume
    void makeInitVol(int dim);
    void makeInitVol(shared_ptr<SplineVolume> vol);
    void makeInitVol(int dim, 
                     int ncoef_u, int order_u, 
                     int ncoef_v, int order_v,
                     int ncoef_w, int order_w,
                     double domain[6]);
    void makeInitVol(int dim, 
                     int ncoef_u, int order_u, 
                     int ncoef_v, int order_v,
                     int ncoef_w, int order_w );
    void makeInitVol(int dim, 
                     int ncoef_u, int order_u,
                     int ncoef_v, int order_v, 
                     int ncoef_w, int order_w, 
                     double *knots_u, double *knots_v, double *knots_w);

    /// Create spline volume
    shared_ptr<SplineVolume> createVol(double* points, int nmb_pts,
                                       int dim, 
                                       int ncoef_u, int order_u, 
                                       int ncoef_v, int order_v, 
                                       int ncoef_w, int order_w, 
                                       double *knots_u, double *knots_v,
				       double *knots_w);
    
    /// Parameter domain surrounding the parameter values of all data points
    void computeParDomain(int dim, 
                          double& umin, double& umax, 
                          double& vmin, double& vmax,
                          double& wmin, double& wmax );

    void defineRefs(LRBSpline3D* bspline, double average_out,
		    std::vector<LRSplineVolume::Refinement3D>& refs_x,
		    std::vector<LRSplineVolume::Refinement3D>& refs_y,
		    std::vector<LRSplineVolume::Refinement3D>& refs_z,
		    std::vector<Element3D*>& elem_div);

    void checkFeasibleRef(Element3D* elem, 
			  std::vector<LRSplineVolume::Refinement3D>& refs_x,
			  std::vector<LRSplineVolume::Refinement3D>& refs_y,
			  std::vector<LRSplineVolume::Refinement3D>& refs_z,
			  std::vector<Element3D*>& affected);

      void appendRef(std::vector<LRSplineVolume::Refinement3D>& refs,
		     LRSplineVolume::Refinement3D& curr_ref,
		     double tol);
    
    void adaptVolumeToConstraints();

    // Turn function into a 3D volume
    void turnTo3D();
};

}  // namespace Go

#endif
