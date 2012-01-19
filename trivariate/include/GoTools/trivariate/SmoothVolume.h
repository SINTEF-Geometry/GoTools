//===========================================================================
//
// File : SmoothVolume.h
//
// Created: Tue Dec 15 10:10:53 2009
//
// Author: Kjell Fredrik Pettersen
//
// Revision: $Id:$
//
// Description:
//
//===========================================================================



#ifndef __SMOOTHVOLUME_H
#define __SMOOTHVOLUME_H



#include "GoTools/trivariate/SplineVolume.h"

#include <memory>
#include <vector>



namespace Go
{


  enum CoefStatus
  {
    CoefFree,      // Coefficient is not known, and free to be determined by the equation system
    CoefKnown,     // Coefficient is not to be altered
    CoefAvoid,     // Coefficient is not of interrest (i.e. assumed known)
    CoefOther      // Coefficient is identical to another, the other should be determined first
  };

  /// \brief This class modifies a NURBS volume with respect to smoothness, editing
  /// constraints and boundary conditions. Specified coefficients are adjusted
  /// to minimize a functional combining the different modification conditions.

  class SmoothVolume
  {

  public:

    /// Default constructor. Initializes class variable.
    SmoothVolume();

    /// Constructor. Initializes class variable.
    /// \param copy_coefs true if coefficients on attached volume are not to be modified.
    SmoothVolume(bool copy_coefs);

    // Destructor
    virtual ~SmoothVolume();


    /// Initializes data given by an intermediate volume.
    /// \param in_vol the initial volume.
    /// \param status array with the freedom for each coefficient of the volume
    ///                   CoefFree: Not known,
    ///                   CoefKnown: Known,
    ///                   CoefAvoid: Not of interrest
    void attach(shared_ptr<SplineVolume>& in_vol,
		std::vector<CoefStatus> status);


    /// Reset all smoothing parameters, but keep the input volume
    void reset();


    /// Add weights for the smoothing part of the equation system.
    /// \param weight1 contribution weight with respect to the 1st derivative.
    /// \param weight2 contribution weight with respect to the 2nd derivative.
    /// \param weight3 contribution weight with respect to the 3rd derivative.
    void setOptimize(const double weight1, const double weight2,
		     const double weight3);


    /// Add data for the least squares approximation.
    /// \param pnts points in volume to be approximated.
    ///             Stored as (for 3D): (x0, y0, z0, x1, y1, z1, ...)
    /// \param param_pnts the corresponding 3-dimensional parametric points.
    ///                   Stored as: (u0, v0, w0, u1, v1, w1, ...)
    /// \param pnt_weights each of the pnts is assigned a weight lying in
    ///                    the unit interval. 1.0 if all pnts are equally important.
    /// \param weight the contribution of the approximation of the pnts in the system.
    ///               weight should lie in the unit interval.
    void setLeastSquares(std::vector<double>& pnts,
			 std::vector<double>& param_pnts,
			 std::vector<double>& pnt_weights,
			 const double weight);

    /// Add periodicity constraints in one par. dir.
    /// \param pardir the direction of the periodicity, 0, 1 or 2
    /// \param cont the continuity across the seem, 0 = C0, 1 = C1, 2 = C2.
    /// \param weight1 C1 continuity contribution used in approximation.
    /// \param weight2 C2 continuity contribution used in approximation.
    void setPeriodicity(int pardir,
		       int cont,
		       double weight1,
		       double weight2);

    /// Set up and solve equation system, and produce output volume.
    /// If failing to solve the routine may throw an exception.
    /// \param vol the output surface.
    /// \return 0 = OK, negative = failed solving system.
    int equationSolve(shared_ptr<SplineVolume>& vol);

  private:

    shared_ptr<SplineVolume> input_volume_;    // Input volume
    shared_ptr<SplineVolume> bspline_volume_;  // 1-dimension rational volume function with
                                                      // 0's as control values and weights from input
                                                      // volume. Used to calculate B-spline products
                                                      // for rational case.

    std::vector<CoefStatus> coef_status_;     // The status of each coefficient in the equation system
    std::vector<int> coef_other_;   // Coefficient i and coef_other_[i] should be identical, the second
                                    // to be determined first. Only used for coef_status_[i] == CoefOther

    // Integrals of product of two B-spline functions of same derivative order.
    // The value of D^k(B_i(u)) * D^k(B_j(u)) integrated over all u is stored in
    //
    //     bsplineintegral_u_[ (j-i+o-1) + i*(2*o-1) + k*n*(2*o-1)]
    //
    // where o and n are the order and number of coefficients in the B-spline space in first parameter
    // direction. bsplineintegral_v_ and bsplineintegral_w_ are used the same way for 2. and 3. par.dir.
    std::vector<double> bsplineintegral_u_, bsplineintegral_v_, bsplineintegral_w_;

    // Integrals of product of two B-spline functions of derivative order differing by 2.
    // The value of D^(k+2)(B_i(u)) * D^k(B_j(u)) integrated over all u is stored in
    //
    //     bsplineintegral_skew_u_[ (j-i+o-1) + i*(2*o-1) + k*n*(2*o-1)]
    std::vector<double> bsplineintegral_skew_u_, bsplineintegral_skew_v_, bsplineintegral_skew_w_;

    bool copy_coefs_;    // If true, copy coefficients into coef_array_, and leave the coefficients
                         // of the input volume untouched when solving. If false, store the new
                         // coefficients into the input volume after solving the linear system
    std::vector<double> coef_array_;   // Holds the coefficients only in the case copy_coefs_ == true
    std::vector<double>::iterator it_coefs_;   // Iterator over volume coefficients

    double weight_der_1_;  // Contribution weight when optimizing 1st derivative
    double weight_der_2_;  // Contribution weight when optimizing 2st derivative
    double weight_der_3_;  // Contribution weight when optimizing 3st derivative

    double weight_least_sq_;   // Contribution weight for least squares approximation.
    std::vector<double> least_sq_pts_;   // Point coordinates used for least squares approxiamtion (x0, y0, z0, x1 etc
    std::vector<double> least_sq_params_;   // Parameters of each point during least squares approxiamtion (u0, v0, w0, u1, etc)
    std::vector<double> least_sq_wgt_;   // Weights of each point during least squares approxiamtion

    int seem_cont_[3];     // Continuity at seem for each parameter direction. -1 = no periodicity, 0 = C0, 1 = C1, 2 = C2
    double seem_weight_[6];   // Weights for C1 and C2 continuity at seem. In order C1 at u-dir, C2 at u-dir, C1 at v-dir, etc

    /// Storage of the equation system.
    int nmb_free_;     // Number of free variables in equation system
    std::vector<double> gmat_;       // Matrix at left side of equation system
    std::vector<double> gright_;     // Right side of equation system
    std::vector<int> pivot_;         // Array giving the position of the free coefficients



    // Calculated "skew product" integrals, from integrals of splines of same derivative order,
    // using the chain rule for integration
    void build_skew_integral(int new_depth, int old_depth,
			     const BsplineBasis& basis,
			     const std::vector<double>& integral,
			     std::vector<double>& skew_integral);

    // Set relations to get periodic C0-continuity when desired. Must be done before building the equation matrices
    void setPeriodicityConstraints();

    // Set relations between four edges that are supposed to coincide during double periodicity.
    // Input is the first four coefficients of the edges, the number of coefficient quadruples to
    // coincide (= the number of coefficients along each edge), the incrementation at each iteration
    // to get the indices of the next quadruple, and a boolean telling if the start- and end coefficients
    // are to coincide as a consequence of a third periodicity
    void setEdgePeriodicity(int corner0, int corner1, int corner2, int corner3,
			    int nmb, int step, bool per_end);

    // Set relation between two opposite edges to coincide because of periodicity in one direction.
    // Input is positions of the first pair of coefficients to coincide, the number of coefficient pairs
    // in the two parameter directions, the step in each direction to get from the coefficient indices
    // of one pair to the next, and wether the volume has an additional periodicity in avy of the two directions
    void setFacePeriodicity(int bottom, int top, int nmb_0, int nmb_1,
			    int step_0, int step_1, bool per_0, bool per_1);

    // For a sequence of sets of control points, set relations between them so that they coincide.
    // Input is the first set, the length of the sequence, and the incrementation at each set to get
    // the next set.
    void setPeriodicityLocally(const std::vector<int>& coefs_in, int nmb, int step);

    void resetPivotAndMatrices();

    // Extend (or build for the first time) the integrals of products of B-spline functions.
    // Only used for the non-rational case. For rational cases, our integrals will be on
    // a function with a denominator, then we can not split into one-dimensional integrals
    // over B-spline products.
    void buildIntegrals();

    void buildBsplineVolume();

    // Add contributions to equation system for least squares approximation
    void addLeastSquares();

    // Add contribution to equation system for smoothness optimizations, non-rational case
    void addOptimizeNonrational();

    // Add contribution to equation system for smoothness optimizations, rational case
    void addOptimizeRational();

    // Add contribution to equation system for periodicity optimizations, nonrational case
    void addNonrationalContinuityAtSeem(int pardir);

    // Add contribution to equation system for periodicity optimizations, rational case
    void addRationalContinuityAtSeem(int pardir);

    // Solve system
    int solve();

    int geoDim() const;     // Dimension of geometry space

    int homogDim() const;   // Dimension of homogeneous space

    bool rational() const;  // Is volume rational?

    int numCoefs() const;             // Product of numbers of coefficients in all parameter directions
    int numCoefs(int pardir) const;   // Number of coefficients in parameter direction (pardir is 0, 1 or 2)

    int order(int pardir) const;   // Polynomial order in parameter direction (pardir is 0, 1 or 2)

    BsplineBasis basis(int pardir) const;


  };    // Class SmoothVolume


} // namespace Go

#endif    // #ifndef __SMOOTHVOLUME_H
