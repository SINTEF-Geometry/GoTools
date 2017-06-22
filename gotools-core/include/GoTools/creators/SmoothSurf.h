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

#ifndef _SMOOTHSURF_H_
#define _SMOOTHSURF_H_


//   -----------------------------------------------------------------------
//      Interface file for class SmoothSurf
//   -----------------------------------------------------------------------
//
//       Used to modify a tensor product B-spline surface with respect
//       to conditions on smoothness, editing constraints and boundary
//       conditions.
//
//       Implementation of the member functions are given in the
//       following files:
//
//          1. SmoothSurf.C
//
//   -----------------------------------------------------------------------
//    Written by: Vibeke Skytt                            08-11.1993.
//    Revised by: Vibeke Skytt                               09.1995
//    Revised by: Vibeke Skytt                               04.1998
//   -----------------------------------------------------------------------

#include "GoTools/geometry/SplineSurface.h"
#include "GoTools/creators/ConstraintDefinitions.h"

#include <vector>


namespace Go
{

    /// This class modifies a tensor product B-spline surface with respect
    /// to conditions on smoothness, editing constraints and boundary conditions.
class GO_API SmoothSurf
{

public:
    /// Default constructor. Initializes class variable.
    SmoothSurf();

    /// Constructor. Initializes class variable.
    /// \param copy_coefs true if coefficients on attached surface are not to be modified.
    SmoothSurf(bool copy_coefs);

    /// Destructor.
    virtual
    ~SmoothSurf();

    /// Initializes data given by an intermediate surface.
    /// \param insf the initial surface.
    /// \param seem continuity across opposite edges.
    ///             0 not specified, 1 = C0, 2 = C1, 3 = C2.
    /// \param coef_known whether the coefs are free to be altered.
    ///                   0: not known, 1: known
    ///                   >= kpointer_ = 3: the coefficients indexed by ki &
    ///                                     coef_known[ki] - kpointer_ should be equal.
    /// \param num_side_constraints the number of linear side constraints in the system.
    /// \param has_normal_cond whether the system must fulfill normal conditions.
    void attach(shared_ptr<SplineSurface>& insf, int seem[], int coef_known[],
		int num_side_constraints = 0, int has_normal_cond = 0);

    /// Compute the smoothing part of the equation system.
    /// \param weight1 contribution weight with respect to the 1st derivative.
    /// \param weight2 contribution weight with respect to the 2nd derivative.
    /// \param weight3 contribution weight with respect to the 3rd derivative.
    virtual
    void setOptimize(const double weight1, const double weight2,
		     const double weight3);

    /// Compute matrices for least squares approximation.
    /// \param pnts points on surface to be approximated.
    ///             Stored as (for 3D): (x0, y0, z0, x1, y1, z1, ...)
    /// \param param_pnts the corresponding 2-dimensional parametric points.
    ///                   Stored as: (u0, v0, u1, v1, ...)
    /// \param pnt_weights each of the pnts is assigned a weight lying in
    ///                    the unit interval. 1.0 if all pnts are equally important.
    /// \param weight the contribution of the approximation of the pnts in the system.
    ///               weight should lie in the unit interval.
    void setLeastSquares(const std::vector<double>& pnts,
			 const std::vector<double>& param_pnts,
			 const std::vector<double>& pnt_weights,
			 const double weight);

    /// Compute matrices for approximation of normal directions.			 
    /// \param pnts normals in sf to be approximated.
    ///             Stored as (for 3D): (x0, y0, z0, x1, y1, z1, ...)
    /// \param param_pnts the corresponding 2-dimensional parametric points.
    ///                   Stored as: (u0, v0, u1, v1, ...)
    /// \param pnt_weights each of the pnts is assigned a weight lying in
    ///                    the unit interval. 1.0 if all pnts are equally important.
    /// \param weight the contribution of the approximation of the normals in the system.
    ///               weight should lie in the unit interval.
    /// \return status value: 0 = OK, 1 = warning (system not prepared for normal conditions).
    int setNormalCond(const std::vector<double>& pnts,
		      const  std::vector<double>& param_pnts,
		      const std::vector<double>&  pnt_weights,
		      const double weight);

    /// Compute the contribution to the equation system from the approximation
    /// of an original surface, i.e. the contribution of the original coefficients.
    /// \param weight the relative contribution of the original coefs.
    ///               Should lie in the unit interval.
    void approxOrig(double weight);

    /// Set periodicity constraints in one par. dir.
    /// \param pardir the direction of the periodicity, 1 == udir && 2 == vdir.
    /// \param cont the continuity across the seem, 0 = C0, 1 = C1, 2 = C2.
    /// \param weight1 C1 continuity contribution used in approximation.
    /// \param weight2 C2 continuity contribution used in approximation.
    virtual
    void setPeriodicity(int pardir,
		       int cont,
		       double weight1,
		       double weight2);


    /// Add linear side constraints to the system equation.
    /// \param constraints the linear side constraints between surface coefficients.
    void setSideConstraints(std::vector<sideConstraint>& constraints);

    /// Solve equation system, and produce output surface.
    /// If failing to solve the routine may throw an exception.
    /// \param surf the output surface.
    /// \return 0 = OK, negative = failed solving system.
    int equationSolve(shared_ptr<SplineSurface>& surf);

    /// Set the relaxation parameter for the RILU preconditioner.
    void setRelaxParam(double omega);

protected:
    int norm_dim_;         // If the problem has normal-conditions: 3, otherwise: 1
    int idim_;             // Dimension of geomtry space. 
    int kdim_;             // Dimension of homogeneous space.
    int idim1_;            // Dimension of projective space (for rationals).
    int ider_;             // Maximum derivative involved in the computations. 
    bool integralset_; //
    int ider_scratch_;     // Number of derivatives used in scratch allocation

    // Rational case
    bool rational_;    // Is surface rational?
    shared_ptr<SplineSurface> bspline_surface_; // 1-dimensional rational surface with
                                                       // 0's as control values, and weights from 
                                                       // input surface. Used to calculate B-spline
                                                       // products for the rational case

    int cont_seem_[2];     // Number of rows affected by continuity at a seem
    // for each parameter direction.
    int kncond_;           // Size of matrix system (# unknown coefficents + # side constraints)
    int knconstraint_;     // Number of side constraints.
    const int kpointer_; // If coefknown_[ki] >= kpointer_, the coefficients indexed by
                   // ki & coefknown_[ki] - kpointer_ should be equal.

    int *coefknown_; // Array indicating the status of coefficients, i.e.
                     // free, fixed, not involved, equal to a given coefficients.
                     // 0: not known, 1: known, 2: not of interest (i.e. assumed known).
                     // >= kpointer_: the coefficients indexed by
                     //               ki & coefknown_[ki] - kpointer_ should be equal.
    std::vector<int> pivot_;    // Array giving the position of the free coefficients
    // in the equation system.

    // Parameters defining the spline space.
    shared_ptr<SplineSurface> srf_;        // Pointer to input surface.
    int kk1_, kk2_;         // Order of surface in both parameter directions.    
    int kn1_, kn2_;         // Number of coefficients of surface.                
    std::vector<double>::const_iterator st1_;    // Pointer to knot vector of 
    // surface in 1. par. dir. 
    std::vector<double>::const_iterator st2_;    // Pointer to knot vector of 
    // surface in 2. par. dir. 

    const bool copy_coefs_;

    /// Parameters used to define the specific input spline surface.
    std::vector<double> coef_array_; // Only used if copy_coefs_ == true
    std::vector<double>::iterator scoef_;   // Pointer to surface coefficients.      

    /// Storage of the equation system.
    std::vector<double> gmat_;         // Matrix at left side of equation system.  
    std::vector<double> gright_;       // Right side of equation system.      

    ///   Free all memory allocated for class members.
    virtual
    void releaseScratch(); 

    /// Prepare storage for integrals of inner products of basis functions.
    virtual
    void prepareIntegral();

    /// Given the value of non-zero B-spline functions, compute the value
    /// of the corresponding surface basis function (i.e. the products of
    /// the u- and v-basis functions).
    /// \param sb1 the basis values in the first parameter direction (u).
    /// \param sb2 the basis values in the second parameter direction (v).
    /// \param kleft1 index of the first knot interval in u-dir.
    /// \param kleft2 index of the first knot interval in v-dir.
    /// \param ider the number of derivatives to compute.
    /// \param sbasis the computed basis values in sf.
    ///               size = order_u()*order_v()*(ider + 1)*(ider + 1).
    ///               The space must be allocated on the outside.
    virtual
    void getBasis(const double *sb1, const double *sb2, int kleft1, int kleft2,
		  int ider, double *sbasis);

    /// Set pointers between identical coefficients at a periodic seem.
    /// If possible, update fixed coefficients at the seem.
    /// \param seem continuity across the seem. Array size = 2.
    void preparePeriodicity(int seem[]);

    /// Compute the smoothing part of the equation system, non-rational case.
    /// \param weight1 contribution weight with respect to the 1st derivative.
    /// \param weight2 contribution weight with respect to the 2nd derivative.
    /// \param weight3 contribution weight with respect to the 3rd derivative.
    void setOptimizeNonrational(const double weight1, const double weight2,
				const double weight3);

    /// Compute the smoothing part of the equation system, rational case.
    /// \param weight1 contribution weight with respect to the 1st derivative.
    /// \param weight2 contribution weight with respect to the 2nd derivative.
    /// \param weight3 contribution weight with respect to the 3rd derivative.
    void setOptimizeRational(const double weight1, const double weight2,
			     const double weight3);

    /// Compute the contribution to the equation system from the approximation
    /// of an original surface, i.e. the contribution of the original coefficients
    /// for the non-rational case
    /// \param weight the relative contribution of the original coefs.
    ///               Should lie in the unit interval.
    void approxOrigNonrational(double weight);

    /// Compute the contribution to the equation system from the approximation
    /// of an original surface, i.e. the contribution of the original coefficients
    /// for the rational case
    /// \param weight the relative contribution of the original coefs.
    ///               Should lie in the unit interval.
    void approxOrigRational(double weight);

    /// Set constraints on C1-continuity across a seem and add these to the
    /// equation system.
    void setC1AtSeem(int pardir, double weight);

    /// Set constraints on C2-continuity across a seem and add these to the
    /// equation system.
    void setC2AtSeem(int pardir, double weight);

    /// Set constraints on C1- and C2-continuity across a seem and add these to the
    /// equation system when surface is rational.
    void setRationalCnAtSeem(int pardir, int cn, double weight1, double weight2);

    /// Ensure that the expected C1 or C2 continuity at the seem is satisfied.
    /// Only applicable if seem_[0] > 1 || seem_[1] > 1.
    virtual
    int adjustAtSeem();


private:

    // Parameters used in integration
    std::vector<double> vec1_, vec2_;
    double ***integral1_;  // Array used to store integrals of inner product
    // of derivatives of B-splines in 1. par. dir.  
    double ***integral2_;  // Array used to store integrals of inner product
    // of derivatives of B-splines in 2. par. dir.  

    double omega_;

}; // end of class SmoothSurf


} // end of namespace Go

#endif
